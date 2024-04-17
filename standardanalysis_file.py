# pampro-manager-processing-scripts - processing scripts for the pampro pipeline managed by 'pampro-manager'
# Copyright (C) 2019  MRC Epidemiology Unit, University of Cambridge
#   
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
#   
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#   
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Script to calibrate data and produce statistical results
## pampro processing pipeline 1

from datetime import datetime, timedelta
import sys
import os
import json
from pampro import data_loading, channel_inference, batch_processing, triaxial_calibration, time_utilities, pampro_utilities, pampro_fourier, Bout, Channel
import pandas as pd
import numpy as np
import traceback
import re
from  multiprocessing import Pool, cpu_count
import argparse
from functools import partial

#######################################################################################################################

#def standardanalysis(job_details, settings):
def standardanalysis(setting, filename=None, job_details=None):
    # number of iterations to be used when optimising during calibration
    num_iterations = 500
    
    # threshold of end error for data validity
    end_error_threshold = 10
    
    results_folder = settings.get("results_folder")[0]
    plots_folder = settings.get("plots_folder")[0]
    noise_cutoff_mg = settings.get("noise_cutoff_mg")[0]
    processing_epoch = settings.get("processing_epoch")[0]
    epochs = settings.get("epochs")[0]
    collapse = settings.get("collapse_data")[0].strip("()'',")
    whole_file = settings.get("whole_file")[0].strip("()'',")
    temperature_calibration = settings.get("temperature_calibration")[0].strip("()'',")
    
    if job_details is not None:
      # extract the multi-file calibration errors from the job details
      mf_start_error = job_details["start_error"]
      mf_end_error = job_details["end_error"]
    else:
      mf_start_error = None
      mf_end_error = None

    # generate the epochs dictionary
    epoch_dict = pampro_utilities.json_epochs_to_dict(epochs)

    # generate the list of cut points
    intensities = settings.get("cutpoints")[0]
    intensities_list = pampro_utilities.json_cutpoints_to_list(intensities)

    # generate the list of angles
    angles = settings.get("angles")[0]
    angles_list = pampro_utilities.json_cutpoints_to_list(angles)

    # generate a list of stats that have been selected (e.g. have value == 1)
    stats_list = ['integrity']
    for stat_name in ['enmo', 'hpfvm', 'pitch', 'roll', 'temperature', 'battery']:
        if settings.get(stat_name)[0] == '1':
            stats_list.append(stat_name)

    stats = pampro_utilities.define_statistics(stats_list, intensities_list, angles_list)

    # find which time increments (epochs), if any, are to be plotted
    epochs_plot = settings.get("epochs_plot")[0]
    
    # create a list of increments to be used later
    plot_inc = []
    
    # if epochs_plot is empty, leave plotting list empty
    if epochs_plot is np.nan:
        pass

    else:
        incs = json.loads(epochs_plot)
        for i in incs:
            name = i["name"]
            if i["plot"] == 1:
                plot_inc.append(name)
    
    # find which statistics are to be plotted, and create a dictionary to be used later
    plotting_dict = dict()
    for var, plot_name in zip(["enmo_plot", "hpfvm_plot", "pitch_plot", "roll_plot", "temperature_plot", "battery_plot"],["ENMO_sum", "HPFVM_sum", "PITCH_mean", "ROLL_mean", "Temperature_mean", "Battery_mean"]):
        # only add to plotting dict if selected for processing, i.e. is in stats_list
        if settings.get(var)[0] == '1' and (var.split("_")[0]) in stats_list:
            plotting_dict[plot_name] = "{}_{}_" + plot_name + ".png"
   
    if job_details is not None: 
        # get job details
        pid = str(job_details["pid"])
        filename = str(job_details["filename"])
        hdf5_filename = str(job_details["hdf5_filename"])
    else:
        hdf5_filename = filename
        pid = '000'
    

    filename_short = os.path.basename(filename).split('.')[0]

    analysis_meta = os.path.join(results_folder, "analysis_meta_{}.csv".format(filename_short))
    # check if analysis_meta already exists...
    if os.path.isfile(analysis_meta):
        os.remove(analysis_meta)
        
    # Load the data from the hdf5 file
    ts, header = data_loading.load(hdf5_filename, "HDF5", hdf5_mode="r+", hdf5_group=("Resampled"))

    # some monitors have manufacturers parameters applied to them, let's preserve these but rename them:
    var_list = ["x_gain", "x_offset", "y_gain", "y_offset", "z_gain", "z_offset", "calibration_date"] 
    for var in var_list:
        if var in header.keys():
            header[("manufacturers_%s" % var)] = header[var]
            header.pop(var)

    x, y, z, temperature, battery, integrity = ts.get_channels(["X", "Y", "Z", "Temperature", "Battery", "Integrity"])

    # Calculate the start and end points to use, depending on whether the whole file is required
    if whole_file == "YES":
        start = time_utilities.start_of_day(x.timeframe[0])
        end = time_utilities.end_of_day(x.timeframe[-1])
        # calculate the time period to use
        tp = (start, end)
        num_days = (end - start).days

    else:
        num_days = settings.get("days_of_data")[0].strip("(),")
        start = time_utilities.start_of_day(x.timeframe[0])
        end = start + timedelta(days=int(num_days))
        # calculate the time period to use
        tp = (start, end)

    # check the integrity channel for "1", which shows that anomalies have been fixed BEFORE applying calibration factors
    if 1 in integrity.data:
        # extract the bouts of the integrity channel where the data == 1 (the "integrity compromised" value)
        exclusion_bouts = integrity.bouts(1, 1)
    else:
        exclusion_bouts = []
    
    epoch_list = []
    for e in epoch_dict.keys():
        epoch_list.append(e)
    
    ################ CALIBRATION #######################
    
    ## FIRST CALIBRATE FILE INDIVIDUALLY
    
    # extracting the still bouts from the data
    if temperature_calibration == "YES":
        calibration_ts, calibration_header = triaxial_calibration.calibrate_stepone(x, y, z, temperature, noise_cutoff_mg=noise_cutoff_mg)
    else:
        calibration_ts, calibration_header = triaxial_calibration.calibrate_stepone(x, y, z, noise_cutoff_mg=noise_cutoff_mg)
    
    # Calibrate the acceleration to local gravity
    file_cal_diagnostics = triaxial_calibration.calibrate_steptwo(calibration_ts, calibration_header, calibration_statistics=False, num_iterations=num_iterations)
    
    # extract the file-level calibration errors from the calibration_results
    file_start_error = file_cal_diagnostics["start_error"]
    file_end_error = file_cal_diagnostics["end_error"]
    
    ## SECOND EXAMINE FILE END ERROR 
    
    # compare individual file error to error threshold
    if file_end_error < end_error_threshold:
        calibration_type = "single"
        
        # calibrate the acceleration to local gravity using calibration factors derived from the file itself:
        if temperature_calibration == "YES":
            try:
                triaxial_calibration.do_calibration(x, y, z, temperature, file_cal_diagnostics)
            # catch if ValueError occurs, if it's to do with broadcasting arrays, calibrate without using temperature
            except ValueError as Argument:
                assert (str(Argument).startswith("operands could not be broadcast together")), "different type of ValueError; not a broadcasting error"
                triaxial_calibration.do_calibration(x, y, z, temperature=None, cp=file_cal_diagnostics)
                header["calibration_temperature_fail"] = "True"
        else:
            triaxial_calibration.do_calibration(x, y, z, temperature=None, cp=file_cal_diagnostics)
    
        metadata = {**header, **file_cal_diagnostics}
        for key in ["start_error", "end_error"]:
            metadata[("file_%s" % key)] = metadata[key]
            metadata.pop(key)
    
    else:
      
      if job_details is not None:

        # compare multi-file error to error threshold
        if float(mf_end_error) < end_error_threshold:
            calibration_type = "multi"
        
            # Calibrate the acceleration to local gravity using multi file calibration factors from the calibration database
            calibration_dict = {"x_scale": float(job_details["x_scale"]),
                                "x_offset": float(job_details["x_offset"]),
                                "x_temp_offset": float(job_details["x_temp_offset"]),
                                "y_scale": float(job_details["y_scale"]),
                                "y_offset": float(job_details["y_offset"]),
                                "y_temp_offset": float(job_details["y_temp_offset"]),
                                "z_scale": float(job_details["z_scale"]),
                                "z_offset": float(job_details["z_offset"]),
                                "z_temp_offset": float(job_details["z_temp_offset"])                        
                                }
            
            try:
                triaxial_calibration.do_calibration(x, y, z, temperature, calibration_dict)
            # catch if ValueError occurs, if it's to do with broadcasting arrays, calibrate without using temperature
            except ValueError as Argument:
                assert (str(Argument).startswith("operands could not be broadcast together")), "different type of ValueError; not a broadcasting error"
                triaxial_calibration.do_calibration(x, y, z, temperature=None, cp=calibration_dict)
                header["calibration_temperature_fail"] = "True"

            metadata = {**header, **calibration_dict}
            metadata["calibration_method"] = job_details["calibration_method"]
            cal_files = job_details["files_used"]
            metadata["files_used_in_calibration"] = cal_files
            metadata["number_files_used"] = cal_files.count(".hdf5")
            metadata["calibration_date"] = job_details["calibration_date"]
            
        else:
            # neither calibration type has been successful
            metadata = {**header}
            calibration_type = "fail"
    
    metadata["hdf5_filename"] = hdf5_filename
    metadata["days_of_data_processed"] = num_days
    metadata["processing_epoch"] = processing_epoch
    metadata["analysis_resolutions"] = epoch_list
    metadata["analysis_statistics"] = stats_list
    metadata["noise_cutoff"] = noise_cutoff_mg
    metadata["mf_start_error"] = mf_start_error
    metadata["mf_end_error"] = mf_end_error
    metadata["file_start_error"] = file_start_error
    metadata["file_end_error"] = file_end_error
    metadata["calibration_type"] = calibration_type

    if job_details is not None:
       metadata["device"] = job_details["monitor"]
    else:
       metadata["device"] = re.findall(r"_(\d{5})_", filename)[0] 

    # Writing out the metadata to a file
    pampro_utilities.dict_write(analysis_meta, pid, metadata)

    ## THIRD ONLY PROCESS IF CALIBRATION TYPE IS NOT "FAIL"
    
    if calibration_type == "fail":
        results_files = None
        charts = None
        
        #os.system("chgrp {} {} & chmod 770 {}".format(group, analysis_meta, analysis_meta))
          
    else:
        results_files = [os.path.join(results_folder, "{}_{}.csv".format(name, filename_short)) for name in epoch_dict.keys()]
        files = [open(file, "w") for file in results_files]

        # Write the column headers to the created files
        for f in files:
            f.write(pampro_utilities.design_file_header(stats) + "\n")
        
        # delete the exclusion bouts, if any:
        x.delete_windows(exclusion_bouts)
        y.delete_windows(exclusion_bouts)
        z.delete_windows(exclusion_bouts)
        temperature.delete_windows(exclusion_bouts)
        battery.delete_windows(exclusion_bouts)

        # Derive some signal features
        vm = channel_inference.infer_vector_magnitude(x, y, z)
        # delete the exclusion bouts from vm too
        vm.delete_windows(exclusion_bouts)

        if 'hpfvm' in stats_list:
            vm_hpf = channel_inference.infer_vm_hpf(vm)
        else:
            vm_hpf = None

        if 'enmo' in stats_list:
            enmo = channel_inference.infer_enmo(vm)
        else:
            enmo = None

        if 'pitch' in stats_list or 'roll' in stats_list:
            pitch, roll = channel_inference.infer_pitch_roll(x, y, z)
        else:
            pitch = roll = None

        # Infer nonwear and mask those data points in the signal
        nonwear_bouts = channel_inference.infer_nonwear_triaxial(x, y, z, noise_cutoff_mg=noise_cutoff_mg)

        annotation_bouts = []
        for bout in nonwear_bouts:
            # Show non-wear bouts in orange
            bout.draw_properties = {'lw': 0, 'alpha': 0.5, 'facecolor': '#ffb366'}
            annotation_bouts.append(bout)

        if collapse == 'YES':
            for channel, channel_name in zip([enmo, vm_hpf, pitch, roll, temperature, battery], ["ENMO", "HPFVM", "PITCH", "ROLL", "Temperature", "Battery"]):
                if channel_name in stats:
                    # Collapse the sample data to a processing epoch (in seconds) so data is summarised
                    epoch_level_channel = channel.piecewise_statistics(timedelta(seconds=int(processing_epoch)), time_period=tp)[0]
                    epoch_level_channel.name = channel_name
                    if channel_name in ["Temperature", "Battery"]:
                        pass
                    else:
                        epoch_level_channel.delete_windows(nonwear_bouts)
                    epoch_level_channel.delete_windows(exclusion_bouts)
                    ts.add_channel(epoch_level_channel)
            
            # collapse binary integrity channel
            epoch_level_channel = integrity.piecewise_statistics(timedelta(seconds=int(processing_epoch)), statistics=[("binary", ["flag"])], time_period=tp)[0]
            epoch_level_channel.name = "Integrity"
            ts.add_channel(epoch_level_channel)

        else:
            # If do not want to collapse data to epoch level
            for channel, channel_name in zip([enmo, vm_hpf, pitch, roll, temperature, battery],
                                             ["ENMO", "HPFVM", "PITCH", "ROLL", "Temperature", "Battery"]):
                if channel_name in stats:
                    channel.name = channel_name
                    if channel_name in ["Temperature", "Battery"]:
                        pass
                    else:
                        channel.delete_windows(nonwear_bouts)
                    channel.delete_windows(exclusion_bouts)    
                    ts.add_channel(channel)

        # Piecewise output
        charts = []
        # loop through the epochs (time resolutions) for results
        for epoch, name, file in zip(epoch_dict.values(), epoch_dict.keys(), files):

            results_ts = ts.piecewise_statistics(epoch, statistics=stats, time_period=tp, name=pid)
            results_ts.write_channels_to_file(file_target=file)
            file.flush()

            # if epoch is to be plotted:
            if name in plot_inc:
            # for each statistic in the plotting dictionary, produce a plot in the charts folder
                for stat, plot in plotting_dict.items():
                    results_ts[stat].add_annotations(annotation_bouts)
                    results_ts[stat].add_annotations(exclusion_bouts)
                    chart_file = os.path.join(plots_folder, plot.format(filename_short, name))
                    results_ts.draw([[stat]], file_target=chart_file)
                    charts.append(chart_file)

        all_files = charts + results_files
        
        # change group and permissions of files
        #for f in all_files:
        #    os.system("chgrp {} {} & chmod 770 {}".format(group, f, f))

    for c in ts:
        del c.data
        del c.timestamps
        del c.indices
        del c.cached_indices
    
    
    
    return {"results": results_files, "analysis_meta_file": analysis_meta, "visualisation_file": charts}


#######################################################################################################################

def file_task(analysis_function, settings, filename, **kwargs):

    task = analysis_function.__name__
    #submission_id = settings.get("submission_id")[0]
    logs_folder = settings.get("logs_folder")[0]
    # archive = os.path.join(logs_folder, "archive")
    output_string = "_completed"
    error_string = "_unsuccessful"

    assert os.path.isfile(filename), f'File {filename} does not exist'
    head, tail = os.path.split(filename)
    job_name = tail.split('.')[0]
    job_start_time = datetime.now()

    try:
        output_dict = analysis_function(settings, filename, **kwargs)
        job_end_time = datetime.now()
        job_duration = job_end_time - job_start_time
        output_dict["job_duration"] = str(job_duration)

        output_log = logs_folder + os.sep + job_name + "_" + task + output_string + ".csv"
        #dict_write(output_log, pid, output_dict)
        pampro_utilities.dict_write(output_log, None, output_dict)
    except Exception:

        tb = traceback.format_exc()


        # Create the error file only if an error has occurred
        #with open(logs_folder + os.sep + job_name + "_" + task + error_string + submission_id + ".csv", "w") as error_log:
        with open(logs_folder + os.sep + job_name + "_" + task + error_string + ".csv", "w") as error_log:

            error_log.write("Error log at " + str(datetime.now()) + "\n")
            #for k, v in job.iteritems():
            #    error_log.write(str(k) + ": " + str(v) + "\n")
            error_log.write("Exception:" + str(sys.exc_info()) + "\n")
            error_log.write(tb + "\n\n")
            error_log.flush()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--filenames', dest='filenames', nargs="*", type=str,required=True)
    parser.add_argument('--settings_file', dest='settings_file', type = str, required=True)
    
    args = parser.parse_args() 
    settings_file = str(args.settings_file)
    settings = pd.read_csv(settings_file, dtype=str)

    filenames = args.filenames
    num_files = len(filenames)
    num_cores = cpu_count()

    assert num_files > 0, "the file list is empty !"
    assert num_cores >= num_files, "the number of files should be lower than the number of cores !"
   
    if num_files == 1:
        filename = str(filenames[0])
        file_task(standardanalysis, settings, filename)
    elif num_cores > num_files:
      filenames = [str(f) for f in filenames]
      with Pool(processes=num_files) as pool:
        pool.map(partial(file_task, standardanalysis, settings), filenames)

