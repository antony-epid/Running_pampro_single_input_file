# pampro-manager-processing-scripts - processing scripts for the pampro pipeline managed by 'pampro-manager'
# Copyright (C) 2019  MRC Epidemiology Unit, University of Cambridge
#   
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
#   
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#   
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Script to use still bouts from multiple files to derive calibration factors
## pampro processing pipeline 1

from datetime import datetime
import pandas as pd
import sys
import os
from pampro import data_loading, batch_processing, triaxial_calibration, pampro_utilities
from collections import OrderedDict

today = datetime.now().strftime('%d%b%Y')

settings_file = str(sys.argv[1])
jobs_file = str(sys.argv[2])
group = str(sys.argv[3])
job_num = int(sys.argv[4])
num_jobs = int(sys.argv[5])

#######################################################################################################################


def calibratemonitor(job_details, settings):

    results_folder = settings.get("results_folder")[0]

    monitor = job_details["monitor_id"]

    files = (job_details["stillbouts_files"]).strip("[]").replace("'", "").replace(" ", "").split(",")

    calibration_output = os.path.join(results_folder, "{}_calibration_{}.csv".format(monitor, today))
    
    # number of iterations to be used when optimising during calibration
    num_iterations = 500

    # Creates an empty ordered dictionary for the header
    calibration_header = OrderedDict()

    # flag to indicate the first file processed
    first_file = True

    for f in files:

        # check that the filepath 'exists':
        if os.path.isfile(f):

            stillbouts_ts, stillbouts_header = data_loading.load(f, source_type="HDF5", hdf5_mode="r+", hdf5_group="Still")

            # These items of the header must be the same for all files
            meta_check = ("budget", "noise_cutoff_mg", "generic_num_channels")

            # These items of the header must be added together for each file
            meta_append = ("num_final_bouts", "num_final_seconds")

            # A list of all meta items of interest

            meta_all = meta_check + meta_append

            # For the first file preserve all the meta data of interest
            if first_file:

                for item in meta_all:
                    calibration_header[item] = stillbouts_header[item]

                    calibration_ts = stillbouts_ts

                # set 'first file' flag to false...
                first_file = False

            # for subsequent files check that the parameters and number of channels in the time series match, and append the number of final bouts and final seconds
            else:
                for item in meta_check:
                    try:
                        calibration_header[item] == stillbouts_header[item]
                    except:
                        raise ValueError("Different parameters have been used.")

                for item in meta_append:
                    calibration_header[item] = int(calibration_header[item]) + int(stillbouts_header[item])

                # Append stillbouts_ts from file to calibrate_ts
                calibration_ts.append(stillbouts_ts)

        # else if file is not accessible, remove file from list of files:
        else:
            files.remove(f)

    calibration_diagnostics = triaxial_calibration.calibrate_steptwo(calibration_ts, calibration_header, calibration_statistics=False, num_iterations=num_iterations)

    calibration_diagnostics["files_used"] = files
    
    cal_dict = {**calibration_diagnostics}
    
    cal_dict["num_files_used"] = len(files)
    cal_dict["num_iterations"] = num_iterations
    
    for key in ["start_error", "end_error"]:
        cal_dict[("mf_%s" % key)] = cal_dict[key]
        cal_dict.pop(key)

    pampro_utilities.dict_write(calibration_output, monitor, cal_dict, other_index="monitor")

    # change group and permissions of calibration output file
    os.system("chgrp {} {} & chmod 770 {}".format(group, calibration_output, calibration_output))
    
    calibration_diagnostics["monitor"] = str(monitor)

    return calibration_diagnostics

#######################################################################################################################


# parse config file
settings = pd.read_csv(settings_file, dtype=str)

# parse jobs list file
jobs_df = pd.read_csv(jobs_file, dtype=str)

batch_processing.batch_process_wrapper(calibratemonitor, jobs_df, settings, job_num, num_jobs)
