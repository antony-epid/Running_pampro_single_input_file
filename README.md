# Running_pampro_single_input_file

This repo gives an example of how to perform an analysis of accelerometer data (in HDF5 format) using pampro library in python3.

The pampro source code can be obtained from the repository https://github.com/MRC-Epid/pampro

The source code contains the modules needed to process raw accelerometer files, convert them into HDF5 format, and carry out the computation for the analyses.

A script [standardanalysis_file.py](standardanalysis_file.py) is provided in this repository and used as an example of how to perform an analysis of an HDF5 accelerometer data.

The following are the most simple steps to carry out such an analysis:  
* Download pampro source code and install the package following the instruction in the repository above. Make sure your python libraries have the dependencies such as pandas, numpy, etc.
* In your chosen work directory, create three subdirectories with the following names : _plots, _results, _logs
* Copy the script [standardanalysis_file.py](standardanalysis_file.py) from this repository to the work directory. 
* Copy the setting file [standardanalysis_settings.csv](standardanalysis_settings.csv) from this repository to the work directory. 
* Copy the input file (HDF5 file) from the google bucket [g-watch-data-in-dev/test-upload/PSA-0001-9_Axivity_52029_1643871614_20220322-160559_e89b4390-a9b6-11ec-9901-b38bc08f2f77_100Hz.hdf5](https://console.cloud.google.com/storage/browser/_details/g-watch-data-in-dev/test-upload/PSA-0001-9_Axivity_52029_1643871614_20220322-160559_e89b4390-a9b6-11ec-9901-b38bc08f2f77_100Hz.hdf5?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&authuser=1&cloudshell=false&project=dark-gateway-416214) to the work directory.
* In your work directory, run the script using the setting file and the input HDF5 file with the following command line  
`python standardanalysis_file.py --settings_file standardanalysis_settings.csv --filenames PSA-0001-9_Axivity_52029_1643871614_20220322-160559_e89b4390-a9b6-11ec-9901-b38bc08f2f77_100Hz.hdf5`

If the run is successful, it will produce outputs all the subdirectories (_plots, _results, _logs).


# Multi-cores run
We can test the script for a simple multi-process run using the following steps :
* Copy the following input files (HDF5 files) to the work directory from the google bucket: [g-watch-data-in-dev/test-upload/PSA-0002-7_Axivity_53826_1643869024_20220322-161106_9f7a4700-a9b7-11ec-bc9d-4de31ddeb071_100Hz.hdf5](https://console.cloud.google.com/storage/browser/_details/g-watch-data-in-dev/test-upload/PSA-0002-7_Axivity_53826_1643869024_20220322-161106_9f7a4700-a9b7-11ec-bc9d-4de31ddeb071_100Hz.hdf5?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&authuser=1&cloudshell=false&project=dark-gateway-416214) 
and [g-watch-data-in-dev/test-upload/PSA-0003-5_Axivity_52549_1643869085_20220322-161054_9885bca0-a9b7-11ec-a15e-739a873c17f5_100Hz.hdf5](https://console.cloud.google.com/storage/browser/_details/g-watch-data-in-dev/test-upload/PSA-0003-5_Axivity_52549_1643869085_20220322-161054_9885bca0-a9b7-11ec-a15e-739a873c17f5_100Hz.hdf5?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&authuser=1&cloudshell=false&project=dark-gateway-416214)

* Run it with the following command line:
`python standardanalysis_file.py --settings_file standardanalysis_settings.csv --filenames PSA-0002-7_Axivity_53826_1643869024_20220322-161106_9f7a4700-a9b7-11ec-bc9d-4de31ddeb071_100Hz.hdf5 PSA-0003-5_Axivity_52549_1643869085_20220322-161054_9885bca0-a9b7-11ec-a15e-739a873c17f5_100Hz.hdf5`

We may add more input files to the command line argument as long as the total does not exceed the number of available cores in your system.  
