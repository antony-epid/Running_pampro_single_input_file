# Running_pampro_single_input_file

This repo gives an example of how to perform an analysis of accelerometer data (in HDF5 format) using the pamrpo script in python3.
The pampro source code can be obtained from the repository.... It contains the libraries needed to process the files and carry out the computation for the analyses.
In this repository, a script standardanalysis_file.py is provided as an example of how to perform an analysis of accelerometer data.

There are somes steps to 
1. Download pampro source code and install the package following the instruction in the repository above. Make sure your python libraries have the dependencies such as pandas, numpy, etc.
2. In your chosen work directory, create 3 subdirectories with the following names : _plots, _results, _logs
3. Copy the script standardanalysis_file.py and the setting file standardanalysis_settings.csv from this repository to the work directory. 
4. Copy the setting file standardanalysis_settings.csv from this repository to the work directory. 
5. Copy the input file (HDF5 file) from the google bucket .... to the work directory.
6. You can run the script from the command line  

python standardanalysis_file.py  standardanalysis_settings.csv  PSA-0006-8_Axivity_25158_1643951574_20220228-163006_a1e05710-9870-11ec-9b91-295965cdd629_100Hz.hdf5
The first argument in the command line is the setting file whereas the second is the input HDF5 file.

If the run is successful, it will produce outputs all the subdirectories (_plots, _results, _logs).
