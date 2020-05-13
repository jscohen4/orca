# Operation of Reservoirs in California (ORCA)

A simulation mode incorporating CVP and SWP Delta exports, reservoir releases, and snow-pack to streamflow forecasts. This branch simulates both historical and projected cmip5 scenarios.

## Requirements:
[NumPy](http://www.numpy.org/), [Pandas](http://pandas.pydata.org/), [Matplotlib](http://matplotlib.org/), [Scipy](http://www.scipy.org/), [Scikit-learn](http://scikit-learn.org/), [Seaborn](https://seaborn.pydata.org/), [MPI for Python](https://mpi4py.readthedocs.io/en/stable/), and [Platypus](https://platypus.readthedocs.io/en/latest/).

## Running ORCA

### Running ORCA in historial mode:
1. In main.py, set projection to False.
2. Choose whether to calculate R2s and plot resuts (set variables to True or False. If running any data processing scripts, set process_hist_data to True. If downloading up-to-date historical data, set cdec to True. Set hist_indices to True if re-processing new data. Set hist_forcast to True if re-running forecasts.
3. Run main.py and results will be in historical_runs_data folder.

### Running ORCA in projection mode for individual cmip5 scenario:
1. Choose cmip5 scenario to run. Its input data file can be found in the git repository [orca_cmip5_inputs](https://github.com/jscohen4/orca_cmip5_inputs). Copy desired input file from this repository to orca's input_climate_files folder.
2. In main.py, set variable sc to desired scenario to read. Set options for plot, process_climate_data, climate_indices, and climate_forecasts. Historical  data can also be processed in same execution of script if updating forecast and gains regression coefficients is desired.
3. After executing main.py for the projection will be in the individual_projection_runs folder.

### Running ORCA in projection mode for multiple cmip5 scenarios:
1. Copy input files for all desired cmip5 scenarios from the [orca_cmip5_inputs](https://github.com/jscohen4/orca_cmip5_inputs) repository to ORCA's scenario_runs folder.
2. Write names of all desired cmip5 scenarios to scenario_names.txt file in data folder.
3. Set climate_indices , climate_forecasts, and run_projection options in run_all_climate_projections.py. Set consolidate_outputs to format results for each scenario in the same csv files.
4. After executing run_all_climate_projections.py, results will be in the scenario_runs and climate_results folders.

## License
Copyright (C) 2018 CALFEWS Team. Released under the [MIT license](LICENSE.md).
