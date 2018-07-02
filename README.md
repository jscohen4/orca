## Operation of Reservoirs in California (ORCA)

A simulation mode incorporating CVP and SWP Delta exports, reservoir releases, and snow-pack to streamflow forecasts. Simulates both historical and projected scenarios.

**Requirements:** [NumPy](http://www.numpy.org/), [pandas](http://pandas.pydata.org/), [Matplotlib](http://matplotlib.org/), [and scikit-learn](http://scikit-learn.org/).

Instructions to run model in historial mode:
In main.py, set projection to False. Choose whether to calculate R2s and plot resuts. If running any data processing scripts, set process_hist_data to True. If downloading up-to-date historical data, set cdec to True. Set hist_indices to True if re-processing new data. Set hist_forcast to True if re-running forecasts. Run main.py and results will be in historical_runs_data folder.

Instructions to run in projection mode for individual cmip5 scenario:
Choose cmip5 scenario to run. Its input data file can be found in the git repository orca_cmip5_inputs. Copy desired input file from this repository to orca's input_climate_files folder. In main.py, set variable sc to desired scenario to read. Set options for plot, process_climate_data, climate_indices, and climate_forecasts. Historical  data can also be processed in same execution of script if updating forecast and gains regression coefficients is desired. After executing main.py for the projection will be in the individual_projection_runs folder.

Instructions to run in projection mode for multiple cmip5 scenarios:
Copy input files for all desired cmip5 scenarios from the orca_cmip5_inputs repository to orcas scenario_runs folder. Write names of all desired cmip5 scenarios to scenario_names.txt file in data folder.  Set climate_indices , climate_forecasts, and run_projection options in run_all_climate_projections.py. Set consolidate_outputs to put results for each scenario in the same csv files. After executing run_all_climate_projections.py, results will be in the scenario_runs and climate_results folders. 
### License
Copyright (C) 2018 CALFEWS Team. Released under the [MIT license](LICENSE.md).
