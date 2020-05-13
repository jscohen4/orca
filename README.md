# Operation of Reservoirs in California (ORCA)
This repository contains all code corresponding to methods used to generate data for figures in the following paper:

Cohen, J.S., Zeff, H.B., Herman, J.D.,  (2020). How do the hydrologic properties of training scenarios influence the robustness of reservoir policy adaptations to climate uncertainty? (submitted manuscript).

## Requirements:
[NumPy](http://www.numpy.org/), [Pandas](http://pandas.pydata.org/), [Matplotlib](http://matplotlib.org/), [Scipy](http://www.scipy.org/), [Scikit-learn](http://scikit-learn.org/), [Seaborn](https://seaborn.pydata.org/), [MPI for Python](https://mpi4py.readthedocs.io/en/stable/), and [Platypus](https://platypus.readthedocs.io/en/latest/).

## Running ORCA
ORCA (Operations of Reservoirs in California): A simulation mode incorporating Central Valley Project and State Water Project Delta exports, reservoir releases, and snowpack to streamflow forecasts.


## Methods:
The following instrutctions correspond to the methods section of the paper. 

### Section 3.2.1 part 1: Baseline Policy Performance:
1. Copy input files for all desired cmip5 scenarios from the [orca_cmip5_inputs](https://github.com/jscohen4/orca_cmip5_inputs) repository to ORCA's scenario_runs folder.
2. Write names of all desired cmip5 scenarios to scenario_names.txt file in data folder.
3. Set climate_indices , climate_forecasts, and run_projection options in run_all_climate_projections.py. Set consolidate_outputs to format results for each scenario in the same csv files.
4. After executing run_all_climate_projections.py, results will be in the scenario_runs and climate_results folders.

## License
Copyright (C) 2020 CALFEWS Team. Released under the [MIT license](LICENSE.md).
