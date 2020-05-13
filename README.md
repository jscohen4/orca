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
  1. Copy input files for all cmip5 scenarios from the [orca_cmip5_inputs](https://github.com/jscohen4/orca_cmip5_inputs) repository to ORCA's scenario_runs folder.
  
  2a. Run `run_all_climate_projections.py` to simulate baseline performance in series. 
  
  2b: Alternatively, simulations can be run in parallel, using 97 processors (1 for each scenario), with the script `main-      parallel-base-cc.py`. Use the following command in a shell script: `mpirun -n 97 python main-parallel-base-cc.py`.
In either case, baseline results will appear in the `orca/data/scenario_runs` and `orca/data/climate_results` folders.

## License
Copyright (C) 2020 CALFEWS Team. Released under the [MIT license](LICENSE.md).
