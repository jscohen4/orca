# Operation of Reservoirs in California (ORCA)
This repository contains all code corresponding to methods used to generate data for figures in the following paper:

Cohen, J.S., Zeff, H.B., Herman, J.D.,  (2020). How do the hydrologic properties of training scenarios influence the robustness of reservoir policy adaptations to climate uncertainty? (submitted manuscript).

## Requirements:
[NumPy](http://www.numpy.org/), [Pandas](http://pandas.pydata.org/), [Matplotlib](http://matplotlib.org/), [Scipy](http://www.scipy.org/), [Scikit-learn](http://scikit-learn.org/), [Seaborn](https://seaborn.pydata.org/), [MPI for Python](https://mpi4py.readthedocs.io/en/stable/), and [Platypus](https://platypus.readthedocs.io/en/latest/).

## Methods:
The following instrutctions correspond to the methods section of the paper. 

### Section 3.2.1 part 1: Baseline Policy Performance:
  1. Copy input files for all cmip5 scenarios from the [orca_cmip5_inputs](https://github.com/jscohen4/orca_cmip5_inputs) repository to ORCA's scenario_runs folder.
  
  2. Run `run_all_climate_projections.py` to simulate baseline performance in series. Alternatively, simulations can be run in parallel, using 97 processors (1 for each scenario), with the script `main-      parallel-base-cc.py`. Use the following command in a shell script: `mpirun -n 97 python main-parallel-base-cc.py`. In either case, baseline results will appear in the `orca/data/scenario_runs` and `orca/data/climate_results` folders.

### Section 3.2.1 part 2: Perfect Foresight Optimization:
`perfect-foresight.py` will run the perfect foresight optimization for each scenario using NSGAIII. The script should be run in parallel with a desired number of processors (num_processors) specified in a shell script with the command `mpirun -n num_processors python perfect-foresight.py`. Pickle files contatininf decision variables and objectives output by each optimization will appear in the folders `data/perfect-foresight\objectives` and `data/perfect-foresight\variables`.

### Section 3.2.2 Hypervolume Metric:
Run `adaptation-potential.py` to calculate the adaptation potential for each scenario. Results will appear in the file `data/adaptation-potential.csv`.

### Section 3.3 Clustering:
`clustering.py` generates the scenario groupings for the high-potential, low-potential wet, and low-potential dry clusters. The file outputs appear in folder `data/scenario-clusters`

### Section 3.4 part 1: Training:

## License
Copyright (C) 2020 CALFEWS Team. Released under the [MIT license](LICENSE.md).
