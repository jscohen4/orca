# Operation of Reservoirs in California (ORCA)
This repository contains all code corresponding to methods used to generate data for figures in the following paper:

Cohen, J.S., Zeff, H.B., Herman, J.D.,  (2020). How do the hydrologic properties of training scenarios influence the robustness of reservoir policy adaptations to climate uncertainty? (submitted manuscript).

## Requirements:
[NumPy](http://www.numpy.org/), [Pandas](http://pandas.pydata.org/), [Matplotlib](http://matplotlib.org/), [Scipy](http://www.scipy.org/), [Scikit-learn](http://scikit-learn.org/), [Seaborn](https://seaborn.pydata.org/), [MPI for Python](https://mpi4py.readthedocs.io/en/stable/), and [Platypus](https://platypus.readthedocs.io/en/latest/).

## Methods:
The following instrutctions correspond to the methods section of the paper. 

### Section 3.2.1 part 1: Baseline Policy Performance:
  1. Copy input files for all cmip5 scenarios from the [orca_cmip5_inputs](https://github.com/jscohen4/orca_cmip5_inputs) repository to ORCA's scenario_runs folder.
  
  2. Run `run_all_climate_projections.py` to simulate baseline performance in series. Alternatively, simulations can be run in parallel, using 97 processors (1 for each scenario), with the script `main-parallel-base-cc.py`. Use the following command in a shell script: `mpirun -n 97 python main-parallel-base-cc.py`. In either case, baseline results will appear in the `orca/data/scenario_runs` and `orca/data/climate_results` folders.

### Section 3.2.1 part 2: Perfect Foresight Optimization:
`perfect-foresight.py` will run the perfect foresight optimization for each scenario using NSGAIII. The script should be run in parallel with 128 processors, although this number can be adjusted. Specity in a shell script with the command `mpirun -n num_processors python perfect-foresight.py`. Pickle files contatining decision variable and objective output from each optimization will appear in the folders `data/perfect-foresight\objectives` and `data/perfect-foresight\variables`.

### Section 3.2.2 Hypervolume Metric:
Run `adaptation-potential.py` to calculate the adaptation potential for each scenario. Results will appear in the file `data/adaptation-potential.csv`.

### Section 3.3 Clustering:
`clustering.py` generates the scenario groupings for the high-potential, low-potential wet, and low-potential dry clusters. The file outputs appear in folder `data/scenario-clusters`

### Section 3.4 part 1: Training:
`orca-MOEA-training.py` will run the optimization procedure using one individual training set. The training set can be specified in line 40. Output file names in lines 208 and 209 should be adjusted accordingly. The script is run using 128 processors, although this number can be adjusted. Specify in a shell script with the command `mpirun -n num_processors python orca-MOEA-training.py`.

### Section 3.4 part 2: Testing:
In `orca-testing.py`, specify an individual testing set in lines 59 and 139. Change output file names in line 135 accordingly. The script must be run in parallel with the number of processors equal to the number of scenarios in the test set. These numbers are displayed in Table 1 in the manuscript. Specify training set policy file in line 134. Run with a shell script containing the command `mpirun -n num_processors python orca-testing.py`, with the number of processors specified. 

### Section 3.5.1 Hypervolume Robustness Metric:
 `hypervolume-robustness.py` calculates the hypervolume robustness metric for each scenario in each training-test set combination. The script is run in serial and output files will apear in the folder `data/robustness`. 
 
### Section 3.5.2 Rank-Sum Tests:
`runk-sum.py` prints results from the Mann-Whitney U test. Specify the training set in line 28. Specify the two testing sets being compared in line 37. 

## License
Copyright (C) 2020 CALFEWS Team. Released under the [MIT license](LICENSE.md).
