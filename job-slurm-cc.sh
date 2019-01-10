#!/bin/bash -l
#SBATCH -n 97            # Total number of processors to request
#SBATCH -p high           # Queue name hi/med/lo
#SBATCH -t 70:00:00        # Run time (hh:mm:ss) - 24 hours
#SBATCH --mail-user= <your email>            # address for email notification
#SBATCH --mail-type=ALL                  # email at Begin and End of job

export PATH=<your path>:$PATH
mpirun -n 97 python main-parallel-cc.py
