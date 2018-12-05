#!/bin/bash -l
#SBATCH -n 10            # Total number of processors to request
#SBATCH -p high           # Queue name hi/med/lo
#SBATCH -t 70:00:00        # Run time (hh:mm:ss) - 24 hours
#SBATCH --mail-user=joncohen@ucdavis.edu              # address for email notification
#SBATCH --mail-type=ALL                  # email at Begin and End of job

export PATH=/group/hermangrp/miniconda3/bin:$PATH
mpirun -n 10 python main-moea-delta.py
