import numpy as np
np.warnings.filterwarnings('ignore') #to not display numpy warnings... be careful
import pandas as pd
import pickle
from orca import *
from orca.data import *
from subprocess import call
from datetime import datetime
now = datetime.now().strftime('Last modified %Y-%m-%d %H:%M:%S')
from platypus import NSGAIII, DTLZ2, EpsMOEA, PoolEvaluator, Problem, Real, Integer
from platypus.mpipool import MPIPool
import random
import pickle 
from mpi4py import MPI
import sys
import logging
logging.basicConfig(level=logging.INFO)

with open('orca/data/scenario_names_all.txt') as f:
  scenarios_all = f.read().splitlines()
with open('data/scenario-groupings/high-potential-random.txt') as f:
    scenarios_HP_train = f.read().splitlines()
with open('data/scenario-groupings/upper-random.txt') as f:
    scenarios_upper_train = f.read().splitlines()
with open('data/scenario-groupings/lower-random.txt') as f:
    scenarios_lower_train = f.read().splitlines()
with open('data/scenario-groupings/high_potential.txt') as f:
    scenarios_HP_all = f.read().splitlines()
with open('data/scenario-groupings/upper.txt') as f:
    scenarios_upper_all = f.read().splitlines()
with open('data/scenario-groupings/lower.txt') as f:
    scenarios_lower_all = f.read().splitlines()


scenarios_HP_test = []
scenarios_upper_test = []
scenarios_lower_test = []


for sc in scenarios_HP_all:
    if sc not in scenarios_HP_train:
        scenarios_HP_test.append(sc)

for sc in scenarios_upper_all:
    if sc not in scenarios_upper_train:
        scenarios_upper_test.append(sc)

for sc in scenarios_lower_all:
    if sc not in scenarios_lower_train:
        scenarios_lower_test.append(sc)

scenarios_HP_test_num = []
for i,sc in enumerate(scenarios_all):
  if sc in scenarios_HP_test:
    scenarios_HP_test_num.append(i)

scenarios_upper_test_num = []
for i,sc in enumerate(scenarios_all):
  if sc in scenarios_upper_test:
    scenarios_upper_test_num.append(i)

scenarios_lower_test_num = []
for i,sc in enumerate(scenarios_all):
  if sc in scenarios_lower_test:
    scenarios_lower_test_num.append(i)

with open('orca/data/scenario_names_all.txt') as f:
  scenarios_all = f.read().splitlines()


scenarios_tested = scenarios_HP_train

for sc in scenarios_all: 
  if sc not in scenarios_tested:
    scenarios_cross.append(sc)
for sc in scenarios_cross:
  call(['mkdir','cross_validation_time_series-curtail/%s'%sc])
  
SHA_baseline = pd.read_csv('orca/data/climate_results_baseline/SHA_storage.csv',parse_dates = True, index_col = 0)
SHA_baseline = SHA_baseline[(SHA_baseline.index >= '2069-10-01') & (SHA_baseline.index <= '2099-10-01')]

ORO_baseline = pd.read_csv('orca/data/climate_results_baseline/ORO_storage.csv',parse_dates = True, index_col = 0)
ORO_baseline = ORO_baseline[(ORO_baseline.index >= '2069-10-01') & (ORO_baseline.index <= '2099-10-01')]

FOL_baseline = pd.read_csv('orca/data/climate_results_baseline/FOL_storage.csv',parse_dates = True, index_col = 0)
FOL_baseline = FOL_baseline[(FOL_baseline.index >= '2069-10-01') & (FOL_baseline.index <= '2099-10-01')]

def evaluate(var,scenarios_cross,scenario_num,rank):  
  flood_risk = 0
  carryover_storage = 0
  pumping = 0
  delta_outflow = 0
  for sc,sc_num in zip(scenarios_cross,scenario_num):
    df_forc = pd.read_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc),parse_dates = True, index_col = 0)
    df_forc = df_forc[(df_forc.index >= '2069-10-01') & (df_forc.index <= '2099-10-01')]
    dfhist = pd.read_csv('orca/data/historical_runs_data/results.csv',parse_dates = True, index_col = 0)

    SHA_exceedance = {"W": var[0]/16.6666666667, "AN": var[1]/16.6666666667, "BN": var[2]/16.6666666667, "D": var[3]/16.6666666667, "C": var[4]/16.6666666667}
    ORO_exceedance = {"W": var[5]/16.6666666667, "AN": var[6]/16.6666666667, "BN": var[7]/16.6666666667, "D": var[8]/16.6666666667, "C": var[9]/16.6666666667}
    FOL_exceedance = {"W": var[10]/16.6666666667, "AN": var[11]/16.6666666667, "BN": var[12]/16.6666666667, "D": var[13]/16.6666666667, "C": var[14]/16.6666666667}

    #(btwn 0 and 1)
    SHA_carryover = {"W": max(var[15],-var[15])/60, "AN": max(var[16],-var[16])/60, "BN": max(var[17],-var[17])/60, "D": max(var[18],-var[18])/60, "C": max(var[19],-var[19])/60}
    ORO_carryover = {"W": max(var[20],-var[20])/60, "AN": max(var[21],-var[21])/60, "BN": max(var[22],-var[22])/60, "D": max(var[23],-var[23])/60, "C": max(var[24],-var[24])/60}
    FOL_carryover = {"W": max(var[25],-var[25])/60, "AN": max(var[26],-var[26])/60, "BN": max(var[27],-var[27])/60, "D": max(var[28],-var[28])/60, "C": max(var[29],-var[29])/60}

    SHA_shift = max(var[30],-var[30])
    ORO_shift = max(var[31],-var[31])
    FOL_shift = max(var[32],-var[32])
    wyt_exc = int(min(abs(var[33]),59)/5)
    model = Model(df_forc,dfhist,SHA_baseline[sc], ORO_baseline[sc], FOL_baseline[sc],SHA_shift, ORO_shift, FOL_shift,SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,wyt_exc,sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
    results = model.simulate() # takes a while... save results
    dfresults = pd.DataFrame(index = results.index)
    dfresults['SHA_out'] = results['SHA_out']
    dfresults['SHA_spill'] = results['SHA_spill']
    dfresults['SHA_storage'] = results['SHA_storage']
    dfresults['ORO_out'] = results['ORO_out']
    dfresults['ORO_spill'] = results['ORO_spill']
    dfresults['ORO_storage'] = results['ORO_storage']
    dfresults['DEL_HRO_pump'] = results['DEL_HRO_pump']
    dfresults['DEL_TRP_pump'] = results['DEL_TRP_pump']
    dfresults['DEL_out'] = results['DEL_out']
    dfresults.to_csv('data/tested-policy-timeseries/high_potential/scenario%s/policy%s.csv'%(sc_num,rank))

variables = pickle.load(open("data/training-outputs/high-potential-variables.pkl", "rb"))
comm = MPI.COMM_WORLD # communication object
rank = comm.rank# what number processor am I?
var = variables[test_policies[rank]] #+84
objs = evaluate(var,scenarios_upper_test,scenarios_upper_test_num,rank)

