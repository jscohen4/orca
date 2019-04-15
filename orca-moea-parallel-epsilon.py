import numpy as np
np.warnings.filterwarnings('ignore') #to not display numpy warnings... be careful
import pandas as pd
import pickle
from orca import *
from orca.data import *
from subprocess import call
from datetime import datetime
now = datetime.now().strftime('Last modified %Y-%m-%d %H:%M:%S')
from platypus import NSGAII, DTLZ2, EpsMOEA, PoolEvaluator, Problem, Real, Integer
from platypus.mpipool import MPIPool
import random
import pickle 
from mpi4py import MPI
import sys
import logging
logging.basicConfig(level=logging.INFO)

seed = 0
random.seed(seed)

with open('orca/data/scenario_names_all.txt') as f:
  scenarios = f.read().splitlines()

num_selections = 70
scenarios = random.sample(scenarios, num_selections)
# pickle.dump(scenarios, open("pickle_files/scenarios_run1.pkl", "wb" ))
with open('pickle_files/scenario_list-eps-100-30.txt', 'w') as scenario_list:  
    scenario_list.writelines("%s\n" % scenario for scenario in scenarios)

df_dict = {}
for sc in scenarios:
  df_dict[sc] = pd.read_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc),parse_dates = True, index_col = 0)
  df_dict[sc] = df_dict[sc][(df_dict[sc].index >= '2069-10-01') & (df_dict[sc].index <= '2099-10-01')] # def wrapper(SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,SHA_shift,ORO_shift,FOL_shift):
dfhist = pd.read_csv('orca/data/historical_runs_data/results.csv',parse_dates = True, index_col = 0)

class ORCA_wrapper(Problem):
    
    def __init__(self, ndec, nobj, lb, ub):
      super(ORCA_wrapper, self).__init__(ndec, nobj)
      self.types[:] = Integer(lb,ub)

    def evaluate(self, solution):
      var = solution.variables
      var = tuple(var)
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
      flood_risk = 0
      carryover_storage = 0
      pumping = 0
      delta_outflow = 0

      for key,df in df_dict.items():
        model = Model(df,dfhist,SHA_shift, ORO_shift, FOL_shift,SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
        results = model.simulate() # takes a while... save results
        
        SHA_out = (results.SHA_spill).values
        for out in SHA_out:
          if out > 0:
            flood_risk += out
        
        ORO_out = (results.ORO_spill).values
        for out in ORO_out:
          if out > 0:
            flood_risk += out

        FOL_out = (results.FOL_spill).values
        for out in FOL_out:
          if out > 0:
            flood_risk += out

        carryover_storage += -(results.SHA_storage.resample('AS-OCT').last().sum()  + results.ORO_storage.resample('AS-OCT').last().sum() + results.FOL_storage.resample('AS-OCT').last().sum())
        delta_outflow += -results.DEL_out.sum()
        pumping += -(results.DEL_HRO_pump.sum() + results.DEL_TRP_pump.sum())

      solution.objectives[0] = flood_risk
      solution.objectives[1] = carryover_storage
      solution.objectives[2] = delta_outflow
      solution.objectives[3] = pumping
      #print(solution.objectives)
# def schaffer(x):
#     return [x[0]**2, (x[0]-2)**2]

if __name__ == "__main__":
  problem = ORCA_wrapper(33,4,-60,60)

  # define the problem definition
  pool = MPIPool()

  # only run the algorithm on the master process
  if not pool.is_master():
      pool.wait()
      sys.exit(0)

  # instantiate the optimization algorithm to run in parallel
  with PoolEvaluator(pool) as evaluator:
      algorithm = EpsMOEA(problem, {"epsilons":[30000,60000,30000,30000]}, evaluator=evaluator)
      algorithm.run(100)
  int1 =  Integer(-60,60)
  results_list = []
  variable_list = []
  for results in algorithm.result:
    results_list.append(results.objectives[:])
    var_list = []
    for var in results.variables:
      var_list.append(int1.decode(var))
    variable_list.append(var_list)
  pickle.dump(results_list, open("pickle_files/objective_all-eps-100-30.pkl", "wb" ))
  pickle.dump(variable_list, open("pickle_files/vriable_all-eps-100-30.pkl", "wb" ))


  pool.close()
