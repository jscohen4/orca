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

seed = 0
random.seed(seed)

with open('orca/data/scenario_names_all.txt') as f:
  scenarios = f.read().splitlines()
df_dict = {}

dfhist = pd.read_csv('orca/data/historical_runs_data/results.csv',parse_dates = True, index_col = 0)

SHA_baseline = pd.read_csv('orca/data/climate_results_baseline/SHA_storage.csv',parse_dates = True, index_col = 0)
SHA_baseline = SHA_baseline[(SHA_baseline.index >= '2069-10-01') & (SHA_baseline.index <= '2099-10-01')]

ORO_baseline = pd.read_csv('orca/data/climate_results_baseline/ORO_storage.csv',parse_dates = True, index_col = 0)
ORO_baseline = ORO_baseline[(ORO_baseline.index >= '2069-10-01') & (ORO_baseline.index <= '2099-10-01')]

FOL_baseline = pd.read_csv('orca/data/climate_results_baseline/FOL_storage.csv',parse_dates = True, index_col = 0)
FOL_baseline = FOL_baseline[(FOL_baseline.index >= '2069-10-01') & (FOL_baseline.index <= '2099-10-01')]
i = 0
sc = scenarios[i]
dfsc = pd.read_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc),parse_dates = True, index_col = 0)
dfsc = dfsc[(dfsc.index >= '2069-10-01') & (dfsc.index <= '2099-10-01')]
scenario = sc
class ORCA_wrapper(Problem):

    def __init__(self, ndec, nobj, lb, ub):#,dfsc,scenario):
      super(ORCA_wrapper, self).__init__(ndec, nobj)
      self.types[:] = Integer(lb,ub)
      
      self.SHA_max_volume = 36000 #af/day
      self.SHA_max_capacity = 710 #MW

      self.ORO_max_volume = 33620 #af/day
      self.ORO_max_capacity = 819 #MW
      
      self.FOL_max_volume = 16000 #af/day
      self.FOL_max_capacity = 210 #MW
      self.dfsc = dfsc
      self.scenario = scenario

    def storage_head_FOL(self,S):
      return 3.636e-28 * S**5 - 1.197e-21* S**4 + 1.588e-15* S**3 - 1.108e-09* S**2 + 0.0005281 *S + 294.5

    def storage_head_ORO(self,S):
      return 4.003e-18 *S**3 - 4.245e-11 *S**2+ 0.0002172 *S + 486.1

    def storage_head_SHA(self,S):
      return 7.961e-19 * S**3 - 1.23e-11 *S**2 + 9.659e-05 * S + 807.2

    def Shasta_EF(self,Head):
      return 1.045 * ((0.83522 * (Head-630)) + 30.5324) #Gross_Head = Head - 580 (580 is height of turbines)

    def Oroville_EF(self,Head):
      T = len(Head)
      Gross_Head = Head - 300#turbine elev is 225
      Hyatt_EF = np.zeros(T)
      for t in range(0,T):
        if Gross_Head[t] < 505:
          A = 8096.1
          B = -1.28
          C = 100
        elif Gross_Head[t] >= 505 and Gross_Head[t] < 555:
          A = 2930.2226
          B = -1.13794
          C = 48
        elif Gross_Head[t] >= 555:
          A = 398.667
          B = -0.84085
          C = -15
        Hyatt_EF = 1000/((A * (Gross_Head[t] + C) ** B)) 
      return Hyatt_EF

    def Folsom_EF(self,Head):
      return  (0.92854 * (Head-175)) - 16.282 #Gross_Head = Head - 125 (125 is height of turbines) 


    def calc_energy(self,EF,Q,max_volume,max_capacity):
      T = len(Q)
      energy = np.zeros(T)
      for t in range(0,T):
        energy[t] = min(max_capacity,EF[t]*min(Q[t],max_volume)*(1/(24*1000)))
      return energy #returns MW/day

    def evaluate(self, solution):
      var = solution.variables
      var = tuple(var)
      #between 0 and 3.6
      SHA_exceedance = {"W": var[0]/16.6666666667, "AN": var[1]/16.6666666667, "BN": var[2]/16.6666666667, "D": var[3]/16.6666666667, "C": var[4]/16.6666666667}
      ORO_exceedance = {"W": var[5]/16.6666666667, "AN": var[6]/16.6666666667, "BN": var[7]/16.6666666667, "D": var[8]/16.6666666667, "C": var[9]/16.6666666667}
      FOL_exceedance = {"W": var[10]/16.6666666667, "AN": var[11]/16.6666666667, "BN": var[12]/16.6666666667, "D": var[13]/16.6666666667, "C": var[14]/16.6666666667}

      #(btwn 0 and 1)
      SHA_carryover = {"W": max(var[15],-var[15])/60, "AN": max(var[16],-var[16])/60, "BN": max(var[17],-var[17])/60, "D": max(var[18],-var[18])/60, "C": max(var[19],-var[19])/60}
      ORO_carryover = {"W": max(var[20],-var[20])/60, "AN": max(var[21],-var[21])/60, "BN": max(var[22],-var[22])/60, "D": max(var[23],-var[23])/60, "C": max(var[24],-var[24])/60}
      FOL_carryover = {"W": max(var[25],-var[25])/60, "AN": max(var[26],-var[26])/60, "BN": max(var[27],-var[27])/60, "D": max(var[28],-var[28])/60, "C": max(var[29],-var[29])/60}

      SHA_shift =max(var[30],-var[30])
      ORO_shift = max(var[31],-var[31])
      FOL_shift = max(var[32],-var[32])
      wyt_exc = int(min(abs(var[33]),59)/5)

      flood_risk = 0
      carryover_storage = 0
      pumping = 0
      delta_outflow = 0
      hydropower = 0

      model = Model(self.dfsc, dfhist, SHA_baseline[self.scenario], ORO_baseline[self.scenario], FOL_baseline[self.scenario], SHA_shift, ORO_shift, FOL_shift,SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,wyt_exc,sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
      raise Exception('model runs')
      results = model.simulate() # takes a while... save results

      results['SHA_head'] = self.storage_head_SHA(results.SHA_storage.values*1000)
      results['ORO_head'] = self.storage_head_ORO(results.ORO_storage.values*1000)
      results['FOL_head'] = self.storage_head_FOL(results.FOL_storage.values*1000)

      results['SHA_EF'] = self.Shasta_EF(results.SHA_head.values)
      results['ORO_EF'] = self.Oroville_EF(results.ORO_head.values)
      results['FOL_EF'] = self.Folsom_EF(results.FOL_head.values)

      results['SHA_energy'] = self.calc_energy(results.SHA_EF.values,results.SHA_out.values*1000,self.SHA_max_volume,self.SHA_max_capacity)*24/1000
      results['ORO_energy'] = self.calc_energy(results.ORO_EF.values,results.ORO_out.values*1000,self.ORO_max_volume,self.ORO_max_capacity)*24/1000
      results['FOL_energy'] = self.calc_energy(results.FOL_EF.values,results.FOL_out.values*1000,self.FOL_max_volume,self.FOL_max_capacity)*24/1000

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
      hydropower += -(results.SHA_energy.sum()+results.ORO_energy.sum()+results.FOL_energy.sum())
      solution.objectives[0] = flood_risk
      solution.objectives[1] = carryover_storage
      solution.objectives[2] = delta_outflow
      solution.objectives[3] = pumping
      solution.objectives[4] = hydropower

if __name__ == "__main__":
  # for i,scenario in enumerate(scenarios):
    seed = 0
    random.seed(seed)

    problem = ORCA_wrapper(34,5,-60,60)#,dfsc,sc)

    # define the problem definition
    pool = MPIPool()

    # only run the algorithm on the master process
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    # instantiate the optimization algorithm to run in parallel
    with PoolEvaluator(pool) as evaluator:
        algorithm = NSGAIII(problem, divisions_outer= 8,evaluator=evaluator)
        algorithm.run(50000)

    int1 =Integer(-60,60)
    results_list = []
    variable_list = []
    print(i)
    for results in algorithm.result:
      results_list.append(results.objectives[:])
      var_list = []
      for var in results.variables:
        var_list.append(int1.decode(var))
      variable_list.append(var_list)
    pickle.dump(results_list, open("data/perfect-foresight/objectives/PF-objectives-sc%s.pkl"%i, "wb" ))
    pickle.dump(variable_list, open("data/perfect-foresight/variables/PF-obj-variables-sc%s.pkl"%i, "wb" ))


    pool.close()


