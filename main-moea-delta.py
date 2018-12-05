import numpy as np
np.warnings.filterwarnings('ignore') #to not display numpy warnings... be careful
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from orca import *
from orca.data import *
from subprocess import call
from datetime import datetime
now = datetime.now().strftime('Last modified %Y-%m-%d %H:%M:%S')
from platypus import NSGAII, Problem, Real, Integer
import random
import pickle 
from platypus.mpipool import MPIPool
from mpi4py import MPI

#Each of these booleans determines the actions that will be run by the model 

projection = False #True if running a single climate projection
calc_R2s = True #True if calculating R2s (only relevant for historical scenario)
plot = False #True if plotting outputs, need calc_R2s to also be true if plotting historical results!!!!
change_inflow_exeedance = False

#######Define a few parameters
#(in zscore bounds)
SHA_exceedance = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
ORO_exceedance = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
FOL_exceedance = {"W": 10, "AN": 10, "BN": 5, "D": 2, "C": 1}

#(btwn 0 and 1)
SHA_carryover = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
ORO_carryover = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
FOL_carryover = {"W": 10, "AN": 10, "BN": 5, "D": 2, "C": 1}

SHA_shift = 0
ORO_shift = 0
FOL_shift = 0


sc = 'access1-0_rcp45_r1i1p1' #cmip5 climate scenario to use, if projection = True
#Nothing below here should be changed!
###############################################
###############################################
#############
df = pd.read_csv('orca/data/individual_projection_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc),parse_dates = True, index_col = 0)
df = df[df.index >= '2050-10-01'] # def wrapper(SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,SHA_shift,ORO_shift,FOL_shift):
dfhist = pd.read_csv('orca/data/historical_runs_data/results.csv',parse_dates = True, index_col = 0)


# def wrapper(SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,SHA_shift,ORO_shift,FOL_shift):
def wrapper(vars):
  # SHA_exceedance = {"W": SHA_exceedance[0], "AN": SHA_exceedance[1], "BN": SHA_exceedance[2], "D": SHA_exceedance[3], "C": SHA_exceedance[4]}
  # ORO_exceedance = {"W": ORO_exceedance[0], "AN": ORO_exceedance[1], "BN": ORO_exceedance[2], "D": ORO_exceedance[3], "C": ORO_exceedance[4]}
  # FOL_exceedance = {"W": FOL_exceedance[0], "AN": FOL_exceedance[1], "BN": FOL_exceedance[2], "D": FOL_exceedance[3], "C": FOL_exceedance[4]}

  # #(btwn 0 and 1)
  # SHA_carryover = {"W": SHA_carryover[0], "AN": SHA_carryover[1], "BN": SHA_carryover[2], "D": SHA_carryover[3], "C": SHA_carryover[4]}
  # ORO_carryover = {"W": ORO_carryover[0], "AN": ORO_carryover[1], "BN": ORO_carryover[2], "D": ORO_carryover[3], "C": ORO_carryover[4]}
  # FOL_carryover = {"W": FOL_carryover[0], "AN": FOL_carryover[1], "BN": FOL_carryover[2], "D": FOL_carryover[3], "C": FOL_carryover[4]}

  # SHA_shift = SHA_shift
  # ORO_shift = ORO_shift
  # FOL_shift = FOL_shift

  SHA_exceedance = {"W": vars[0]/16.6666666667, "AN": vars[1]/16.6666666667, "BN": vars[2]/16.6666666667, "D": vars[3]/16.6666666667, "C": vars[4]/16.6666666667}
  ORO_exceedance = {"W": vars[5]/16.6666666667, "AN": vars[6]/16.6666666667, "BN": vars[7]/16.6666666667, "D": vars[8]/16.6666666667, "C": vars[9]/16.6666666667}
  FOL_exceedance = {"W": vars[10]/16.6666666667, "AN": vars[11]/16.6666666667, "BN": vars[12]/16.6666666667, "D": vars[13]/16.6666666667, "C": vars[14]/16.6666666667}

  #(btwn 0 and 1)
  SHA_carryover = {"W": max(vars[15],-vars[15])/60, "AN": max(vars[16],-vars[16])/60, "BN": max(vars[17],-vars[17])/60, "D": max(vars[18],-vars[18])/60, "C": max(vars[19],-vars[19])/60}
  ORO_carryover = {"W": max(vars[20],-vars[20])/60, "AN": max(vars[21],-vars[21])/60, "BN": max(vars[22],-vars[22])/60, "D": max(vars[23],-vars[23])/60, "C": max(vars[24],-vars[24])/60}
  FOL_carryover = {"W": max(vars[25],-vars[25])/60, "AN": max(vars[26],-vars[26])/60, "BN": max(vars[27],-vars[27])/60, "D": max(vars[28],-vars[28])/60, "C": max(vars[29],-vars[29])/60}

  SHA_shift = max(vars[30],-vars[30])
  ORO_shift = max(vars[31],-vars[31])
  FOL_shift = max(vars[32],-vars[32])

  model = Model(df,dfhist,SHA_shift, ORO_shift, FOL_shift,SHA_exceedance,ORO_exceedance,FOL_exceedance,SHA_carryover,ORO_carryover,FOL_carryover,sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
  results = model.simulate() # takes a while... save results
  # print(results.DEL_HRO_pump.sum() + results.DEL_TRP_pump.sum())
  print(results.DEL_out.sum())
  return([-results.DEL_out.sum()])
  # print(results.SHA_storage.resample('AS-OCT').last().sum()  + results.ORO_storage.resample('AS-OCT').last().sum() + results.FOL_storage.resample('AS-OCT').last().sum())
  # print('')
  # return([-results.DEL_HRO_pump.sum() - results.DEL_TRP_pump.sum(),-results.DEL_out.sum(),-(results.SHA_storage.resample('AS-OCT').last().sum()  + results.ORO_storage.resample('AS-OCT').last().sum() + results.FOL_storage.resample('AS-OCT').last().sum())])

  # results.to_csv('orca/data/individual_projection_runs/%s/%s-results.csv'%(sc,sc))
# print(wrapper([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]))
comm = MPI.COMM_WORLD # communication object
rank = comm.rank # what number processor am I?
seed = comm.rank
print(seed)
random.seed(seed)
problem = Problem(33,1)
# problem.types[:] = [Real(-3.6, 3.6), Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),
#                     Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),
#                     Real(-3.6, 3.6)Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),Real(-3.6, 3.6),
#                     Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),
#                     Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),
#                     Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),Real(0.0,1.0),
#                     Real(0.0,60.0),Real(0.0,60.0),Real(0.0,60.0),Real(0.0,60.0),Real(0.0,60.0)]
problem.types[:] = Integer(-60,60)
problem.function = wrapper

algorithm = NSGAII(problem)
algorithm.run(1000)
pickle.dump(algorithm, open("pickle_files/delta_results%s.pkl"%seed, "wb" ))

print(algorithm.result)


