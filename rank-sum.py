import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.stats import norm,skew
import seaborn as sns
import inspyred.ec.analysis as pyred
import pickle
import math
from scipy import stats


baseline = pd.read_csv('baseline-objectives.csv',index_col = 0)
with open('scenario_names_all.txt') as f:
  scenarios = f.read().splitlines()
with open('data/scenario-groupings/high_potential.txt') as f:
    scenarios_HP = f.read().splitlines()
with open('data/scenario-groupings/upper.txt') as f:
    scenarios_upper = f.read().splitlines()
with open('data/scenario-groupings/lower.txt') as f:
    scenarios_lower = f.read().splitlines()


dfHPset = pd.read_csv('hypervolume-percents/dfHPrandomset-leaveout.csv', index_col = 0)
dfUset = pd.read_csv('hypervolume-percents/dfUrandomset-leaveout.csv', index_col = 0)
dfLset = pd.read_csv('hypervolume-percents/dfLrandomset-leaveout.csv', index_col = 0)
highpot = dfUset['high-potential-random'].values
upper =dfHPset['upper-random'].values
lower =dfHPset['lower-random'].values
highpot_upper = dfHPset['high-potential-upper-random'].values
highpot_lower = dfHPset['high-potential-lower-random'].values
upper_lower = dfHPset['upper-lower-random'].values
allvar = dfHPset['all-random'].values

U,p = stats.mannwhitneyu(highpot,lower,alternative='greater')
if p < 0.05:
    print('Reject the null hypothesis')# 
else:
    print('Fail to reject')
print(p)
print(U)
test = allvar
print(np.mean(test))
print(np.median(test))
print(np.std(test))
