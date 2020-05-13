import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.stats import norm,skew
import seaborn as sns
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
dfallset = dfHPset + dfUset + dfLset

testset = dfHPset
highpot = testset['high-potential-random'].values
upper =testset['upper-random'].values
lower =testset['lower-random'].values
highpot_upper = testset['high-potential-upper-random'].values
highpot_lower = testset['high-potential-lower-random'].values
upper_lower = testset['upper-lower-random'].values
allvar = testset['all-random'].values

U,p = stats.mannwhitneyu(trainset,testset,alternative='greater')
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
