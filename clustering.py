import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.stats import norm,skew
import seaborn as sns
from sklearn.cluster import MiniBatchKMeans, KMeans
from sklearn.metrics.pairwise import pairwise_distances_argmin
from sklearn.datasets import make_blobs
from scipy.stats import pearsonr

with open('scenario_names_all.txt') as f:
  scenarios = f.read().splitlines()

dfgroup = pd.read_csv('data/scenario-groupings.csv',index_col = 0)
dfhv = pd.read_csv('data/adaptation-potential.csv',index_col = 0)
dfgroup['HV'] = dfhv['hypervolume']
dfgroup=(dfgroup-dfgroup.min())/(dfgroup.max()-dfgroup.min())




X = dfgroup[['HV','FNF','SWE']].values

k_means  = KMeans(n_clusters=3,random_state=0)
y_pred = k_means.fit_predict(X)

high_pot = np.where(y_pred == 0)
upper = np.where(y_pred == 1)
lower = np.where(y_pred == 2)
sc_highpot = []
sc_upper = []
sc_lower = []

for i in high_pot[0]:
	sc_highpot.append(scenarios[i])

for i in upper[0]:
	sc_upper.append(scenarios[i])

for i in lower[0]:
	sc_lower.append(scenarios[i])


with open('data/scenario-clusters/high_potential.txt','w') as f:
  f.write('\n'.join(sc_highpot))

with open('data/scenario-clusters/upper.txt','w') as f:
  f.write('\n'.join(sc_upper))

with open('data/scenario-clusters/lower.txt','w') as f:
  f.write('\n'.join(sc_lower))
