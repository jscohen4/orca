import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cord import *

model = Model('cord/data/cord-data.csv', sd='10-01-1999')
results = model.simulate() # takes a while... save results
results.to_csv('cord/data/results.csv')
# # results = pd.read_csv('cord/data/results.csv', index_col=0, parse_dates=True)

# calibration points (lists of pandas series)
sim = [results['DEL_HRO_pump'] / cfs_tafd,
       results['DEL_TRP_pump'] / cfs_tafd,
       (results['DEL_HRO_pump'] + results['DEL_TRP_pump']) / cfs_tafd,
       results['SHA_storage'], 
       results['SHA_out'] / cfs_tafd,
       results['FOL_storage'],
       results['FOL_out'] / cfs_tafd,
       results['ORO_storage'],
       results['ORO_out'] / cfs_tafd,
       results['DEL_in'] / cfs_tafd,
       results['DEL_out'] / cfs_tafd]

obs = [model.df['HRO_pump'],
       model.df['TRP_pump'],
       (model.df['HRO_pump'] + model.df['TRP_pump']),
       model.df['SHA_storage'],
       model.df['SHA_out'],
       model.df['FOL_storage'],
       model.df['FOL_out'],
       model.df['ORO_storage'],
       model.df['ORO_out'],
       model.df['DeltaIn'],
       model.df['DeltaOut']]
plotter.Rsquares(sim,obs)

for f in ['D','W','M','AS-OCT']:
  i = 0
  for s,o in zip(sim,obs):
    plotter.compare(s, o, freq=f)
    plt.savefig('cord/figs/%s%d.png' % (f,i), dpi=150)
    plt.close()
    i += 1
