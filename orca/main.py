import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cord import *

model = Model('cord/data/cord-data.csv', sd='10-01-1999')
model.create_flow_cdf()
model.find_running_WYI()
results = model.simulate() # takes a while... save results
results.to_csv('cord/data/results.csv')
results = pd.read_csv('cord/data/results.csv', index_col=0, parse_dates=True)

# calibration points (lists of pandas series)
sim = [results['DEL_HRO_pump'] / cfs_tafd,
       results['DEL_TRP_pump'] / cfs_tafd,
       results['SHA_storage'] * 1000.0,
       results['SHA_out'] / cfs_tafd,
       results['FOL_storage'] * 1000.0,
       #results['FOL_out'] / cfs_tafd,
       results['ORO_storage'] * 1000.0,
       results['ORO_out'] / cfs_tafd,
       results['YRS_storage'] * 1000.0,
       results['NML_storage'] * 1000.0,
       results['DNP_storage'] * 1000.0,
       results['EXC_storage'] * 1000.0]

obs = [model.df['HRO_pump'],
       model.df['TRP_pump'],
       model.df['SHA_storage'],
       model.df['SHA_otf'],
       model.df['FOL_storage'],
       #model.df['FOL_otf'],
       model.df['ORO_storage'],
       model.df['ORO_otf'],
       model.df['YRS_storage'],
       model.df['NML_storage'],
       model.df['DNP_storage'],
       model.df['EXC_storage']]

for f in ['W','AS-OCT']:
  i = 0
  for s,o in zip(sim,obs):
    plotter.compare(s, o, freq=f)
    plt.savefig('cord/figs/%s%d.png' % (f,i), dpi=150)
    i += 1
