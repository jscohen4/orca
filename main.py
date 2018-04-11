import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from orca import *

scenario = True

if not scenario:
  model = Model('orca/data/orca-data-forecasted.csv', 'orca/data/orca-data-forecasted.csv',sd='10-01-1999',scenario = False, sim_gains = True) #beacuse of rolling calc in gains, we start on 10th day of
  results = model.simulate() # takes a while... save results
  results.to_csv('orca/data/results.csv')

if scenario:
  model = Model('orca/data/orca-data-climate-forecasted.csv', 'orca/data/results.csv',sd='10-01-1999',scenario = True, sim_gains = True) #climate scenario test
  results = model.simulate() # takes a while... save results

# calibration points (lists of pandas series)
# results = pd.read_csv('orca/data/results.csv', index_col=0, parse_dates=True)
  results['Combined_pump'] = results['DEL_HRO_pump'] + results['DEL_TRP_pump']
  sim = [results['DEL_HRO_pump'] / cfs_tafd,
       results['DEL_TRP_pump'] / cfs_tafd, 
       # (results['DEL_HRO_pump'] + results['DEL_TRP_pump']) / cfs_tafd,
       results['Combined_pump'] / cfs_tafd,
       results['SHA_storage'], 
       results['SHA_out'] / cfs_tafd,
       results['FOL_storage'],
       results['FOL_out'] / cfs_tafd,
       results['ORO_storage'],
       results['ORO_out'] / cfs_tafd,
       results['DEL_in'] / cfs_tafd,
       results['DEL_out'] / cfs_tafd]

  calibr_pts = ['HRO_pump','TRP_pump','Combined_pump','SHA_storage','SHA_out','FOL_storage','FOL_out','ORO_storage','ORO_out','DeltaIn','DeltaOut']
  for f in ['D','W','M','AS-OCT']:
    for s,c in zip(sim,calibr_pts):
      plotter.plotting(s, freq=f)
      plt.savefig('orca/figs/%s_%s.pdf' % (f,c), dpi=150)
      plt.close()


if not scenario:
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

  calibr_pts = ['HRO_pump','TRP_pump','Combined_pump','SHA_storage','SHA_out','FOL_storage','FOL_out','ORO_storage','ORO_out','DeltaIn','DeltaOut']
  for f in ['D','W','M','AS-OCT']:
    for s,o,c in zip(sim,obs,calibr_pts):
      plotter.compare(s, o, freq=f)
      plt.savefig('orca/figs/%s_%s.pdf' % (f,c), dpi=150)
      plt.close()

# elif scenario:
  