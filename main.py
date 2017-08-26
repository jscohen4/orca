import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cord import *

model = Model('cord/data/cord-data.csv', sd='10-01-1999')
model.create_flow_cdf('cord/data/cord-release-cdf-data.csv', sd='10-01-1996')
model.read_delta_inflow('cord/data/cord-delta-data.csv', sd = '10-01-1996')
model.find_running_WYI(0, 0)
results = model.simulate() # takes a while... save results
results.to_csv('cord/data/results.csv')
results = pd.read_csv('cord/data/results.csv', index_col=0, parse_dates=True)

# calibration points (lists of pandas series)
sim = [results['DEL_HRO_pump'][0:6209] / cfs_tafd,
       results['DEL_TRP_pump'][0:6209] / cfs_tafd,
       results['SHA_storage'][0:6209],
       results['SHA_out'][0:6209] / cfs_tafd,
       results['FOL_storage'][0:6209],
       results['FOL_out'][0:6209] / cfs_tafd,
       results['ORO_storage'][0:6209],
       results['ORO_out'][0:6209] / cfs_tafd]

obs = [model.df['HRO_pump'][0:6209],
       model.df['TRP_pump'][0:6209],
       model.df['SHA_storage'][0:6209],
       model.df['SHA_out'][0:6209],
       model.df['FOL_storage'][0:6209],
       model.df['FOL_out'][0:6209],
       model.df['ORO_storage'][0:6209],
       model.df['ORO_out'][0:6209]]

# for f in ['W','AS-OCT']:
#   i = 0
#   for s,o in zip(sim,obs):
#     plotter.compare(s, o, freq=f)
#     plt.savefig('cord/figs/%s%d.png' % (f,i), dpi=150)
#     i += 1

plotter.Rsquares(sim,obs)
