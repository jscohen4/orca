import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cord import *

datafile = 'cord/data/cord-data.csv'

df = pd.read_csv(datafile, index_col=0, parse_dates=True)
T = len(df)

shasta = Reservoir(df, 'SHA')
delta = Delta(df, 'DEL')

for t in range(1,T):
  sumin = shasta.Q[t] + df.FOL_in[t] * cfs_tafd + df.ORO_in[t] * cfs_tafd
  delta.calc_flow_bounds(t, sumin)
  shasta.step(t, delta.dmin[t], delta.sodd_cvp[t])
  # folsom.step(t, delta.dmin[t], delta.sodd_cvp[t])
  # oroville.step(t, delta.dmin[t], delta.sodd_swp[t])
  # delta.step(t, gains, shasta.R[t], folsom.R[t], oroville.R[t])

results = shasta.results_as_df(df.index)

# (results.SHA_rexport / cfs_tafd).plot()
# df.TRP_pump.plot()

# results = results.resample('M').mean()
# df = df.resample('M').mean()

plt.subplot(2,1,1)
h = results.SHA_storage.plot()
df.SHA_storage.plot(ax=h)

plt.subplot(2,1,2)
h = results.SHA_out.plot()
(df.SHA_out * cfs_tafd).plot(ax=h)

plt.show()

# plt.scatter(results.SHA_out.values, (df.SHA_out*cfs_tafd).values)
# plt.show()
