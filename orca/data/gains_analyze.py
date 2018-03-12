import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000
df = pd.read_csv('orca-data.csv', index_col = 0)
df.netgains = df.netgains*cfs_tafd 
# df.gains_sim = df.gains_sim*cfs_tafd
df.gains_sim = df.gains_sim.shift(periods = -75)*cfs_tafd

df[['netgains','gains_sim']].plot(legend = True)
plt.show()