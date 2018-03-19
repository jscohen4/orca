import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000

df = pd.read_csv('orca-data-forecasted.csv', index_col=0, parse_dates=True)
df.netgains = df.netgains *cfs_tafd
df.gains_sim = df.gains_sim *cfs_tafd
for index, row in df.iterrows():
	ix = index.month
	# if (ix <= 6) & (ix >= 4):
	# 	df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim']*0.01
	# if (ix >= 4) & (ix <= 9):
	# 	df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 0.2
	if (ix >= 3) & (ix <= 9):
			df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 0.2
	if (ix >= 5) & (ix <= 8):
			df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] - 10
    # self.gains[t] = gains 
plt.plot(df.netgains)
plt.plot(df.gains_sim)
plt.show() 