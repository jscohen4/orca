import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000

# df = pd.read_csv('cdec-data.csv',parse_dates = True, index_col = 0)
# # for m in range(1,13):
# for m in [1,4,7,10]:
# 	print(m)
# 	# df = pd.read_csv('DEL_HRO_pump.csv',parse_dates = True, index_col = 0)
# 	# df2 = pd.read_csv('DEL_TRP_pump.csv',parse_dates = True, index_col = 0)
# 	dfp = df['HRO_pump'] + df['TRP_pump']
# 	# df = df + df2	
# 	dfp = dfp/cfs_tafd
# 	# dfp = dfp.loc[(df.index.month==m)]
# 	dfp = dfp.resample('QS-OCT').sum()
# 	print(dfp)
# 	dfp = dfp.loc[(dfp.index.month==m)]
# 	dfp = dfp.rolling(5).mean()
# 	plt.plot(dfp,label = '%s'%m)
# plt.legend()
# plt.grid()
# plt.show()

df = pd.read_csv('DEL_HRO_pump.csv',parse_dates = True, index_col = 0)
df2 = pd.read_csv('DEL_TRP_pump.csv',parse_dates = True, index_col = 0)
df = df + df2

# for m in range(1,13):
for m in [1,4,7,10]:
	print(m)
	df = pd.read_csv('DEL_HRO_pump.csv',parse_dates = True, index_col = 0)
	# df2 = pd.read_csv('DEL_TRP_pump.csv',parse_dates = True, index_col = 0)
	# df = df + df2	     
	df = df/cfs_tafd
	# df = df.loc[(df.index.month==m)]
	df = df.resample('QS-OCT').sum()
	df = df.loc[(df.index.month==m)]
	df = df.rolling(20).mean().mean(axis = 1)
	plt.plot(df, label = '%s'%m)
plt.legend()
plt.grid()
plt.show()
