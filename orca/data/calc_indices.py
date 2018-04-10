import numpy as np 
import scipy.stats as sp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from util import *
import matplotlib.pyplot as plt
from itertools import * 
from sklearn import linear_model
from sklearn.linear_model import (
    LinearRegression, TheilSenRegressor, RANSACRegressor, HuberRegressor)
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.mlab as mlab
from write_json import modify

sns.set_style('whitegrid')
# pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# calc WYT and 8RI. add columns to datafile from cdec_scraper.
# confirm against http://cdec.water.ca.gov/cgi-progs/iodir/WSIHIST
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000

water_year = lambda d: d.year+1 if d.dayofyear >= 274 else d.year
def water_year_day(d): 
	if d.is_leap_year:
		if d.dayofyear >= 275:
			return d.dayofyear - 274
		elif d.dayofyear <= 274 and d.dayofyear >= 59:	
			return d.dayofyear + 92
		else:
			return d.dayofyear + 92
	elif not d.is_leap_year:
		if d.dayofyear >= 274:
			return d.dayofyear - 273
		else:
			return d.dayofyear + 92
winter = lambda y: (y.index.month >= 10) | (y.index.month <= 3)
summer = lambda y: (y.index.month >= 4) & (y.index.month <= 7)
SR_pts = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
SJR_pts = ['NML_fnf', 'TLG_fnf', 'MRC_fnf', 'MIL_fnf']

# don't change this data
df = pd.read_csv('cdec-data.csv', index_col=0, parse_dates=True)
df['WY'] = pd.Series([water_year(d) for d in df.index], index=df.index)
df['DOWY'] = pd.Series([water_year_day(d) for d in df.index], index=df.index)
# estimate delta inflow from this (ignores GCD and direct precip)
# df.drop(['DeltaIn', 'netgains'], axis=1, inplace=True)

#historical net gains- will use regression later
df['DeltaIn'] = df['DeltaOut'] + df['HRO_pump'] + df['TRP_pump']
df['netgains'] = (df.DeltaIn - 
                  df.SHA_out.shift(5) - 
                  df.ORO_out.shift(3) - 
                  df.FOL_out.shift(1))
df.netgains.fillna(method='bfill', inplace=True)

def water_month(m):
  return m - 9 if m >= 9 else m + 3


def WYI_to_WYT(WYI, thresholds, values):
  for t,v in zip(thresholds,values):
    if WYI > t:
      return v

# Sacramento Water Year Index (historical)
get_SR_WYI = lambda x,p: 0.3*x[winter(x)].sum() + 0.4*x[summer(x)].sum() + 0.3*p

df['SR_WYI'] = pd.Series(index=df.index)

prev_year = 9.8 # WY 1996, 1999 was 9.8
for y,g in df.groupby('WY'):
  flow = (g[SR_pts] * cfsd_mafd).sum(axis=1)
  # plt.plot(flow.cumsum().values)
  WYI = get_SR_WYI(flow, prev_year)
  df.loc[df.WY==y, 'SR_WYI'] = WYI
  prev_year = np.min((10.0,WYI))


df['SR_WYT'] = df.SR_WYI.apply(WYI_to_WYT,
                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
                               values=['W', 'AN', 'BN', 'D', 'C'])

df['SR_WYT_rolling'] = (df.SR_WYI
                          .rolling(120).mean()
                          .apply(WYI_to_WYT,
                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
                               values=['W', 'AN', 'BN', 'D', 'C']))

df['SR_WYT_rolling'].fillna(method='bfill', inplace=True)

# San Joaquin Water Year Type #only using historical now--- may end up predicting to help with delta regressions
thresholds = [3.8, 3.1, 2.5, 2.1, 0.0]
values = ['W', 'AN', 'BN', 'D', 'C']
prev_year = 4.12 # WY 1996, 3.59 in 1999
get_SJR_WYI = lambda x,p: 0.2*x[winter(x)].sum() + 0.6*x[summer(x)].sum() + 0.2*p

df['SJR_WYI'] = pd.Series(index=df.index)
df['SJR_WYT'] = pd.Series(index=df.index)

for y,g in df.groupby('WY'):
  flow = (g[SJR_pts] * cfsd_mafd).sum(axis=1)
  WYI = get_SJR_WYI(flow, prev_year)
  prev_year = np.min((4.5,WYI))  
  for t,v in zip(thresholds,values):
    if WYI > t:
      df.loc[df.WY==y, 'SJR_WYT'] = v
      break

df['8RI'] = ((df[SR_pts + SJR_pts] * cfsd_mafd)
             .sum(axis=1)
             .resample('M')
             .sum())

df['8RI'].fillna(method='bfill', inplace=True)

# flood control indices
def rolling_fci(inflow, k, start):
  T = len(inflow)
  x = np.zeros(T)
  for i,t in enumerate(inflow.index):
    if t.month==10 and t.day==1:
      x[i] = start # cfs, start each WY here
    else:
      x[i] = inflow[t] + k*x[i-1]

  return pd.Series(x, index=inflow.index)

df['SHA_fci'] = rolling_fci(df['SHA_in_fix'], k=0.95, start=100000)
df.SHA_fci.fillna(method='bfill', inplace=True)

df['ORO_fci'] = rolling_fci(df['ORO_precip'], k=0.97, start=0)
df.ORO_fci.fillna(method='bfill', inplace=True)

df['FOL_fci'] = rolling_fci(df['FOL_precip'], k=0.97, start=0)
df.ORO_fci.fillna(method='bfill', inplace=True)

# folsom is different
# FMD = 136.4 - df.FMD_storage
# FMD[FMD > 45] = 45
# UNV = 266.4 - df.UNV_storage
# UNV[UNV > 80] = 80
# HHL = 207.6 - df.HHL_storage
# HHL[HHL > 75] = 75
# df['FOL_fci'] = FMD + UNV + HHL

# # folsom is different
# FMD = 136.4 - df.FMD_storage
# FMD[FMD > 45] = 45
# UNV = 266.4 - df.UNV_storage
# UNV[UNV > 80] = 80
# HHL = 207.6 - df.HHL_storage
# HHL[HHL > 75] = 75
# df['FOL_fci'] = FMD + UNV + HHL

### evap regression
res_ids = ['SHA','ORO','FOL']
for r in res_ids:
	dfe = df[['%s_storage'%r,'%s_evap'%r,'%s_tas'%r]]
	dfe[['storage','evap','tas']] = dfe[['%s_storage'%r,'%s_evap'%r,'%s_tas'%r]]
	dfe = dfe.dropna(axis = 0)

	storage = dfe.storage.values
	evap = dfe.evap.values
	temp = dfe.tas.values
	storage2 = storage**2
	temp2 = temp**2
	X2Y2 = (temp**2)*(storage**2)
	XY = temp*storage
	reg = linear_model.LinearRegression()
	X = np.vstack([temp, storage,XY,temp2, storage2])
	# reg.fit(X.T , evap)
	# R2 = reg.score(X.T, evap)
	# print('R2: %f' %(R2))
	coeffs = reg.coef_
	intercept = reg.intercept_
	modify('evap_regression.json',"%s_evap_coeffs" %r, coeffs.tolist())
	modify('evap_regression.json',"%s_evap_int"%r, intercept)

	# print(intercept)
	# print('Coefficients: \n', reg.coef_)

##clean up snowpack data and resample monthly 
snow_ids = ['GOL_swe','CSL_swe','HYS_swe','SCN_swe','RBB_swe','CAP_swe','RBP_swe','KTL_swe',
				'HMB_swe','FOR_swe','RTL_swe','GRZ_swe','SDF_swe','SNM_swe','SLT_swe','MED_swe']
dfs = df[snow_ids] #working with only snow for these calculations
num = dfs._get_numeric_data()
num[num < 0 ] = np.NaN
#num[num > 150 ] = np.NaN#oroville,folsom,shast,new bullards
num[num > 150 ] = np.NaN
dfs = dfs.interpolate(method = 'linear')
dfs = dfs.resample('M').mean()
df = df.drop(df[snow_ids],axis = 1)
df = df.join(dfs).fillna(method = 'ffill') #snow stations now cleaned up and back in main datafile 

df = df[(df.index > '1996-09-30')]#start at 2000 water year
#sum of stations for each basins
df['BND_swe'] = df[['GOL_swe','CSL_swe']].mean(axis=1)
df['ORO_swe'] = df[['HYS_swe', 'SCN_swe', 'RBB_swe', 'CAP_swe']].mean(axis = 1) #taking out RBP (for this time), also test taking out RBB later
df['YRS_swe'] = df[['KTL_swe', 'HMB_swe', 'FOR_swe', 'RTL_swe', 'GRZ_swe']].mean(axis = 1)
df['FOL_swe'] = df[['SDF_swe', 'SNM_swe', 'SLT_swe']].mean(axis = 1)
# plt.plot(df[['BND_swe','ORO_swe','YRS_swe','FOL_swe']])
# plt.show()

# forcast_pts = ['BND','ORO','YRS','FOL']
# snow_basins = ['BND_swe', 'ORO_swe', 'YRS_swe','FOL_swe']
# fnfs = ['BND_fnf','ORO_fnf','YRS_fnf','FOL_inf']

BND = (df['BND_fnf'].to_frame(name='inf'))
ORO = (df['ORO_fnf'].to_frame(name='inf'))
YRS = (df['YRS_fnf'].to_frame(name='inf'))
FOL = (df['FOL_fnf'].to_frame(name='inf'))

df.to_csv('orca-data-processed.csv')

############### for plotting WYI
	# fig, axes = plt.subplots(4,2) 
	# fig.subplots_adjust(hspace = 1)
	# month_val = []
	# titles = ['Oct','Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May']
	# month_ind = [10,11,12,1,2,3,4,5]
	# for ax, m, title in zip(axes.flat, month_ind, titles):
	# 	for index, row in Qm.iterrows():
	# 		ix = index.month
	# 		if ix == m: 
	# 			month_val.append(row['aprjul_cumulative'])
	# 	month_val.sort()
	# 	data_hist, bins = np.histogram(month_val)
	# 	(mu, sigma) = sp.norm.fit(month_val)
	# 	mean = np.mean(month_val)
	# 	median = np.median(month_val)
	# 	ax.set_title(title)
	# 	pdf = sp.norm.pdf(month_val, mu, sigma)
	# 	shape, loc, scale = sp.lognorm.fit(month_val)
	# 	logpdf = sp.lognorm.pdf(month_val, shape, loc, scale)
	# 	ax.plot(month_val,pdf)
	# 	ax.plot(month_val,logpdf, 'r--')
	# 	#ax.hist(month_val, edgecolor='black', linewidth=1.2, normed = 0)
	# 	ax.hist(month_val, edgecolor='black', linewidth=1.2, normed = 1, bins = 14)#, bins = 10)
	# 	ax.axvline(mean, color = 'm')
	# 	ax.axvline(median, color = 'y')

	# 	month_val = []
 #    #################plotting wyi timeseries

	# plt.figure()
	# Qm['Sac_WYT'] = Qm.WYI.apply(WYI_to_WYT,
 #                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
 #                               values=['W', 'AN', 'BN', 'D', 'C'])
	# b120_wyi = pd.read_csv('wyi_data.csv', index_col=0, parse_dates=True)
	# Qm = Qm.join(b120_wyi.loc[:,['50%']])
	# plt.axhline(9.2, color = 'g', label = 'Water Year Type Thresholds')
	# plt.axhline(7.8, color = 'g')
	# plt.axhline(6.5, color = 'g')
	# plt.axhline(5.4, color = 'g')
	# plt.plot(Qm['WYI'], 'b', label = 'Simulated Forecasts')#.plot()
	# plt.plot(Qm['50%'], 'r', label = 'Historical Forecasts')#.plot() 
	# plt.title('Sacramento Valley Water Year Type Index Forecasting',size = 18)
	# plt.xlabel('Date',size=16)
	# plt.ylabel('Water Year Type Index (million acre-ft)',size=16)
	# plt.legend(frameon=True,fontsize = 12)
	# df = df.join(Qm.loc[:,['WYI','Sac_WYT']])
	# df = df.fillna(method = 'bfill')


	# #plt.plot(Qm['WYI'])

