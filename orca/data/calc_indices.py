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
sns.set_style('whitegrid')

# calc WYT and 8RI. add columns to datafile from cdec_scraper.
# confirm against http://cdec.water.ca.gov/cgi-progs/iodir/WSIHIST
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000

water_year = lambda d: d.year+1 if d.dayofyear >= 274 else d.year
water_year_day = lambda d: d.dayofyear - 274 if d.dayofyear >= 274 else d.dayofyear + 91 
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
df.FOL_fci.fillna(method='bfill', inplace=True)

# def get_evap(tas,storage,coeffs,intercept):
#   T = len(storage)
#   y = np.zeros(T)
#   for i,(t,s) in enumerate(zip(tas,storage)):
#     X = [t,s,t*s,t**2,s**2]
#     y[i] = coeffs * X + intercept
#   return pd.Series(x, index=tas.index)
# print(df.SHA_tas)
# df['SHA_evap_sim'] = get_evap(df['SHA_tas'],df['SHA_storage'],[-6.12695687e+00,-5.56242336e-02,1.24250417e-03,5.85546703e-02,2.44535628e-06],175.162983004)
# df.SHA_fci.fillna(method='bfill', inplace=True)

# df['ORO_evap_sim'] = rolling_fci(df['ORO_precip'], k=0.97, start=0)
# df.ORO_fci.fillna(method='bfill', inplace=True)

# df['FOL_evap_sim'] = rolling_fci(df['FOL_precip'], k=0.97, start=0)
# df.ORO_fci.fillna(method='bfill', inplace=True)

# df.ORO_fci.plot()
# plt.show()

# df.SHA_fci.plot()
# plt.show()

# df.FOL_fci.plot()
# plt.show()
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

df = df[(df.index > '1999-09-30')]#start at 2000 water year
#sum of stations for each basins
df['BND_swe'] = df[['GOL_swe','CSL_swe']].mean(axis=1)
df['ORO_swe'] = df[['HYS_swe', 'SCN_swe', 'RBB_swe', 'CAP_swe']].mean(axis = 1) #taking out RBP (for this time), also test taking out RBB later
df['YRS_swe'] = df[['KTL_swe', 'HMB_swe', 'FOR_swe', 'RTL_swe', 'GRZ_swe']].mean(axis = 1)
df['FOL_swe'] = df[['SDF_swe', 'SNM_swe', 'SLT_swe']].mean(axis = 1)



# forcast_pts = ['BND','ORO','YRS','FOL']
# snow_basins = ['BND_swe', 'ORO_swe', 'YRS_swe','FOL_swe']
# fnfs = ['BND_fnf','ORO_fnf','YRS_fnf','FOL_inf']

BND = (df['BND_fnf'].to_frame(name='inf'))
ORO = (df['ORO_fnf'].to_frame(name='inf'))
YRS = (df['YRS_fnf'].to_frame(name='inf'))
FOL = (df['FOL_fnf'].to_frame(name='inf'))

# def octmar_cumulative(x):
# 	ix = (x.index.month >= 10) | (x.index.month <= 3)
# 	octmar = (x[ix]. sum() - x[ix].cumsum())
# 	ix = (x.index.month >= 4) & (x.index.month <= 7)
# 	x[ix] = octmar[5]
# 	aprjul = x[ix]
# 	return pd.concat([octmar,aprjul])

# def aprjul_cumulative(x):
# 	ix = (x.index.month >= 4) & (x.index.month <= 7)
# 	aprjul = (x[ix].sum() - x[ix].cumsum())
# 	ix = (x.index.month >= 10) | (x.index.month <= 3)
# 	x[ix] = aprjul[0]
# 	octmar = x[ix]
# 	ix = (x.index.month == 8) | (x.index.month ==9)
# 	x[ix] = aprjul[3]
# 	augsep = x[ix]

# 	return pd.concat([octmar,aprjul,augsep]) 

# def octmar_flow_to_date(x):
# 	ix = (x.index.month >= 10) | (x.index.month <= 3)
# 	octmar = x[ix].cumsum()
# 	ix = (x.index.month >= 4) & (x.index.month <= 7)
# 	x[ix] = octmar[5]
# 	aprjul =  x[ix]
# 	ix = (x.index.month == 8) | (x.index.month == 9)
# 	x[ix] = octmar[5]
# 	augsep = x[ix]
# 	return pd.concat([octmar,aprjul,augsep])

# def aprjul_flow_to_date(x):
# 	ix = (x.index.month >= 4) & (x.index.month <= 7)
# 	aprjul =  x[ix].cumsum()
# 	return pd.concat([aprjul])

# def get_forecast_WYI(df, zscore1, zscore2): 
# 	flow_sites = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
# 	snow_sites = ['BND_swe', 'ORO_swe', 'YRS_swe', 'FOL_swe']

# 	Qm = (df[flow_sites].sum(axis=1).resample('M').sum().to_frame(name='flow') * cfsd_mafd)

# 	Qm['WY'] = df.WY.resample('M').first() # need this for grouping

# 	Qm['snow'] = df[snow_sites].sum(axis=1).resample('M').first()
# 	Qm = Qm[Qm.WY != Qm.WY[-1]] 
# 	Qm['octmar_cumulative'] = (Qm.groupby('WY').flow
# 	                          .apply(octmar_cumulative)
# 	                          .reset_index(level=0)
# 	                          .drop('WY', axis=1)) 

# 	Qm['aprjul_cumulative'] = (Qm.groupby('WY').flow
# 	                           .apply(aprjul_cumulative))

# 	Qm['octmar_flow_to_date'] = (Qm.groupby('WY').flow
# 	                           .apply(octmar_flow_to_date))

# 	Qm['aprjul_flow_to_date'] = (Qm.groupby('WY').flow
# 	                           .apply(aprjul_flow_to_date)
# 	                           .reset_index(level=0).drop('WY', axis=1)) 

# 	aprjul_slopes = np.zeros(12)
# 	aprjul_intercepts = np.zeros(12)
# 	aprjul_means = np.zeros(12)
# 	aprjul_stds = np.zeros(12)
# 	octmar_means = np.zeros(12)
# 	octmar_stds = np.zeros(12)
# 	octmar_slopes = np.zeros(12)
# 	octmar_intercepts = np.zeros(12)

# 	for m in range(1,13): # calendar months
# 		oct_index = (Qm.index.month == 10)
# 		ix = (Qm.index.month == m)
# 		coeffs = np.polyfit(Qm.snow[ix],Qm.aprjul_cumulative[ix],1)
# 		aprjul_slopes[m-1] = coeffs[0]
# 		aprjul_intercepts[m-1] = coeffs[1]
# 		aprjul_means[m-1] = np.mean(Qm.aprjul_cumulative[ix]) 
# 		aprjul_stds[m-1] =  np.std(Qm.aprjul_cumulative[ix])
# 	for m in list(range(10,13))+list(range(1,3)):  
# 		flow_coeffs = np.polyfit(Qm.octmar_flow_to_date[oct_index], Qm.octmar_cumulative[ix], 1)
# 		octmar_slopes[m-1] = flow_coeffs[0]
# 		octmar_intercepts[m-1] = flow_coeffs[1]
# 		octmar_means[m-1] = np.mean(Qm.octmar_cumulative[ix])
# 		octmar_stds[m-1] = np.std(Qm.octmar_cumulative[ix])
# 	stats =  {'aprjul_slope': aprjul_slopes,
#                     'aprjul_intercept': aprjul_intercepts, 
#                     'aprjul_mean': aprjul_means,
#                     'aprjul_std':  aprjul_stds, 
#                     'octmar_mean': octmar_means,
#                     'octmar_std':  octmar_stds, 
#                     'octmar_intercept': octmar_intercepts, 
#                     'octmar_slope': octmar_slopes}

# 	stats = pd.DataFrame(stats, columns = ['aprjul_slope',
#                     'aprjul_intercept', 
#                     'aprjul_mean',
#                     'aprjul_std', 
#                     'octmar_mean',
#                     'octmar_std', 
#                     'octmar_intercept', 
#                     'octmar_slope'])
# 	for s in stats: 
# 		df[s] = pd.Series(index=df.index)
# 		for m in range(1,13):
# 			Qm.loc[Qm.index.month==m, s] = stats[s][m-1]
#   #Qm.octmar_cumulative = Qm.octmar_cumulative.shift(periods=1)
#   #Qm.aprjul_cumulative = Qm.aprjul_cumulative.shift(periods=1)

# 	Qm = Qm.fillna(0)
# 	Qm['WYI'] = pd.Series(index=Qm.index)
# 	prev = 10
# 	for index, row in Qm.iterrows():
# 		ix = index.month
# 		if (ix == 10) | (ix == 11):
# 			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'])\
# 			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + Qm.loc[index, 'aprjul_mean'])
# 		elif (ix == 12) | (ix <= 4):
# 			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'])\
# 			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] +  (Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept'])*zscore1)
# 		elif ix == 5: 
# 			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'])\
# 			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + (Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept'])*zscore1)
# #right now using zscore- more conservative- will discuss later
# 		if (ix == 9) | (ix == 8):
# 			WYI = np.NaN
# 		Qm.loc[index, 'WYI'] = WYI

# 	Qm.WYI = Qm.WYI.shift(periods=-1)
# 	# Qm.WYI = Qm.WYI.fillna(method = 'bfill')
# 	return(Qm.WYI)

# WYI_sim = get_forecast_WYI(df,0.4,0) #wyt
# df['WYI_sim'] = WYI_sim
# df.WYI_sim = df.WYI_sim.fillna(method = 'bfill')
# df.loc[df['WYI_sim'].isnull(),'WYI_sim'] = df['SR_WYI']

# df['WYT_sim'] = df.WYI_sim.apply(WYI_to_WYT,
#                                thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
#                                values=['W', 'AN', 'BN', 'D', 'C'])

# ### delta gains calculations
# dfg = df[['MIL_fnf','NML_fnf','YRS_fnf','TLG_fnf','MRC_fnf','MKM_fnf','NHG_fnf','netgains','WYI_sim']] #gains datafile
# stations = ['MIL','NML','YRS','TLG','MRC','MKM','NHG']
# for station in stations:
# 	dfg['%s_fnf' %station] = df['%s_fnf' %station]
# 	df['%s_rol' %station] = df['%s_fnf' %station].rolling(10).sum()
# 	dfg['%s_rol' %station] = df['%s_rol' %station]
# 	df['%s_prev' %station] = df['%s_fnf' %station].shift(1)
# 	dfg['%s_prev' %station] = df['%s_prev' %station]
# 	df['%s_prev2' %station] = df['%s_fnf' %station].shift(2)
# 	dfg['%s_prev2' %station] = df['%s_prev2' %station]
# dfg = dfg.drop(pd.Timestamp('2012-07-01'), axis = 0)
# dfg = dfg.dropna()
# month_arr = np.arange(1,13)
# R2_arr = np.zeros(12)
# coeffs = []
# intercepts = []
# for m in month_arr:
# 	dfm = dfg[dfg.index.month==m] #for each month
# 	gains = dfm.netgains.values
# 	WYI = dfm.WYI_sim.values
# 	X = np.vstack([WYI])
# 	for station in stations:
# 		V = np.vstack([dfm['%s_fnf' %station],dfm['%s_rol' %station].values, dfm['%s_prev' %station].values,dfm['%s_prev2' %station].values])
# 		X = np.vstack([X,V])
# 	X = X.T
# 	reg = linear_model.LinearRegression()
# 	reg.fit(X, gains.T)
# 	# R2 = reg.score(X, gains.T)
# 	# R2_arr[m-1] = R2
# 	# print('Coefficients: \n', reg.coef_)
# 	coeffs.append(reg.coef_)
# 	intercepts.append(reg.intercept_)

# df['gains_sim'] = pd.Series(index=df.index)
# for index, row in df.iterrows():
# 	m = index.month
# 	X=[]
# 	b = coeffs[m-1]
# 	e = intercepts[m-1]
# 	for station in stations:
# 		X.append(df.loc[index,'%s_fnf' %station])
# 		X.append(df.loc[index,'%s_rol' %station])
# 		X.append(df.loc[index,'%s_prev' %station]) 
# 		X.append(df.loc[index,'%s_prev2' %station])
# 	X.append(df.loc[index,'WYI_sim'])
# 	X = np.array(X)
# 	gains = np.sum(X * b)*cfs_tafd + e
# 	df.loc[index, 'gains_sim'] = gains
# df['gains_sim'] = df.gains_sim.fillna(method = 'bfill') #fill in missing beggining values (because of rolling)
# # df[['netgains','gains_sim']].plot(legend = True)
# # plt.show()

# ###making forcasts for reservoir carryover:
# def flow_to_date(x):
# 	ix = (x.index.month >= 1) 
# 	cum_flow = x[ix].cumsum()
# 	return cum_flow

# def rem_flow(x):
# 	ix = (x.index.month >= 1)
# 	remaining_flow = (x[ix].sum() - x[ix].cumsum())
# 	return remaining_flow

# snow_sites = ['BND_swe', 'ORO_swe','FOL_swe']
# res_ids = ['SHA','ORO','FOL']

# SHA_inf = (df['BND_fnf'].to_frame(name='inf') * cfsd_mafd)  
# ORO_inf = (df['ORO_fnf'].to_frame(name='inf') * cfsd_mafd)
# FOL_inf = (df['FOL_fnf'].to_frame(name='inf') * cfsd_mafd)


# res_frames = [SHA_inf,ORO_inf,FOL_inf]

# for r, swe, res_id in zip(res_frames, snow_sites, res_ids):

# 	r['WY'] = df.WY
# 	r['DOWY'] = df.DOWY
# 	r['snowpack'] = df[swe]
# 	r = r[r.WY != r.WY[-1]]

# 	r['cum_flow_to_date'] = (r.groupby('WY').inf
# 	                           .apply(flow_to_date))
# 	r['remaining_flow'] = (r.groupby('WY').inf
# 	                            .apply(rem_flow))
# 	r.cum_flow_to_date.fillna(method='ffill', inplace=True)
# 	slopes = np.zeros(365)
# 	intercepts = np.zeros(365)
# 	means = np.zeros(365)
# 	stds = np.zeros(365)
# 	for d in range(1,366):
# 		ix = (r.DOWY == d-1)
# 		coeffs = np.polyfit(r.snowpack[ix],r.remaining_flow[ix],1)
# 		slopes[d-1] = coeffs[0]
# 		intercepts[d-1] = coeffs[1]
# 		means[d-1] = np.mean(r.remaining_flow[ix]) 
# 		stds[d-1] = np.std(r.remaining_flow[ix])

# 	stats =  {'%s_slope'%res_id: slopes,
# 	                  '%s_intercept'%res_id: intercepts, 
# 	                  '%s_mean'%res_id: means,
# 	                  '%s_std'%res_id:  stds}


# 	stats = pd.DataFrame(stats, columns = ['%s_slope'%res_id,
# 	                  '%s_intercept'%res_id, 
# 	                  '%s_mean'%res_id,
# 	                  '%s_std'%res_id]) 
# 	for s in stats: 
# 		r[s] = pd.Series(index=r.index)
# 		for d in range(1,366):
# 			r.loc[r.DOWY==d-1, s] = stats[s][d-1]
# 	r.rename(columns = {'cum_flow_to_date':'%s_cum_flow_to_date'%res_id}, inplace=True)
# 	r.rename(columns = {'remaining_flow':'%s_remaining_flow'%res_id}, inplace=True)
# 	r.rename(columns = {'snowpack':'%s_snowpack'%res_id}, inplace=True)

# 	r.drop(['inf','WY','DOWY'], axis=1, inplace=True)
# 	df = pd.concat([df, r], axis=1, join_axes=[df.index])

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

