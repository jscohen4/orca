import numpy as np 
import scipy.stats as sp
import pandas as pd
from util import *
from sklearn import linear_model
from write_json import modify

# calc WYT and 8RI. add columns to datafile from cdec_scraper.
# confirm against http://cdec.water.ca.gov/cgi-progs/iodir/WSIHIST
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000
pd.options.mode.chained_assignment = None  # default='warn'

water_year = lambda d: d.year+1 if d.dayofyear >= 274 else d.year
def water_year_day(d):  #gets day of water year, accoutnting for leap years
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

#historical net Deltagains
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

prev_year = 9.8 # WY 1999 was 9.8
for y,g in df.groupby('WY'):
  flow = (g[SR_pts] * cfsd_mafd).sum(axis=1)
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

# San Joaquin Water Year Type 
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

df['8RI'] = ((df[SR_pts + SJR_pts] * cfsd_mafd) #8 station river index
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

### evaporation  regression
res_ids = ['SHA','ORO','FOL']
for r in res_ids:
	dfe = df[['%s_storage'%r,'%s_evap'%r,'%s_tas'%r]] #evaporation datafile
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
	reg.fit(X.T , evap)
	# R2 = reg.score(X.T, evap)
	# print('R2: %f' %(R2))
	coeffs = reg.coef_
	intercept = reg.intercept_
	modify('evap_regression.json',"%s_evap_coeffs" %r, coeffs.tolist())
	modify('evap_regression.json',"%s_evap_int"%r, intercept)

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

df = df[(df.index > '1996-09-30')]#start at 1997 water year
#sum of stations for each basin
df['BND_swe'] = df[['GOL_swe','CSL_swe']].mean(axis=1)
df['ORO_swe'] = df[['HYS_swe', 'SCN_swe', 'RBB_swe', 'CAP_swe']].mean(axis = 1) #taking out RBP (for this time), also test taking out RBB later
df['YRS_swe'] = df[['KTL_swe', 'HMB_swe', 'FOR_swe', 'RTL_swe', 'GRZ_swe']].mean(axis = 1)
df['FOL_swe'] = df[['SDF_swe', 'SNM_swe', 'SLT_swe']].mean(axis = 1)

BND = (df['BND_fnf'].to_frame(name='inf'))
ORO = (df['ORO_fnf'].to_frame(name='inf'))
YRS = (df['YRS_fnf'].to_frame(name='inf'))
FOL = (df['FOL_fnf'].to_frame(name='inf'))

#gains regression with perfect foresight WYI
# ### delta gains calculations
dfg = df[['MIL_fnf','NML_fnf','YRS_fnf','TLG_fnf','MRC_fnf','MKM_fnf','NHG_fnf','netgains','SR_WYI']] #gains datafile
stations = ['MIL','NML','YRS','TLG','MRC','MKM','NHG']
# dfg = df[['MIL_fnf','NML_fnf','YRS_fnf','TLG_fnf','MRC_fnf','MKM_fnf','NHG_fnf','netgains','WYI_sim']] #gains datafile
# stations = ['MIL','NML','YRS','TLG','MRC','MKM','NHG']

for station in stations:
  dfg['%s_fnf' %station] = df['%s_fnf' %station].shift(2)
  # dfg['%s_rol' %station] = df['%s_fnf' %station].rolling(10).sum()
  # dfg['%s_prev' %station] = df['%s_fnf' %station].shift(3)
  # dfg['%s_prev2' %station] = df['%s_fnf' %station].shift(4)
dfg = dfg.drop(pd.Timestamp('2012-07-01'), axis = 0)
dfg = dfg.dropna()
month_arr = np.arange(1,13)
R2_arr = np.zeros(12)
coeffs = []
intercepts = []
for m in month_arr:
  dfm = dfg[dfg.index.month==m] #for each month
  gains = dfm.netgains.values
  WYI = dfm.SR_WYI.values
  X = np.vstack([WYI])
  for station in stations:
    V = np.vstack([dfm['%s_fnf' %station]])
    X = np.vstack([X,V])
  X = X.T
  reg = linear_model.LinearRegression()
  reg.fit(X, gains.T)
  coeffs.append(reg.coef_)
  intercepts.append(reg.intercept_)
i = 0
for c in coeffs:
  i += 1
  modify('gains_regression.json',"month_%s" %i, c.tolist())
modify('gains_regression.json',"intercepts", intercepts)

df['gains_sim'] = pd.Series(index=df.index)
for index, row in df.iterrows():
  m = index.month
  X=[]
  b = coeffs[m-1]
  e = intercepts[m-1]
  for station in stations:
    X.append(df.loc[index,'%s_fnf' %station])
    # X.append(df.loc[index,'%s_rol' %station])
    # X.append(df.loc[index,'%s_prev' %station]) 
    # X.append(df.loc[index,'%s_prev2' %station])
  X.append(df.loc[ index,'SR_WYI'])
  X = np.array(X)
  gains = (np.sum(X * b) + e) * cfs_tafd
  df.loc[index, 'gains_sim'] = gains
df['gains_sim'] = df.gains_sim.fillna(method = 'bfill') * cfs_tafd #fill in missing beggining values (because of rolling)
df['netgains'] = df.netgains.fillna(method = 'bfill') * cfs_tafd #fill in missing beggining values (because of rolling)
for index, row in df.iterrows():
  ix = index.month
  d = index.day
  if ix == 10:
    df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 35
  if ix == 11:
    df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 4.5
  if ix == 12:
      df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] *3.5
  if ix == 1:
    df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 1.4 
  if (ix == 2):
    df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 1.7
  if ix == 3:
      df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 1.2
  if ix == 4:
      df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] *0.4
  if ix == 5: 
    df.loc[index, 'gains_sim'] = (df.loc[index, 'gains_sim'] - 12- d*0.4)*0.5 -20
  if ix ==6:
    df.loc[index, 'gains_sim'] = (df.loc[index, 'gains_sim'] - 15)*0.5
  if ix ==7:
    df.loc[index, 'gains_sim'] = (df.loc[index, 'gains_sim']) * 3 -20
  if (ix == 8):
      df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * 0.2 + d*0.55 -12
  if ix == 9:
    df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim'] * -10 
  df.loc[index, 'gains_sim'] = df.loc[index, 'gains_sim']*0.9


df.to_csv('orca-data-processed.csv')