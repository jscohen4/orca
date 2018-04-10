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

df = pd.read_csv('access1-0_rcp45_r1i1p1_input_data.csv', index_col=0, parse_dates=True) #climate projection datafile
df = df[(df.index > '1996-09-30')]
df['WY'] = pd.Series([water_year(d) for d in df.index], index=df.index)
df['DOWY'] = pd.Series([water_year_day(d) for d in df.index], index=df.index)

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

df['SHA_fci'] = rolling_fci(df['SHA_fnf'], k=0.95, start=100000)
df.SHA_fci.fillna(method='bfill', inplace=True)

df['ORO_fci'] = rolling_fci(df['ORO_pr'], k=0.97, start=0)
df.ORO_fci.fillna(method='bfill', inplace=True)

df['FOL_fci'] = rolling_fci(df['FOL_pr'], k=0.97, start=0)
df.ORO_fci.fillna(method='bfill', inplace=True)



##clean up snowpack data and resample monthly 
snow_ids = ['GOL_swe','CSL_swe','HYS_swe','SCN_swe','RBB_swe','CAP_swe','RBP_swe',
				'HMB_swe','FOR_swe','RTL_swe','GRZ_swe','SDF_swe','SLT_swe','MED_swe']
dfs = df[snow_ids] #working with only snow for these calculations
num = dfs._get_numeric_data()
num[num < 0 ] = np.NaN
#num[num > 150 ] = np.NaN#oroville,folsom,shast,new bullards
num[num > 150 ] = np.NaN
dfs = dfs.interpolate(method = 'linear')
dfs = dfs.resample('M').mean()
df = df.drop(df[snow_ids],axis = 1)
df = df.join(dfs).fillna(method = 'ffill') #snow stations now cleaned up and back in main datafile 

# df = df[(df.index > '1999-09-30')]#start at 2000 water year
df = df[(df.index > '1999-09-30')]#start at 2000 water year

#sum of stations for each basins
df['BND_swe'] = df[['GOL_swe','CSL_swe']].mean(axis=1)
df['ORO_swe'] = df[['HYS_swe', 'SCN_swe', 'RBB_swe', 'CAP_swe']].mean(axis = 1) #taking out RBP (for this time), also test taking out RBB later
df['YRS_swe'] = df[['HMB_swe', 'FOR_swe', 'RTL_swe', 'GRZ_swe']].mean(axis = 1)
df['FOL_swe'] = df[['SDF_swe','SLT_swe']].mean(axis = 1)



BND = (df['BND_fnf'].to_frame(name='inf'))
ORO = (df['ORO_fnf'].to_frame(name='inf'))
YRS = (df['YRS_fnf'].to_frame(name='inf'))
FOL = (df['FOL_fnf'].to_frame(name='inf'))

df.to_csv('orca-data-processed-climate.csv')



