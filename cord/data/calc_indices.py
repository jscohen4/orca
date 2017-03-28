import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')

# calc WYT and 8RI. add columns to datafile from cdec_scraper.
# confirm against http://cdec.water.ca.gov/cgi-progs/iodir/WSIHIST
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
winter = lambda y: (y.index.month >= 10) | (y.index.month <= 3)
summer = lambda y: (y.index.month >= 4) & (y.index.month <= 7)
SR_pts = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
SJR_pts = ['NML_fnf', 'TLG_fnf', 'MRC_fnf', 'MIL_fnf']

# don't change this data
df = pd.read_csv('cord-data.csv', index_col=0, parse_dates=True)

# Sacramento Water Year Type
thresholds = [9.2, 7.8, 6.5, 5.4, 0.0]
values = ['W', 'AN', 'BN', 'D', 'C']
prev_year = 9.8 # WY 1999

def get_SR_WYT(y,index=None):
  global prev_year
  WYI = 0.3*y[winter(y)].sum() + 0.4*y[summer(y)].sum() + 0.3*prev_year
  prev_year = np.min((10.0,WYI))
  for t,v in zip(thresholds,values):
    if WYI > t:
      return v

df['SR_WYT'] = ((df[SR_pts] * cfsd_mafd)
                .sum(axis=1)
                .resample('AS-OCT')
                .apply(get_SR_WYT))

# df['SR_WYT_rolling'] = ((df[SR_pts] * cfsd_mafd)
#                         .sum(axis=1)
#                         .rolling(365)
#                         .sum().clip(lower=None, upper=10))

df.SR_WYT.fillna(method='ffill', inplace=True)

# df.SR_WYT_rolling.fillna(method='ffill', inplace=True)
# df.SR_WYT_rolling.plot()
# plt.show()

# San Joaquin Water Year Type
thresholds = [3.8, 3.1, 2.5, 2.1, 0.0]
values = ['W', 'AN', 'BN', 'D', 'C']
prev_year = 3.59 # WY 1999

def get_SJR_WYT(y):
  global prev_year
  WYI = 0.2*y[winter(y)].sum() + 0.6*y[summer(y)].sum() + 0.2*prev_year
  prev_year = np.min((4.5,WYI))
  
  for t,v in zip(thresholds,values):
    if WYI > t:
      return v

df['SJR_WYT'] = ((df[SJR_pts] * cfsd_mafd)
                 .sum(axis=1)
                 .resample('AS-OCT')
                 .apply(get_SJR_WYT))

df.SJR_WYT.fillna(method='ffill', inplace=True)

df['8RI'] = ((df[SR_pts + SJR_pts] * cfsd_mafd)
             .sum(axis=1)
             .resample('M')
             .sum())

df['8RI'].fillna(method='bfill', inplace=True)

# flood control indices
def rolling_fci(inflow):
  T = len(inflow)
  x = np.zeros(T)
  for i,t in enumerate(inflow.index):
    if t.month==10 and t.day==1:
      x[i] = 100000 # cfs, start each WY here
    else:
      x[i] = inflow[t] + 0.95*x[i-1]

  return pd.Series(x, index=inflow.index)

df['SHA_fci'] = rolling_fci(df['SHA_in_fix'])
df.SHA_fci.fillna(method='bfill', inplace=True)

df.to_csv('cord-data.csv')
