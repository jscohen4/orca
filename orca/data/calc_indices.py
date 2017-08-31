import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')

# calc WYT and 8RI. add columns to datafile from cdec_scraper.
# confirm against http://cdec.water.ca.gov/cgi-progs/iodir/WSIHIST
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
water_year = lambda d: d.year+1 if d.dayofyear >= 274 else d.year
winter = lambda y: (y.index.month >= 10) | (y.index.month <= 3)
summer = lambda y: (y.index.month >= 4) & (y.index.month <= 7)
SR_pts = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
SJR_pts = ['NML_fnf', 'TLG_fnf', 'MRC_fnf', 'MIL_fnf']

# don't change this data
df = pd.read_csv('orca-data.csv', index_col=0, parse_dates=True)

df['WY'] = pd.Series([water_year(d) for d in df.index], index=df.index)

# estimate delta inflow from this (ignores GCD and direct precip)
# df.drop(['DeltaIn', 'netgains'], axis=1, inplace=True)
df['DeltaIn'] = df['DeltaOut'] + df['HRO_pump'] + df['TRP_pump']
df['netgains'] = (df.DeltaIn - 
                  df.SHA_out.shift(5) - 
                  df.ORO_out.shift(3) - 
                  df.FOL_out.shift(1))
df.netgains.fillna(method='bfill', inplace=True)

def WYI_to_WYT(WYI, thresholds, values):
  for t,v in zip(thresholds,values):
    if WYI > t:
      return v

# Sacramento Water Year Index
get_SR_WYI = lambda x,p: 0.3*x[winter(x)].sum() + 0.4*x[summer(x)].sum() + 0.3*p

df['SR_WYI'] = pd.Series(index=df.index)

prev_year = 9.8 # WY 1999
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

# df.SR_WYI.plot()
# df.SR_WYI.rolling(120).mean().plot()
# plt.show()

# San Joaquin Water Year Type
thresholds = [3.8, 3.1, 2.5, 2.1, 0.0]
values = ['W', 'AN', 'BN', 'D', 'C']
prev_year = 3.59 # WY 1999
get_SJR_WYI = lambda x,p: 0.2*x[winter(x)].sum() + 0.6*x[summer(x)].sum() + 0.2*p
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

# df.ORO_fci.plot()
# plt.show()

# folsom is different
FMD = 136.4 - df.FMD_storage
FMD[FMD > 45] = 45
UNV = 266.4 - df.UNV_storage
UNV[UNV > 80] = 80
HHL = 207.6 - df.HHL_storage
HHL[HHL > 75] = 75
df['FOL_fci'] = FMD + UNV + HHL

df.to_csv('orca-data.csv')




  # def find_release_func(self, data): #still will need to work on this for the purposes of speeding things up- that'll be done when merging to master branch
  #   ##this function is used to make forecasts when calculating available storage for export releases from reservoir
  #   ##using data from 1996 to 2016 (b/c data is available for all inputs needed), calculate total flows in oct-mar period and apr-jul period
  #   ##based on linear regression w/snowpack (apr-jul) and w/inflow (oct-mar)
  #   ##this function is called before simulation loop, and the linear regression coefficient & standard deviation of linear regresion residuals
  #   ##is used in the find_available_storage function
  #   ## data used is from release-cdf-data.csv
  #   self.cum_snow = data['%s_cdf_snow'% self.key].values##cumulative yearly snowpack on each dayin data set
  #   self.daily_inflow = data['%s_cdf_inf'% self.key].values##cumulative oct-mar inflow (Based on how the data look, I'm inclined to think this is daily inflow- it matches that data in the master branch)
  #   #time_of_year = 0;##1 is oct-mar, 2 is apr-jul, 3 is aug-sept
  #   current_year = 0;
  #   complete_year = 0;##complete year counts all the years that have complete oct-jul data (partial years not used for lin. regression)
  #   self.oct_mar_cum_inflows = np.zeros(data.index.year[len(data)-1]-data.index.year[0]) #culmulative yearly october- march inflows (to be calculated), one per year
  #   self.apr_jul_cum_inflows = np.zeros(data.index.year[len(data)-1]-data.index.year[0]) #culmulative yearly april-july inflows (to be calculated), one per year
  #   self.regression_ceoffs = np.zeros((365,4))##coefficients for linear regression: 2 for oct-mar, 2 for apr-jul
  #   self.snow_stds = np.zeros(365) #standard deviations for snowpack regressions
  #   self.flow_stds = np.zeros(182) #standard deviations for oct-march flow regressions
  #   self.cum_snow_matrix = np.zeros((365,(data.index.year[len(data)-1]-data.index.year[0]))) #will become matrix with cumulative snowpack for each date, structured by year (365 x # of years matrix)
  #   self.daily_inflow_matrix = np.zeros((365,(data.index.year[len(data)-1]-data.index.year[0])))#will become matrix with daily inflows for each oct-mar date, structured by year (365 x # of years matrix) -only oct-mar inflows are included
    
  #   for t in range(1,len(data)): #I'll see if I can make this more consise, although I think it's good for now
  #     d = int(data.index.dayofyear[t-1])
  #     dowy = water_day(d)
  #     m = int(data.index.month[t-1])
  #     y = int(data.index.year[t-1])
  #     day = int(data.index.day[t-1])
  #     if dowy == 1:
  #       current_year +=1
    
  #     if m == 10: #october- beggining of water year
  #       time_of_year = "oct-mar"
  #     elif m == 4: # april- start of summer
  #       time_of_year = "apr-jul"
  #     elif m == 8 and day == 1: # august- end of apr-july regression
  #       time_of_year = "aug-sept"
  #       complete_year += 1##if data exists through jul, counts as a 'complete year' for linear regression purposes
      
  #     if time_of_year == "oct-mar":
  #       self.oct_mar_cum_inflows[current_year-1] += self.daily_inflow[t-1] * cfs_tafd##total oct-mar inflow (one value per year - Y vector in lin regression)
  #     elif time_of_year == "apr-jul":
  #       self.apr_jul_cum_inflows[current_year-1] += self.daily_inflow[t-1] * cfs_tafd##total apr-jul inflow (one value per year - Y vector in lin regression)
      
  #   current_year = 0;
  #   for t in range(1,len(data)):
  #     d = int(data.index.dayofyear[t-1])
  #     dowy = water_day(d)
  #     m = int(data.index.month[t-1])
  #     y = int(data.index.year[t-1])
  #     #cum oct-mar inflows through each day(X vector in oct-mar lin regression - each day has a unique 21 value (year) vector giving us 365 seperate lin regressions)
  #     if dowy == 1:
  #       current_year += 1;
  #       self.daily_inflow_matrix[dowy-1][current_year-1] = self.daily_inflow[t-1] * cfs_tafd
  #     elif dowy < 182:
  #       self.daily_inflow_matrix[dowy-1][current_year-1] = self.daily_inflow_matrix[dowy-2][current_year-1] + self.daily_inflow[t-1] * cfs_tafd
      
  #     ##cum snowpack through each day (X vector in apr-jul lin regression - each day has a unique 21 value (year) vector giving us 365 sepearte lin regressions)    
  #     self.cum_snow_matrix[dowy-1][current_year-1] = self.cum_snow[t-1]
      
  #   for x in range(1,182):
  #     flow_each_year = self.daily_inflow_matrix[x-1]##this days set of cumulative flow values, one value for each year(X vector) #still not seeing why this is culmalative
  #     coef = np.polyfit(flow_each_year[0:(complete_year-1)],self.oct_mar_cum_inflows[0:(complete_year-1)],1) #coef is the set of two regression coeffiients for this year's flow regression
  #     self.regression_ceoffs[x-1][0] = coef[0]
  #     self.regression_ceoffs[x-1][1] = coef[1]
  #     pred_dev = np.zeros(complete_year)
  #     for y in range(1,complete_year):
  #       pred_dev[y-1] = self.oct_mar_cum_inflows[y-1] - coef[0]*flow_each_year[y-1] - coef[1]##how much was the linear regression off actual observations

  #     self.flow_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals 
  #     ##for conservative estimate, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
  #   for x in range(1,365):
  #     snow_each_year = self.cum_snow_matrix[x-1]##this days set of cumulative snowpack values (X vector)  one value for each year
  #     coef = np.polyfit(snow_each_year[0:(complete_year-1)],self.apr_jul_cum_inflows[0:(complete_year-1)],1) #coef is the set of two regression coeffiients for this year's snowpack regression
  #     self.regression_ceoffs[x-1][2] = coef[0]
  #     self.regression_ceoffs[x-1][3] = coef[1]
  #     pred_dev = np.zeros(complete_year)
  #     for y in range(1,complete_year):
  #       pred_dev[y-1] = self.apr_jul_cum_inflows[y-1] - coef[0]*snow_each_year[y-1] - coef[1]##how much was the linear regression off actual observations

  #     self.snow_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals 
  #     ##for conservative estimate, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
