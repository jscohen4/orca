import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from util import *

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
df = pd.read_csv('orca-data-HB.csv', index_col=0, parse_dates=True)

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





############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################



def find_running_WYI(df,zscore1, zscore2): #water year type was calculated in scraper in master. I think we should determine wyt in a different script, just to keep model.py more concise. df is orca-data.csv
  fnf_SAC = df.BND_fnf.values * cfsd_mafd #Sacramento River flow at Bend Bridge, which is south of Shasta (fnf is full natural flow)
  fnf_FEA = df.ORO_fnf.values * cfsd_mafd #Feather River flow, from Oroville
  fnf_YUB = df.YRS_fnf.values * cfsd_mafd #Yuba River flow near Smartville. The Yuba River drains to the the Feather River downstream of Oroville
  fnf_AME = df.FOL_fnf.values * cfsd_mafd #American River flow downstream of Folsom
  snow_SAC = df.SHA_snowpack.values  #culmulative snowpack values
  snow_FEA = df.ORO_snowpack.values
  snow_AME = df.FOL_snowpack.values
  num_days = len(df)
  print num_days
  d_array = np.zeros(num_days)
  dowy_array = np.zeros(num_days)
  m_array = np.zeros(num_days)
  y_array = np.zeros(num_days)
  for t in range(0,num_days):
    d_array[t] = int(df.index.dayofyear[t]) 
    dowy_array[t] = water_day(d_array[t])
    m_array[t] = int(df.index.month[t])
    y_array[t] = int(df.index.year[t])  
  num_years = int(y_array[num_days-1] - y_array[0])
  apr_jul_fnf = np.zeros(num_years) #length of vector is latest year in dataset - earliest year (#of years). Maybe there's a more consise way to get this that will save simulation time
  oct_mar_fnf = np.zeros(num_years)
  remaining_flow = np.zeros((365,num_years)) #365 x number of years
  total_snowpack = np.zeros((365,num_years)) #365 x number of years
  forecastWYI_expt = np.zeros(num_days) #expt vs. envr? (export vs. environmental?)   #wyi is water year index
  forecastWYI_envr = np.zeros(num_days)


  current_year = 0

  for t in range(0,num_days): 
    dowy = int(dowy_array[t-1])
    m = m_array[t-1]
    if dowy == 0: #try 0 if worse
      current_year = current_year + 1
      oct_mar_fnf[current_year-1] = 0.0
      apr_jul_fnf[current_year-1] = 0.0
    if dowy < 182:    # oct-mar     #adding to flows for the three sections of the year.
      oct_mar_fnf[current_year - 1] += fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1]
      total_snowpack[dowy][current_year-1] = snow_SAC[t-1] + snow_FEA[t-1] + snow_AME[t-1]
    elif dowy < 305:    #apr-jul
      apr_jul_fnf[current_year - 1] += fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1]
      total_snowpack[dowy][current_year-1] = snow_SAC[t-1] + snow_FEA[t-1] + snow_AME[t-1]
    else: #aug-sept
      total_snowpack[dowy][current_year-1] = snow_SAC[t-1] + snow_FEA[t-1] + snow_AME[t-1]
  current_year = 0
  complete_year = 0 

  for t in range(0,num_years):   #for calculating october-march remaining flow?
    d = d_array[t-1]  
    dowy = int(dowy_array[t-1])
    if dowy == 0:
      current_year = current_year + 1
      remaining_flow[dowy][current_year-1] = oct_mar_fnf[current_year-1]
    elif dowy < 182: #oct-mar
      remaining_flow[dowy][current_year-1] = remaining_flow[dowy-1][current_year-1] - (fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1])
    else: #apr-jul
      remaining_flow[dowy][current_year-1] = 0.0
    
    if dowy == 304: #stop at august
      complete_year = complete_year + 1   

  #this is where the regression comes in. 
  snowfall_std = np.zeros(365)
  earlyflow_std = np.zeros(365)
  earlyflow_mean = np.zeros(365)
  snowfallCoef = np.zeros((365,2))  
  pred_devi = np.zeros(complete_year)
  for x in range(1,365):
    flow_each_year = remaining_flow[x-1]##this days set of  flow values, one value for each year(X vector) 
    snow_each_year = total_snowpack[x-1] ##this days set of cumulative snowpack values (X vector)  one value for each year
    coef = np.polyfit(snow_each_year[0:(complete_year-1)], apr_jul_fnf[0:(complete_year-1)],1)
    snowfallCoef[x-1][0] = coef[0]
    snowfallCoef[x-1][1] = coef[1]
    for y in range(1,complete_year):
      pred_devi[y-1] = apr_jul_fnf[y-1] - coef[0]*snow_each_year[y-1] - coef[1]
  
    snowfall_std[x-1] = np.std(pred_devi)
    if x < 183: #oct-mar
      earlyflow_std[x-1] = np.std(flow_each_year)
    else: #rest of year (no values needed)
      earlyflow_std[x-1] = 0.0

    earlyflow_mean[x-1] = np.mean(flow_each_year)

  prevValue = 10.0
  for t in range(1,num_days):
    d = d_array[t-1]   
    dowy = int(dowy_array[t-1])
    if dowy == 0: #begin october
      current_mar_oct = 0.0
      current_apr_jul = 0.0
      if t > 1: #after october 1st
        prevValue = min(forecastWYI_envr[t-2],10.0) #need documentation on this value
    
    if dowy < 182: #oct-mar
      current_mar_oct += (fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1])  #add flow of the four rivers each day
    elif dowy < 305: #apr-jul
      current_apr_jul += (fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1])
    snowpackEstimate = snowfallCoef[dowy][0]*(snow_SAC[t-1] + snow_FEA[t-1] + snow_AME[t-1]) + snowfallCoef[dowy][1] #using regression for snow esitmate
    forecastWYI_envr[t-1] = 0.3 * prevValue + 0.3 * (current_mar_oct + earlyflow_std[dowy]*zscore1 + earlyflow_mean[dowy]) + 0.4 *( current_apr_jul) #+ max(snowpackEstimate + snowfall_std[dowy]*zscore1 - current_apr_jul,0.0) ) #look more into this, why do we have zscores
    forecastWYI_expt[t-1] = 0.3 * prevValue + 0.3 * (current_mar_oct + earlyflow_std[dowy]*zscore2 + earlyflow_mean[dowy]) + 0.4 *( current_apr_jul) #+ max(snowpackEstimate + snowfall_std[dowy]*zscore2 - current_apr_jul,0.0) )
  SRIforecast = forecastWYI_expt
  #self.shasta.SRIforecast = self.forecastWYI_expt
  #self.folsom.SRIforecast = self.forecastWYI_expt
  #self.oroville.SRIforecast = self.forecastWYI_expt
  forecastWYT = []
  for x in range(0,num_days):
    if forecastWYI_envr[x-1] <= 5.4:
      forecastWYT.append("C")
    elif forecastWYI_envr[x-1] <= 6.5:
      forecastWYT.append("D")
    elif forecastWYI_envr[x-1] <= 7.8:
      forecastWYT.append("BN")
    elif forecastWYI_envr[x-1] <= 9.2:
      forecastWYT.append("AN")
    else:
      forecastWYT.append("W")
  return forecastWYT
  #self.shasta.wyt = self.forecastWYT  
  #self.oroville.wyt = self.forecastWYT
  #self.folsom.wyt = self.forecastWYT
    
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


df['WYT'] = pd.Series(find_running_WYI(df, 0, 0),index=df.index)

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

df.to_csv('orca-data-HB.csv')


