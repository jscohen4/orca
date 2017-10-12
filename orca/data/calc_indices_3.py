import numpy as np 
import scipy.stats as sp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from util import *
import matplotlib.pyplot as plt
from itertools import * 
import matplotlib.mlab as mlab
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

def water_month(m):
  return m - 9 if m >= 9 else m + 3

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



def get_forecast_WYI(df, zscore1, zscore2): 
  flow_sites = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
  snow_sites = ['SHA_snowpack', 'ORO_snowpack', 'FOL_snowpack']

  Qm = (df[flow_sites].sum(axis=1)
         .resample('M').sum()
         .to_frame(name='flow') * cfsd_mafd)
  
  Qm['WY'] = df.WY.resample('M').first() # need this for grouping
  Qm[flow_sites] = (df[flow_sites].resample('M').sum() * cfsd_mafd)
  Qm['snow'] = df[snow_sites].sum(axis=1).resample('M').first()
  Qm[snow_sites] =  df[snow_sites].resample('M').first()

  def octmar_cumulative(x):
    ix = (x.index.month >= 10) | (x.index.month <= 3)
    octmar = (x[ix]. sum() - x[ix].cumsum())
    ix = (x.index.month >= 4) & (x.index.month <= 7)
    x[ix] = octmar[5]
    aprjul = x[ix]
    return pd.concat([octmar,aprjul])

  def aprjul_cumulative(x):
    ix = (x.index.month >= 4) & (x.index.month <= 7)
    aprjul = (x[ix].sum() - x[ix].cumsum())
    ix = (x.index.month >= 10) | (x.index.month <= 3)
    x[ix] = aprjul[0]
    octmar = x[ix]
    return pd.concat([octmar,aprjul]) 

  def octmar_flow_to_date(x):
    ix = (x.index.month >= 10) | (x.index.month <= 3)
    octmar = x[ix].cumsum()
    ix = (x.index.month >= 4) & (x.index.month <= 7)
    x[ix] = octmar[5]
    aprjul =  x[ix]
    return pd.concat([octmar,aprjul])
  

  def aprjul_flow_to_date(x):
    ix = (x.index.month >= 4) & (x.index.month <= 7)
    aprjul =  x[ix].cumsum()
    return pd.concat([aprjul])

  # def octmar_tot_period_flow(x):
  #   ix = (x.index.month >= 10) | (x.index.month <= 3)
  #   octmar = (x[ix]. sum() - x[ix].cumsum())
  #   ix = (x.index.month >= 4) & (x.index.month <= 7)
  #   x[ix] = octmar[5]
  #   aprjul = x[ix]
  #   return pd.concat([octmar,aprjul])

  Qm = Qm[Qm.WY != Qm.WY[-1]] 
  Qm['octmar_cumulative'] = (Qm.groupby('WY').flow.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
  Qm['aprjul_cumulative'] = (Qm.groupby('WY').flow.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
  Qm['octmar_flow_to_date'] = (Qm.groupby('WY').flow.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
  Qm['aprjul_flow_to_date'] = (Qm.groupby('WY').flow.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1))

  # Qm['BND_fnf_octmar_cumulative'] = (Qm.groupby('WY').BND_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['BND_fnf_aprjul_cumulative'] = (Qm.groupby('WY').BND_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['BND_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').BND_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
  # Qm['BND_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').BND_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1))
  # Qm['ORO_fnf_octmar_cumulative'] = (Qm.groupby('WY').ORO_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['ORO_fnf_aprjul_cumulative'] = (Qm.groupby('WY').ORO_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['ORO_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').ORO_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
  # Qm['ORO_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').ORO_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1)) 
  # Qm['YRS_fnf_octmar_cumulative'] = (Qm.groupby('WY').YRS_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['YRS_fnf_aprjul_cumulative'] = (Qm.groupby('WY').YRS_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['YRS_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').YRS_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
  # Qm['YRS_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').YRS_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1)) 
  # Qm['FOL_fnf_octmar_cumulative'] = (Qm.groupby('WY').FOL_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['FOL_fnf_aprjul_cumulative'] = (Qm.groupby('WY').FOL_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
  # Qm['FOL_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').FOL_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
  # Qm['FOL_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').FOL_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1)) 

  for river, snowpack in zip(flow_sites,['SHA_snowpack', 'ORO_snowpack', 'ORO_snowpack','FOL_snowpack']):
    if river == 'BND_fnf':
      Qm['BND_fnf_octmar_cumulative'] = (Qm.groupby('WY').BND_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['BND_fnf_aprjul_cumulative'] = (Qm.groupby('WY').BND_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['BND_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').BND_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
      Qm['BND_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').BND_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1))
    elif river == 'ORO_fnf':
      Qm['ORO_fnf_octmar_cumulative'] = (Qm.groupby('WY').ORO_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['ORO_fnf_aprjul_cumulative'] = (Qm.groupby('WY').ORO_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['ORO_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').ORO_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
      Qm['ORO_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').ORO_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1)) 
    elif river == 'YRS_fnf':
      Qm['YRS_fnf_octmar_cumulative'] = (Qm.groupby('WY').YRS_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['YRS_fnf_aprjul_cumulative'] = (Qm.groupby('WY').YRS_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['YRS_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').YRS_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
      Qm['YRS_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').YRS_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1))
    elif river == 'FOL_fnf':
      Qm['FOL_fnf_octmar_cumulative'] = (Qm.groupby('WY').FOL_fnf.apply(octmar_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['FOL_fnf_aprjul_cumulative'] = (Qm.groupby('WY').FOL_fnf.apply(aprjul_cumulative).reset_index(level=0).drop('WY', axis=1))
      Qm['FOL_fnf_octmar_flow_to_date'] = (Qm.groupby('WY').FOL_fnf.apply(octmar_flow_to_date).reset_index(level=0).drop('WY', axis=1))
      Qm['FOL_fnf_aprjul_flow_to_date'] = (Qm.groupby('WY').FOL_fnf.apply(aprjul_flow_to_date).reset_index(level=0).drop('WY', axis=1))
    aprjul_slopes = np.zeros(12)
    aprjul_intercepts = np.zeros(12)
    aprjul_means = np.zeros(12)
    aprjul_stds = np.zeros(12)
    octmar_means = np.zeros(12)
    octmar_stds = np.zeros(12)
    octmar_slopes = np.zeros(12)
    octmar_intercepts = np.zeros(12)

    for m in range(1,13): # calendar months
      oct_index = (Qm.index.month == 10)
      ix = (Qm.index.month == m)
      if river == 'BND_fnf':
        coeffs = np.polyfit(Qm.SHA_snowpack[ix],Qm.BND_fnf_aprjul_cumulative[ix],1)
        aprjul_means[m-1] = np.mean(Qm.BND_fnf_aprjul_cumulative[ix]) 
        aprjul_stds[m-1] =  np.std(Qm.BND_fnf_aprjul_cumulative[ix])
        octmar_means[m-1] = np.mean(Qm.BND_fnf_octmar_cumulative[ix])
        octmar_stds[m-1] = np.std(Qm.BND_fnf_octmar_cumulative[ix])

      elif river == 'ORO_fnf':
        coeffs = np.polyfit(Qm.ORO_snowpack[ix],Qm.ORO_fnf_aprjul_cumulative[ix],1)
        aprjul_means[m-1] = np.mean(Qm.ORO_fnf_aprjul_cumulative[ix]) 
        aprjul_stds[m-1] =  np.std(Qm.ORO_fnf_aprjul_cumulative[ix])
        octmar_means[m-1] = np.mean(Qm.ORO_fnf_octmar_cumulative[ix])
        octmar_stds[m-1] = np.std(Qm.ORO_fnf_octmar_cumulative[ix])
      elif river == 'YRS_fnf':
        coeffs = np.polyfit(Qm.ORO_snowpack[ix],Qm.YRS_fnf_aprjul_cumulative[ix],1)
        aprjul_means[m-1] = np.mean(Qm.YRS_fnf_aprjul_cumulative[ix]) 
        aprjul_stds[m-1] =  np.std(Qm.YRS_fnf_aprjul_cumulative[ix])
        octmar_means[m-1] = np.mean(Qm.YRS_fnf_octmar_cumulative[ix])
        octmar_stds[m-1] = np.std(Qm.YRS_fnf_octmar_cumulative[ix])
      elif river == 'FOL_fnf':
        coeffs = np.polyfit(Qm.FOL_snowpack[ix],Qm.FOL_fnf_aprjul_cumulative[ix],1)
        aprjul_means[m-1] = np.mean(Qm.FOL_fnf_aprjul_cumulative[ix]) 
        aprjul_stds[m-1] =  np.std(Qm.FOL_fnf_aprjul_cumulative[ix])
        octmar_means[m-1] = np.mean(Qm.FOL_fnf_octmar_cumulative[ix])
        octmar_stds[m-1] = np.std(Qm.FOL_fnf_octmar_cumulative[ix])
      #coeffs = np.polyfit(Qm.snow[ix],Qm.aprjul_cumulative[ix],1)
      aprjul_slopes[m-1] = coeffs[0]
      aprjul_intercepts[m-1] = coeffs[1]
      #aprjul_means[m-1] = np.mean(Qm.aprjul_cumulative[ix]) 
      #aprjul_stds[m-1] =  np.std(Qm.aprjul_cumulative[ix])
      
      #flow_coeffs = np.polyfit(Qm.octmar_flow_to_date[oct_index], Qm.octmar_cumulative[ix], 1)
      #octmar_slopes[m-1] = flow_coeffs[0]
      #octmar_intercepts[m-1] = flow_coeffs[1]
      #octmar_means[m-1] = np.mean(Qm.octmar_cumulative[ix])
      #octmar_stds[m-1] = np.std(Qm.octmar_cumulative[ix])
    stats =  {'%s_aprjul_slope'%river : aprjul_slopes,
                      '%s_aprjul_intercept'%river : aprjul_intercepts, 
                      '%s_aprjul_mean'%river : aprjul_means,
                      '%s_aprjul_std'%river :  aprjul_stds, 
                      '%s_octmar_mean'%river : octmar_means,
                      '%s_octmar_std'%river :  octmar_stds, 
                      '%s_octmar_intercept'%river : octmar_intercepts, 
                      '%s_octmar_slope'%river: octmar_slopes}

    stats = pd.DataFrame(stats, columns = ['%s_aprjul_slope'%river,
                      '%s_aprjul_intercept'%river, 
                      '%s_aprjul_mean'%river,
                      '%s_aprjul_std'%river, 
                      '%s_octmar_mean'%river,
                      '%s_octmar_std'%river, 
                      '%s_octmar_intercept'%river, 
                      '%s_octmar_slope'%river])
    for s in stats: 
    	Qm[s] = pd.Series(index=Qm.index)
    	for m in range(1,13):
    		Qm.loc[Qm.index.month==m, s] = stats[s][m-1]
    if river == 'BND_fnf':
      Qm.BND_fnf_octmar_cumulative = Qm.BND_fnf_octmar_cumulative.shift(periods=1)
      Qm.BND_fnf_aprjul_cumulative = Qm.BND_fnf_aprjul_cumulative.shift(periods=1)
    elif river == 'ORO_fnf':
      Qm.ORO_fnf_octmar_cumulative = Qm.ORO_fnf_octmar_cumulative.shift(periods=1)
      Qm.ORO_fnf_aprjul_cumulative = Qm.ORO_fnf_aprjul_cumulative.shift(periods=1)
    elif river == 'YRS_fnf':
      Qm.YRS_fnf_octmar_cumulative = Qm.YRS_fnf_octmar_cumulative.shift(periods=1)
      Qm.YRS_fnf_aprjul_cumulative = Qm.YRS_fnf_aprjul_cumulative.shift(periods=1)
    elif river == 'FOL_fnf':
      Qm.FOL_fnf_octmar_cumulative = Qm.FOL_fnf_octmar_cumulative.shift(periods=1)
      Qm.FOL_fnf_aprjul_cumulative = Qm.FOL_fnf_aprjul_cumulative.shift(periods=1)     

    Qm.octmar_cumulative = Qm.octmar_cumulative.shift(periods=1)
    Qm.aprjul_cumulative = Qm.aprjul_cumulative.shift(periods=1)
    Qm = Qm.fillna(0)
  print Qm
  Qm['WYI'] = pd.Series(index=Qm.index)
  prev = 10
  for index, row in Qm.iterrows():
    ix = index.month
    if (ix == 10) | (ix == 11):
      octmar_flow = 0
      for river in flow_sites: 
        octmar_flow += (Qm.loc[index, '%s_octmar_flow_to_date'%river] + Qm.loc[index, '%s_octmar_mean'%river])

      aprjul_flow = 0
      for river in flow_sites: 
        aprjul_flow += (Qm.loc[index, '%s_aprjul_flow_to_date'%river] + Qm.loc[index, '%s_aprjul_mean'%river])
      WYI = 0.3 * prev + 0.3 * (octmar_flow)\
      + 0.4 * (aprjul_flow)# === + 0.675*Qm.loc[index,'aprjul_std'])
    elif (ix == 12) | (ix <= 4):
      aprjul_flow = 0
      for river, snowpack in zip(flow_sites, ['SHA_snowpack', 'ORO_snowpack', 'ORO_snowpack','FOL_snowpack']):
        aprjul_flow += (Qm.loc[index, '%s_aprjul_flow_to_date'%river] +  Qm.loc[index, '%s_aprjul_slope'%river] * Qm.loc[index, '%s'%snowpack] + Qm.loc[index,'%s_aprjul_intercept'%river])
      WYI = 0.3 * prev + 0.3 * (octmar_flow)\
      + 0.4 * (aprjul_flow)#+ 0.675*Qm.loc[index,'aprjul_std'])# +  Qm.loc[index,'aprjul_std'])#Qm.loc[index, 'aprjul_mean'])#  + 0.675*Qm.loc[index,'aprjul_std'])
    elif ix == 5: 
      aprjul_flow = 0
      for river, snowpack in zip(flow_sites, ['SHA_snowpack', 'ORO_snowpack', 'ORO_snowpack','FOL_snowpack']):
        aprjul_flow += (Qm.loc[index, '%s_aprjul_flow_to_date'%river] +  Qm.loc[index, '%s_aprjul_slope'%river] * Qm.loc[index, '%s'%snowpack] + Qm.loc[index,'%s_aprjul_intercept'%river])
      WYI = 0.3 * prev + 0.3 * (octmar_flow)\
      + 0.4 * (aprjul_flow)#Qm.loc[index, 'aprjul_mean'])# d + 0.675*Qm.loc[index,'aprjul_std'])
      prev = min(WYI, 10)
    if (ix == 9) | (ix == 8):
      WYI = np.NaN
    Qm.loc[index, 'WYI'] = WYI
  Qm.WYI = Qm.WYI.shift(periods=-1)


################### plotting stats
  fig, axes = plt.subplots(4,2) 
  fig.subplots_adjust(hspace = 1)
  month_val = []
  titles = ['Oct','Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May']
  month_ind = [10,11,12,1,2,3,4,5]
  for ax, m, title in zip(axes.flat, month_ind, titles):
    for index, row in Qm.iterrows():
      ix = index.month
      if ix == m: 
        month_val.append(row['FOL_fnf'])
    month_val.sort()
    data_hist, bins = np.histogram(month_val)
    (mu, sigma) = sp.norm.fit(month_val)
    mean = np.mean(month_val)
    median = np.median(month_val)
    ax.set_title(title)
    pdf = sp.norm.pdf(month_val, mu, sigma)
    shape, loc, scale = sp.lognorm.fit(month_val)
    logpdf = sp.lognorm.pdf(month_val, shape, loc, scale)
    ax.plot(month_val,pdf)
    ax.plot(month_val,logpdf, 'r--')
    #ax.hist(month_val, edgecolor='black', linewidth=1.2, normed = 0)
    ax.hist(month_val, edgecolor='black', linewidth=1.2, normed = 1, bins = 14)#, bins = 10)
    ax.axvline(mean, color = 'm')
    ax.axvline(median, color = 'y')

    month_val = []
    #################plotting wyi timeseries

  plt.figure()
  Qm['Sac_WYT'] = Qm.WYI.apply(WYI_to_WYT,
                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
                               values=['W', 'AN', 'BN', 'D', 'C'])
  b120_wyi = pd.read_csv('wyi_data.csv', index_col=0, parse_dates=True)
  Qm = Qm.join(b120_wyi.loc[:,['50%']])
  #print b120_wyi
  #print(Qm.to_string())
  plt.axhline(9.2, color = 'g')
  plt.axhline(7.8, color = 'g')
  plt.axhline(6.5, color = 'g')
  plt.axhline(5.4, color = 'g')
  plt.plot(Qm['WYI'], 'b-')#.plot()
  plt.plot(Qm['50%'], 'r--')#.plot() 

  df = df.join(Qm.loc[:,['WYI','Sac_WYT']])
  df = df.fillna(method = 'bfill')
  #print Qm


  #plt.plot(Qm['WYI'])
  
  plt.show()
get_forecast_WYI(df,0,0) #wyt

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

  	 
  # #this is where the regression comes in. 
  # snowfall_std = np.zeros(365)
  # earlyflow_std = np.zeros(365)
  # earlyflow_mean = np.zeros(365)
  # snowfallCoef = np.zeros((365,2))  
  # pred_devi = np.zeros(complete_year)
  # for x in range(1,365):
  #   flow_each_year = remaining_flow[x-1]##this days set of  flow values, one value for each year(X vector) 
  #   snow_each_year = monthly_snowpack[x-1] ##this days set of cumulative snowpack values (X vector)  one value for each year
  #   coef = np.polyfit(snow_each_year[0:(complete_year-1)], apr_jul_fnf[0:(complete_year-1)],1)
  #   snowfallCoef[x-1][0] = coef[0]
  #   snowfallCoef[x-1][1] = coef[1]
  #   for y in range(1,complete_year):
  #     pred_devi[y-1] = apr_jul_fnf[y-1] - coef[0]*snow_each_year[y-1] - coef[1]
  
  #   snowfall_std[x-1] = np.std(pred_devi)
  #   if x < 183: #oct-mar
  #     earlyflow_std[x-1] = np.std(flow_each_year)
  #   else: #rest of year (no values needed)
  #     earlyflow_std[x-1] = 0.0

  #   earlyflow_mean[x-1] = np.mean(flow_each_year)

  # prevValue = 10.0
  # for t in range(1,6320):
  #   d = d_array[t-1]   
  #   dowy = int(dowy_array[t-1])
  #   if dowy == 0: #begin october
  #     current_mar_oct = 0.0
  #     current_apr_jul = 0.0
  #     if t > 1: #after october 1st
  #       prevValue = min(forecastWYI_envr[t-2],10.0) 
    
  #   if dowy < 182: #oct-mar
  #     current_mar_oct += (fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1])  #add flow of the four rivers each day
  #   elif dowy < 305: #apr-jul
  #     current_apr_jul += (fnf_SAC[t-1] + fnf_FEA[t-1] + fnf_YUB[t-1] + fnf_AME[t-1])
  #   snowpackEstimate = snowfallCoef[dowy][0]*(snow_SAC[t-1] + snow_FEA[t-1] + snow_AME[t-1]) + snowfallCoef[dowy][1] #using regression for snow esitmate
  #   forecastWYI_envr[t-1] = 0.3 * prevValue + 0.3 * (current_mar_oct + earlyflow_std[dowy]*zscore1 + earlyflow_mean[dowy]) + 0.4 *( current_apr_jul + max(snowpackEstimate + snowfall_std[dowy]*zscore1 - current_apr_jul,0.0) ) 
  #   forecastWYI_expt[t-1] = 0.3 * prevValue + 0.3 * (current_mar_oct + earlyflow_std[dowy]*zscore2 + earlyflow_mean[dowy]) + 0.4 *( current_apr_jul + max(snowpackEstimate + snowfall_std[dowy]*zscore2 - current_apr_jul,0.0) )
  # SRIforecast = forecastWYI_expt #sacramento river index
  # #self.shasta.SRIforecast = self.forecastWYI_expt
  # #self.folsom.SRIforecast = self.forecastWYI_expt
  # #self.oroville.SRIforecast = self.forecastWYI_expt
  # forecastWYT = []
  # for x in range(1,num_days+1):
  #   if forecastWYI_envr[x-1] <= 5.4:
  #     forecastWYT.append("C")
  #   elif forecastWYI_envr[x-1] <= 6.5:
  #     forecastWYT.append("D")
  #   elif forecastWYI_envr[x-1] <= 7.8:
  #     forecastWYT.append("BN")
  #   elif forecastWYI_envr[x-1] <= 9.2:
  #     forecastWYT.append("AN")
  #   else:
  #     forecastWYT.append("W")
  # return forecastWYT, SRIforecast
  #self.shasta.wyt = self.forecastWYT  
  #self.oroville.wyt = self.forecastWYT
  #self.folsom.wyt = self.forecastWYT
    
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
# df['SRI'] = pd.Series(find_running_WYI(df,0,0)
# df['WYT'] = pd.Series(find_running_WYI(df,0,0)[0],index=df.index) #wyt
# df['SRI'] = pd.Series(find_running_WYI(df,0,0)[1],index=df.index) #sacramento river index
# df['SHA_fci'] = rolling_fci(df['SHA_in_fix'], k=0.95, start=100000)
# df.SHA_fci.fillna(method='bfill', inplace=True)

# df['ORO_fci'] = rolling_fci(df['ORO_precip'], k=0.97, start=0)
# df.ORO_fci.fillna(method='bfill', inplace=True)

# # df.ORO_fci.plot()
# # plt.show()

# # folsom is different
# FMD = 136.4 - df.FMD_storage
# FMD[FMD > 45] = 45
# UNV = 266.4 - df.UNV_storage
# UNV[UNV > 80] = 80
# HHL = 207.6 - df.HHL_storage
# HHL[HHL > 75] = 75
# df['FOL_fci'] = FMD + UNV + HHL

# df.to_csv('orca-data-HB.csv')


