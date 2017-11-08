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
water_year_day = lambda d: d.dayofyear - 274 if d.dayofyear >= 274 else d.dayofyear + 91 
winter = lambda y: (y.index.month >= 10) | (y.index.month <= 3)
summer = lambda y: (y.index.month >= 4) & (y.index.month <= 7)
SR_pts = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
SJR_pts = ['NML_fnf', 'TLG_fnf', 'MRC_fnf', 'MIL_fnf']

# don't change this data
df = pd.read_csv('orca-data-HB.csv', index_col=0, parse_dates=True)
df2 = pd.read_csv('orca-release-cdf-data.csv', index_col=0, parse_dates=True)
df['WY'] = pd.Series([water_year(d) for d in df.index], index=df.index)
df2['WY'] = pd.Series([water_year(d) for d in df2.index], index=df2.index)
df2['DOWY'] = pd.Series([water_year_day(d) for d in df2.index], index=df2.index)
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

def get_forecast_WYI(df, zscore1, zscore2): 
  flow_sites = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
  snow_sites = ['SHA_snowpack', 'ORO_snowpack', 'FOL_snowpack']

  Qm = (df[flow_sites].sum(axis=1)
        .resample('M').sum()
        .to_frame(name='flow') * cfsd_mafd)

  Qm['WY'] = df.WY.resample('M').first() # need this for grouping

  Qm['snow'] = df[snow_sites].sum(axis=1).resample('M').first()


  # def octmar_tot_period_flow(x):
  #   ix = (x.index.month >= 10) | (x.index.month <= 3)
  #   octmar = (x[ix]. sum() - x[ix].cumsum())
  #   ix = (x.index.month >= 4) & (x.index.month <= 7)
  #   x[ix] = octmar[5]
  #   aprjul = x[ix]
  #   return pd.concat([octmar,aprjul])

  Qm = Qm[Qm.WY != Qm.WY[-1]] 
  Qm['octmar_cumulative'] = (Qm.groupby('WY').flow
                            .apply(octmar_cumulative)
                            .reset_index(level=0)
                            .drop('WY', axis=1)) 
  
  Qm['aprjul_cumulative'] = (Qm.groupby('WY').flow
                             .apply(aprjul_cumulative)
                             .reset_index(level=0).drop('WY', axis=1))

  Qm['octmar_flow_to_date'] = (Qm.groupby('WY').flow
                             .apply(octmar_flow_to_date)
                             .reset_index(level=0).drop('WY', axis=1))

  Qm['aprjul_flow_to_date'] = (Qm.groupby('WY').flow
                             .apply(aprjul_flow_to_date)
                             .reset_index(level=0).drop('WY', axis=1)) 

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
    coeffs = np.polyfit(Qm.snow[ix],Qm.aprjul_cumulative[ix],1)
    aprjul_slopes[m-1] = coeffs[0]
    aprjul_intercepts[m-1] = coeffs[1]
    aprjul_means[m-1] = np.mean(Qm.aprjul_cumulative[ix]) 
    aprjul_stds[m-1] =  np.std(Qm.aprjul_cumulative[ix])
    
    flow_coeffs = np.polyfit(Qm.octmar_flow_to_date[ix], Qm.octmar_cumulative[ix], 1)
    octmar_slopes[m-1] = flow_coeffs[0]
    octmar_intercepts[m-1] = flow_coeffs[1]
    octmar_means[m-1] = np.mean(Qm.octmar_cumulative[ix])
    octmar_stds[m-1] = np.std(Qm.octmar_cumulative[ix])
  stats =  {'aprjul_slope': aprjul_slopes,
                    'aprjul_intercept': aprjul_intercepts, 
                    'aprjul_mean': aprjul_means,
                    'aprjul_std':  aprjul_stds, 
                    'octmar_mean': octmar_means,
                    'octmar_std':  octmar_stds, 
                    'octmar_intercept': octmar_intercepts, 
                    'octmar_slope': octmar_slopes}

  stats = pd.DataFrame(stats, columns = ['aprjul_slope',
                    'aprjul_intercept', 
                    'aprjul_mean',
                    'aprjul_std', 
                    'octmar_mean',
                    'octmar_std', 
                    'octmar_intercept', 
                    'octmar_slope'])
  for s in stats: 
  	Qm[s] = pd.Series(index=Qm.index)
  	for m in range(1,13):
  		Qm.loc[Qm.index.month==m, s] = stats[s][m-1]

  #Qm.octmar_cumulative = Qm.octmar_cumulative.shift(periods=1)
  #Qm.aprjul_cumulative = Qm.aprjul_cumulative.shift(periods=1)

  Qm = Qm.fillna(0)
  Qm['WYI'] = pd.Series(index=Qm.index)
  prev = 10
  for index, row in Qm.iterrows():
    ix = index.month
    if (ix == 10) | (ix == 11):
      WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'])\
      + 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + Qm.loc[index, 'aprjul_mean'])# === + 0.675*Qm.loc[index,'aprjul_std'])
    elif (ix == 12) | (ix <= 4):
      WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'])\
      + 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] +  Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept'])#+ 0.675*Qm.loc[index,'aprjul_std'])# +  Qm.loc[index,'aprjul_std'])#Qm.loc[index, 'aprjul_mean'])#  + 0.675*Qm.loc[index,'aprjul_std'])
    elif ix == 5: 
      WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'])\
      + 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept'])#Qm.loc[index, 'aprjul_mean'])# d + 0.675*Qm.loc[index,'aprjul_std'])
      prev = min(WYI, 10)
    if (ix == 9) | (ix == 8):
      WYI = np.NaN
    Qm.loc[index, 'WYI'] = WYI
  Qm.WYI = Qm.WYI.shift(periods=-1)
  # Qm = Qm.drop(['aprjul_slope',
  #         'aprjul_intercept', 
  #         'aprjul_mean',
  #         'aprjul_std', 
  #         'octmar_mean',
  #         'octmar_std', 
  #         'octmar_intercept', 
  #         'octmar_slope',
  #         'aprjul_slope',
  #         'aprjul_intercept', 
  #         'aprjul_mean',
  #         'aprjul_std', 
  #         'octmar_mean',
  #         'octmar_std', 
  #         'octmar_intercept', 
  #         'octmar_slope',
  #         'octmar_cumulative',
  #         'aprjul_cumulative',
  #         'octmar_flow_to_date',
  #         'aprjul_flow_to_date'], axis=1) 
  Qm['SV_WYT'] = Qm.WYI.apply(WYI_to_WYT,
                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
                               values=['W', 'AN', 'BN', 'D', 'C'])
  Qm = Qm.fillna(method = 'ffill')
  #df['WYI','Sac_WYT'] = pd.Series(Qm['WYI','Sac_WYT'])
  DWY = Qm[['WYI','SV_WYT']]


  #df = df.join(DWY,how='right',rsuffix = '_new')
  #print df
  # df = df.fillna(method = 'bfill')

 # print df

  #print Qm
# ################### plotting stats
#   fig, axes = plt.subplots(4,2) 
#   fig.subplots_adjust(hspace = 1)
#   month_val = []
#   titles = ['Oct','Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May']
#   month_ind = [10,11,12,1,2,3,4,5]
#   for ax, m, title in zip(axes.flat, month_ind, titles):
#     for index, row in Qm.iterrows():
#       ix = index.month
#       if ix == m: 
#         month_val.append(row['aprjul_cumulative'])
#     month_val.sort()
#     data_hist, bins = np.histogram(month_val)
#     (mu, sigma) = sp.norm.fit(month_val)
#     mean = np.mean(month_val)
#     median = np.median(month_val)
#     ax.set_title(title)
#     pdf = sp.norm.pdf(month_val, mu, sigma)
#     shape, loc, scale = sp.lognorm.fit(month_val)
#     logpdf = sp.lognorm.pdf(month_val, shape, loc, scale)
#     ax.plot(month_val,pdf)
#     ax.plot(month_val,logpdf, 'r--')
#     #ax.hist(month_val, edgecolor='black', linewidth=1.2, normed = 0)
#     ax.hist(month_val, edgecolor='black', linewidth=1.2, normed = 1, bins = 14)#, bins = 10)
#     ax.axvline(mean, color = 'm')
#     ax.axvline(median, color = 'y')

#     month_val = []
#     #################plotting wyi timeseries

#   plt.figure()
#   Qm['Sac_WYT'] = Qm.WYI.apply(WYI_to_WYT,
#                                thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
#                                values=['W', 'AN', 'BN', 'D', 'C'])
#   b120_wyi = pd.read_csv('wyi_data.csv', index_col=0, parse_dates=True)
#   print b120_wyi
#   Qm = Qm.join(b120_wyi.loc[:,['50%']])
#   #print b120_wyi
#   #print(Qm.to_string())
#   plt.axhline(9.2, color = 'g', label = 'Water Year Type Thresholds')
#   plt.axhline(7.8, color = 'g')
#   plt.axhline(6.5, color = 'g')
#   plt.axhline(5.4, color = 'g')
#   plt.plot(Qm['WYI'], 'b', label = 'Simulated Forecasts')#.plot()
#   plt.plot(Qm['50%'], 'r', label = 'Historical Forecasts')#.plot() 
#   plt.title('Sacramento Valley Water Year Type Index Forecasting',size = 18)
#   plt.xlabel('Date',size=16)
#   plt.ylabel('Water Year Type Index (million acre-ft)',size=16)
#   plt.legend(frameon=True,fontsize = 12)
#   df = df.fillna(method = 'bfill')
#   print Qm


#   #plt.plot(Qm['WYI'])
  #df = df.join(DWY.loc[:,['WYI','SV_WYT']],how='right',rsuffix = '_new')
  #df = 
  #df = df.fillna(method = 'bfill')
  df =  pd.concat([df, DWY], axis=1, join_axes=[df.index])

  df[['WYI','WYI','SV_WYT']] = df[['WYI','SV_WYT']].fillna(method = 'bfill')
  #plt.show()
  df.to_csv('orca-data-HB-wyi-generate.csv')
#get_forecast_WYI(df,0,0) #wyt



def res_snow_inflow_regression(df1, df2, zscore1, zscore2): 
  snow_cdfs = ['ORO_cdf_snow', 'SHA_cdf_snow', 'FOL_cdf_snow']
  res_ids = ['SHA','ORO','FOL']
  ORO_inf = (df2['ORO_cdf_inf'].to_frame(name='inf') * cfsd_mafd)
  SHA_inf = (df2['SHA_cdf_inf'].to_frame(name='inf') * cfsd_mafd)
  FOL_inf = (df2['FOL_cdf_inf'].to_frame(name='inf') * cfsd_mafd)

  ORO_inf['snow_cdf'] = df2.ORO_cdf_snow
  SHA_inf['snow_cdf'] = df2.SHA_cdf_snow
  FOL_inf['snow_cdf'] = df2.FOL_cdf_snow

  res_frames = [SHA_inf,ORO_inf,FOL_inf]

  def octmar_flow_to_date(x):
    ix = (x.index.month >= 10) | (x.index.month <= 3)
    octmar = x[ix].cumsum()
    ix = (x.index.month >= 4) & (x.index.month <= 9)
    x[ix] = np.NaN #octmar[181]
    aprjul = x[ix] 
    return pd.concat([octmar,aprjul])
  

  def aprjul_flow_to_date(x):
    ix = (x.index.month >= 10) | (x.index.month <= 3)
    x[ix] = 0
    octmar = x[ix]
    ix = (x.index.month >= 4) & (x.index.month <= 7)
    aprjul =  x[ix].cumsum()
    ix = (x.index.month == 8) | (x.index.month == 9)
    x[ix] = np.NaN
    augsept = x[ix]
    return pd.concat([octmar,aprjul,augsept])


  for r,snow,res_id in zip(res_frames,snow_cdfs,res_ids):
    
    r['WY'] = df2.WY
    r['DOWY'] = df2.DOWY
    r['snow_cfd'] = df2[snow]
    r = r[r.WY != r.WY[-1]] 

    # r['octmar_cumulative'] = (r.groupby('WY').inf
    #                           .apply(octmar_cumulative)
    #                           .reset_index(level=0)
    #                           .drop('WY', axis=1))     
    
    # r['aprjul_cumulative'] = (r.groupby('WY').inf
    #                            .apply(aprjul_cumulative)
    #                            .reset_index(level=0).drop('WY', axis=1))

    r['octmar_flow_to_date'] = (r.groupby('WY').inf
                               .apply(octmar_flow_to_date)
                               .reset_index(level=0).drop('WY', axis=1))

    r['aprjul_flow_to_date'] = (r.groupby('WY').inf
                               .apply(aprjul_flow_to_date)
                               .reset_index(level=0).drop('WY', axis=1)) 
    #print r
    r.octmar_flow_to_date.fillna(method='ffill', inplace=True)
    r.aprjul_flow_to_date.fillna(method='ffill', inplace=True)
    #print(r.to_string())

    aprjul_slopes = np.zeros(365)
    aprjul_intercepts = np.zeros(365)
    aprjul_means = np.zeros(365)
    aprjul_stds = np.zeros(365)
    octmar_means = np.zeros(365)
    octmar_stds = np.zeros(365)
    octmar_slopes = np.zeros(365)
    octmar_intercepts = np.zeros(365)

    for d in range(1,182): # calendar months
      #oct_index = (r.index.month == 10)
      ix = (r.DOWY == d-1)
      coeffs = np.polyfit(r.snow_cdf[ix],r.octmar_flow_to_date[ix],1)
      octmar_slopes[d-1] = coeffs[0]
      octmar_intercepts[d-1] = coeffs[1]
      octmar_means[d-1] = np.mean(r.octmar_flow_to_date[ix])
      octmar_stds[d-1] = np.std(r.octmar_flow_to_date[ix])
    
    for d in range(1,366):
      ix = (r.DOWY == d-1)
      coeffs = np.polyfit(r.snow_cdf[ix],r.aprjul_flow_to_date[ix],1)
      aprjul_slopes[d-1] = coeffs[0]
      aprjul_intercepts[d-1] = coeffs[1]
      aprjul_means[d-1] = np.mean(r.aprjul_flow_to_date[ix]) 
      aprjul_stds[d-1] = np.std(r.aprjul_flow_to_date[ix])
    
    stats =  {'%s_aprjul_slope'%res_id: aprjul_slopes,
                      '%s_aprjul_intercept'%res_id: aprjul_intercepts, 
                      '%s_aprjul_mean'%res_id: aprjul_means,
                      '%s_aprjul_std'%res_id:  aprjul_stds, 
                      '%s_octmar_mean'%res_id: octmar_means,
                      '%s_octmar_std'%res_id:  octmar_stds, 
                      '%s_octmar_intercept'%res_id: octmar_intercepts, 
                      '%s_octmar_slope'%res_id: octmar_slopes}

    stats = pd.DataFrame(stats, columns = ['%s_aprjul_slope'%res_id,
                      '%s_aprjul_intercept'%res_id, 
                      '%s_aprjul_mean'%res_id,
                      '%s_aprjul_std'%res_id, 
                      '%s_octmar_mean'%res_id,
                      '%s_octmar_std'%res_id, 
                      '%s_octmar_intercept'%res_id, 
                      '%s_octmar_slope'%res_id])

    for s in stats: 
      r[s] = pd.Series(index=r.index)
      for d in range(1,366):
        r.loc[r.DOWY==d-1, s] = stats[s][d-1]
    r.rename(columns = {'aprjul_flow_to_date':'%s_aprjul_flow_to_date'%res_id,'octmar_flow_to_date':'%s_octmar_flow_to_date'%res_id}, inplace=True)
    r.drop(['inf','snow_cdf','WY','DOWY','snow_cfd'], axis=1, inplace=True)
    df2 =  pd.concat([df2, r], axis=1, join_axes=[df2.index])
  df2.to_csv('orca-data-cdf-generated.csv')



def res_snow_inflow_regression_2(df1, df2, zscore1, zscore2): 
  snow_cdfs = ['ORO_cdf_snow', 'SHA_cdf_snow', 'FOL_cdf_snow']
  res_ids = ['SHA','ORO','FOL']
  ORO_inf = (df2['ORO_cdf_inf'].to_frame(name='inf') * cfsd_mafd)
  SHA_inf = (df2['SHA_cdf_inf'].to_frame(name='inf') * cfsd_mafd)
  FOL_inf = (df2['FOL_cdf_inf'].to_frame(name='inf') * cfsd_mafd)

  ORO_inf['snow_cdf'] = df2.ORO_cdf_snow
  SHA_inf['snow_cdf'] = df2.SHA_cdf_snow
  FOL_inf['snow_cdf'] = df2.FOL_cdf_snow

  res_frames = [SHA_inf,ORO_inf,FOL_inf]

  def flow_to_date(x):
    ix = (x.index.month >= 1) 
    cum_flow = x[ix].cumsum()
    return cum_flow

  def rem_flow(x):
    ix = (x.index.month >= 1)
    remaining_flow = (x[ix].sum() - x[ix].cumsum())
    return remaining_flow


  for r,snow,res_id in zip(res_frames,snow_cdfs,res_ids):
    
    r['WY'] = df2.WY
    r['DOWY'] = df2.DOWY
    r['snow_cfd'] = df2[snow]
    r = r[r.WY != r.WY[-1]]

    r['cum_flow_to_date'] = (r.groupby('WY').inf
                               .apply(flow_to_date))
    r['remaining_flow'] = (r.groupby('WY').inf
                                .apply(rem_flow))
    r.cum_flow_to_date.fillna(method='ffill', inplace=True)
    slopes = np.zeros(365)
    intercepts = np.zeros(365)
    means = np.zeros(365)
    stds = np.zeros(365)

    # for d in range(1,182): # calendar months
    #   #oct_index = (r.index.month == 10)
    #   ix = (r.DOWY == d-1)
    #   coeffs = np.polyfit(r.snow_cdf[ix],r.octmar_flow_to_date[ix],1)
    #   octmar_slopes[d-1] = coeffs[0]
    #   octmar_intercepts[d-1] = coeffs[1]
    #   octmar_means[d-1] = np.mean(r.octmar_flow_to_date[ix])
    #   octmar_stds[d-1] = np.std(r.octmar_flow_to_date[ix])
    
    for d in range(1,366):
      ix = (r.DOWY == d-1)
      coeffs = np.polyfit(r.snow_cdf[ix],r.remaining_flow[ix],1)
      slopes[d-1] = coeffs[0]
      intercepts[d-1] = coeffs[1]
      means[d-1] = np.mean(r.remaining_flow[ix]) 
      stds[d-1] = np.std(r.remaining_flow[ix])
    
    stats =  {'%s_slope'%res_id: slopes,
                      '%s_intercept'%res_id: intercepts, 
                      '%s_mean'%res_id: means,
                      '%s_std'%res_id:  stds}


    stats = pd.DataFrame(stats, columns = ['%s_slope'%res_id,
                      '%s_intercept'%res_id, 
                      '%s_mean'%res_id,
                      '%s_std'%res_id]) 

    for s in stats: 
      r[s] = pd.Series(index=r.index)
      for d in range(1,366):
        r.loc[r.DOWY==d-1, s] = stats[s][d-1]
    r.rename(columns = {'flow_to_date':'%sflow_to_date'%res_id}, inplace=True)
    r.drop(['inf','snow_cdf','WY','DOWY','snow_cfd'], axis=1, inplace=True)
    df2 =  pd.concat([df2, r], axis=1, join_axes=[df2.index])
  df2.to_csv('orca-data-cdf-generated-new.csv')

res_snow_inflow_regression_2(df, df2,0 , 0)
    #inf_snow.groupby(inf_snow.index.dayofyear).cumsum()

def res_snow_inflow_regression_3(df1, df2, zscore1, zscore2): 
  snow_cdfs = ['ORO_cdf_snow', 'SHA_cdf_snow', 'FOL_cdf_snow']
  res_ids = ['SHA','ORO','FOL']
  ORO_inf = (df2['ORO_cdf_inf'].resample('M').sum().to_frame(name='inf') * cfsd_mafd)
  SHA_inf = (df2['SHA_cdf_inf'].resample('M').sum().to_frame(name='inf') * cfsd_mafd)
  FOL_inf = (df2['FOL_cdf_inf'].resample('M').sum().to_frame(name='inf') * cfsd_mafd)
  ORO_inf['snow_cdf'] = df2.ORO_cdf_snow
  SHA_inf['snow_cdf'] = df2.SHA_cdf_snow
  FOL_inf['snow_cdf'] = df2.FOL_cdf_snow

  res_frames = [SHA_inf,ORO_inf,FOL_inf]

  def flow_to_date(x):
    ix = (x.index.month >= 1) 
    cum_flow = x[ix].cumsum()

    return cum_flow
  



  for r,snow,res_id in zip(res_frames,snow_cdfs,res_ids):
    
    r['WY'] = df2.WY
    r['DOWY'] = df2.DOWY
    r['snow_cfd'] = df2[snow]
    r = r[r.WY != r.WY[-1]] 

    r['cum_flow_to_date'] = (r.groupby('WY').inf
                               .apply(flow_to_date))

    r.cum_flow_to_date.fillna(method='ffill', inplace=True)

    slopes = np.zeros(12)
    intercepts = np.zeros(12)
    means = np.zeros(12)
    stds = np.zeros(12)

    # for d in range(1,182): # calendar months
    #   #oct_index = (r.index.month == 10)
    #   ix = (r.DOWY == d-1)
    #   coeffs = np.polyfit(r.snow_cdf[ix],r.octmar_flow_to_date[ix],1)
    #   octmar_slopes[d-1] = coeffs[0]
    #   octmar_intercepts[d-1] = coeffs[1]
    #   octmar_means[d-1] = np.mean(r.octmar_flow_to_date[ix])
    #   octmar_stds[d-1] = np.std(r.octmar_flow_to_date[ix])
    
    for m in range(1,13):
      ix = (r.index.month == m)
      coeffs = np.polyfit(r.snow_cdf[ix],r.cum_flow_to_date[ix],1)
      slopes[m-1] = coeffs[0]
      intercepts[m-1] = coeffs[1]
      means[m-1] = np.mean(r.cum_flow_to_date[ix]) 
      stds[m-1] = np.std(r.cum_flow_to_date[ix])
    
    stats =  {'%s_slope'%res_id: slopes,
                      '%s_intercept'%res_id: intercepts, 
                      '%s_mean'%res_id: means,
                      '%s_std'%res_id:  stds}


    stats = pd.DataFrame(stats, columns = ['%s_slope'%res_id,
                      '%s_intercept'%res_id, 
                      '%s_mean'%res_id,
                      '%s_std'%res_id]) 

    for s in stats: 
      r[s] = pd.Series(index=r.index)
      for m in range(2,14):
        r.loc[r.index.month==m-1, s] = stats[s][m-2]
    r.rename(columns = {'flow_to_date':'%sflow_to_date'%res_id}, inplace=True)
    r.drop(['inf','snow_cdf','WY','DOWY','snow_cfd'], axis=1, inplace=True)
    df2 =  pd.concat([df2, r], axis=1, join_axes=[df2.index])
    df2[['%s_slope'%res_id,
                      '%s_intercept'%res_id, 
                      '%s_mean'%res_id,
                      '%s_std'%res_id]] = df2[['%s_slope'%res_id,
                      '%s_intercept'%res_id, 
                      '%s_mean'%res_id,
                      '%s_std'%res_id]].fillna(method = 'ffill')

  df2.to_csv('orca-data-cdf-generated-monthly3.csv')

# res_snow_inflow_regression_3(df, df2,0 , 0)
  #       .to_frame(name='flow') * cfsd_mafd)
          
  # Qm = (df[flow_sites].sum(axis=1)
  #       .resample('M').sum()
  #       .to_frame(name='flow') * cfsd_mafd)








    #cum_snow = df2['%s_cdf_snow'% r].values##cumulative yearly snowpack on each dayin data set
    #daily_inflow = df2['%s_cdf_inf'% r].values##cumulative oct-mar inflow (Based on how the data look, I'm inclined to think this is daily inflow- it matches that data in the master branch) 
    #res_data.append(cum_snow)
    #res_data.append(daily_inflow)
  # print res_data
  # inf_snow = (df2[res_data])#.sum(axis=1)
# for r in res_ids: 

        # .resample('M').sum()
        # .to_frame(name='flow') * cfsd_mafd)
    # Qm['%s_cdf_snow'% r] = cum_snow
    # Qm['%s_cdf_inf'% r] = daily_inflow


#res_snow_inflow_regression(df,df2,0,0)
  # Qm = (df2[flow_sites].sum(axis=1)
  #       .resample('M').sum()
  #       .to_frame(name='flow') * cfsd_mafd)
  #   ##this function is used to make forecasts when calculating available storage for export releases from reservoir
    ##using data from 1996 to 2016 (b/c data is available for all inputs needed), calculate total flows in oct-mar period and apr-jul period
    ##based on linear regression w/snowpack (apr-jul) and w/inflow (oct-mar)
    ##this function is called before simulation loop, and the linear regression coefficient & standard deviation of linear regresion residuals
    ##is used in the find_available_storage function
    ## data used is from release-cdf-data.csv

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


