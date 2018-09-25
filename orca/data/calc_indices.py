import numpy as np 
import scipy.stats as sp
import pandas as pd
from .util import *
from sklearn.neural_network import MLPRegressor
from sklearn import linear_model
from sklearn import tree
from sklearn import preprocessing
from sklearn import utils
from sklearn.datasets import load_iris
from .write_json import modify
import json
import matplotlib.pyplot as plt
import pickle

# import matplotlib.pyplot as plt
# calc WYT and 8RI. add columns to datafile from cdec_scraper.
# confirm against http://cdec.water.ca.gov/cgi-progs/iodir/WSIHIST
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000
pd.options.mode.chained_assignment = None  # default='warn'

water_year = lambda d: d.year+1 if d.dayofyear >= 274 else d.year
# Sacramento Water Year Index (historical)
get_SR_WYI = lambda x,p: 0.3*x[winter(x)].sum() + 0.4*x[summer(x)].sum() + 0.3*p

# def water_year_day(d):  #obtain day of water year, which begins on October 1st
#   if d.is_leap_year:
#     if d.dayofyear >= 275:
#       return d.dayofyear - 274
#     elif d.dayofyear <= 274 and d.dayofyear >= 59:  
#       return d.dayofyear + 92
#     else:
#       return d.dayofyear + 92
#   elif not d.is_leap_year:
#     if d.dayofyear >= 274:
#       return d.dayofyear - 273
#     else:
#       return d.dayofyear + 92

winter = lambda y: (y.index.month >= 10) | (y.index.month <= 3)
summer = lambda y: (y.index.month >= 4) & (y.index.month <= 7)

def water_month(m): #obtain month of water year
  return m - 9 if m >= 9 else m + 3


def WYI_to_WYT(WYI, thresholds, values): #converts water year index to water year type
  for t,v in zip(thresholds,values):
    if WYI > t:
      return v

def rolling_fci(inflow, k, start): #used to obtain reservoir flood control index based on inflow
  T = len(inflow)
  x = np.zeros(T)
  for i,t in enumerate(inflow.index):
    if t.month==10 and t.day==1:
      x[i] = start # cfs, start each WY here
    else:
      x[i] = inflow[t] + k*x[i-1]

  return pd.Series(x, index=inflow.index)

def process(df,evap_regr,gains_regr,inf_regr): #used for historical data processing
  SR_pts = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf'] # flows to calculate Sacramento Valley WYI
  SJR_pts = ['NML_fnf', 'TLG_fnf', 'MRC_fnf', 'MIL_fnf'] #flows to calculate San-Joaquin Valley WYI

  df['WY'] = pd.Series([water_year(d) for d in df.index], index=df.index)
  df['DOWY'] = pd.Series([water_year_day(d) for d in df.index], index=df.index) #day of water year

  #historical net Deltagains, with perfect hindsight in pumping 
  df['DeltaIn'] = df['DeltaOut'] + df['HRO_pump'] + df['TRP_pump']
  df['netgains'] = (df.DeltaIn -  
                    df.SHA_out.shift(5) - 
                    df.ORO_out.shift(3) - 
                    df.FOL_out.shift(1))
  df.netgains.fillna(method='bfill', inplace=True)

  df['SR_WYI'] = pd.Series(index=df.index)

  prev_year = 12.89 # WY 1999 was 9.8
  for y,g in df.groupby('WY'):
    flow = (g[SR_pts] * cfsd_mafd).sum(axis=1)
    WYI = get_SR_WYI(flow, prev_year)
    df.loc[df.WY==y, 'SR_WYI'] = WYI
    prev_year = np.min((10.0,WYI))


  df['SR_WYT'] = df.SR_WYI.apply(WYI_to_WYT,
                                 thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
                                 values=['W', 'AN', 'BN', 'D', 'C'])

  df['SR_WYT_rolling'] = (df.SR_WYI #updating WYT throughout the year, thus limiting perfect foresight
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
 
  
  ###transformed inflows
  res_ids = ['SHA','ORO','FOL']
  for r in res_ids:
    coeffs = []
    intercepts = []
    df['%s_in_tr'%r] = pd.Series()
    for m in range(1,13):
        dfm = df[(df.index.month == m)]
        fnf = dfm['%s_fnf'%r].values
        inf = dfm['%s_in'%r].values
        fit = np.polyfit(fnf,inf,1)

        in_tr = fit[1] + fnf*fit[0]
        dfm['%s_in_tr'%r] = in_tr
        df['%s_in_tr'%r] = df['%s_in_tr'%r].fillna(dfm['%s_in_tr'%r])
        coeffs.append(fit[0])
        intercepts.append(fit[1])
    modify(inf_regr,"%s_inf_coeffs" %r, coeffs)
    modify(inf_regr,"%s_inf_int"%r, intercepts)
  df['ORO_precip'] = df[['SVL_pr','FRD_pr','DAV_pr','SBY_pr','CHS_pr','BRS_pr','QCY_pr','DES_pr']].mean(axis=1)
  df['FOL_precip'] = df[['BLC_pr','GTW_pr','PFH_pr']].mean(axis=1)
  df = df.drop(df[['SVL_pr','FRD_pr','DAV_pr','SBY_pr','CHS_pr','BRS_pr','QCY_pr','DES_pr','BLC_pr','GTW_pr','PFH_pr']],axis = 1)

  # flood control indices
  df['SHA_fci'] = rolling_fci(df['SHA_in_tr'], k=0.95, start=100000)
  df.SHA_fci.fillna(method='bfill', inplace=True)

  df['ORO_fci'] = rolling_fci(df['ORO_precip'], k=0.97, start=0)
  df.ORO_fci.fillna(method='bfill', inplace=True)

  df['FOL_fci'] = rolling_fci(df['FOL_precip'], k=0.97, start=0)
  df.FOL_fci.fillna(method='bfill', inplace=True)
  df.loc['2014-01-14','ORO_tas'] = df.loc['2014-01-13','ORO_tas']
  df.loc['2012-06-21','ORO_tas'] = df.loc['2012-06-20','ORO_tas']
  ### evaporation  regression
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
    reg.fit(X.T, evap)
    # R2 = reg.score(X.T, evap)
    # print('R2: %f' %(R2))
    coeffs = reg.coef_
    intercept = reg.intercept_
    modify(evap_regr,"%s_evap_coeffs" %r, coeffs.tolist())
    modify(evap_regr,"%s_evap_int"%r, intercept)



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
  df['YRS_swe'] = df[['GOL_swe','CSL_swe']].mean(axis=1)
  df['FOL_swe'] = df[['HYS_swe', 'SCN_swe', 'RBB_swe', 'CAP_swe']].mean(axis = 1) #taking out RBP (for this time), also test taking out RBB later
  df['ORO_swe'] = df[['KTL_swe', 'HMB_swe', 'FOR_swe', 'RTL_swe', 'GRZ_swe']].mean(axis = 1)
  df['BND_swe'] = df[['SDF_swe', 'SNM_swe', 'SLT_swe']].mean(axis = 1)

  BND = (df['BND_fnf'].to_frame(name='inf'))
  ORO = (df['ORO_fnf'].to_frame(name='inf'))
  YRS = (df['YRS_fnf'].to_frame(name='inf'))
  FOL = (df['FOL_fnf'].to_frame(name='inf'))

  #gains regression with perfect foresight WYI
  # ### delta gains calculations
  dfg = df[['MIL_fnf','NML_fnf','YRS_fnf','TLG_fnf','MRC_fnf','MKM_fnf','NHG_fnf','netgains','SR_WYI']] #gains datafile
  stations = ['SHA','ORO','YRS','FOL','MIL','NML','YRS','TLG','MRC','MKM','NHG']
  stations_WYI = ['SHA','ORO','YRS','FOL','MIL','NML','YRS','TLG','MRC','MKM','NHG','WYI']
  
  dfg = dfg[dfg.index != '1998-10-09']
  dfg = dfg[dfg.index != '1997-01-03']
  dfg = dfg[dfg.index != '1997-01-04']
  dfg = dfg[dfg.index != '1997-01-05']
  dfg = dfg[dfg.index != '2012-07-01']

  gain_errors = ['1998-10-09','1997-01-03','1997-01-04','1997-01-05','2012-07-01']
  for g in gain_errors:
    df.loc[g, 'netgains'] = np.NaN 
  for station in stations:
    dfg['%s_fnf' %station] = df['%s_fnf' %station].shift(2)
    dfg['%s_rol10' %station] = df['%s_fnf' %station].rolling(10).mean()
    dfg['%s_rol20' %station] = df['%s_fnf' %station].rolling(20).sum()
    dfg['%s_rol30' %station] = df['%s_fnf' %station].rolling(30).sum()

  dfg = dfg.dropna()

  month_arr = np.arange(1,13)
  dfgs = pd.DataFrame(index = df.index) #gains_sim df
  dfgs['gains_sim'] = pd.Series()
  for mth in month_arr:
    dfm = dfg[dfg.index.month == mth] #for each month #hist
    gains = dfm.netgains.values#.rolling(3, center = True, min_periods=1).mean().values
    WYI = dfm.SR_WYI.values#rolling(3, center = True, min_periods=1).mean().values
    m = dfm.index.month.values
    d = dfm.index.day.values

    X = np.vstack([WYI])

    X = np.vstack([X,[m]])

    for station in stations:
      V = np.vstack([dfm['%s_fnf' %station]])
      V = np.vstack([dfm['%s_rol10' %station]])
      V = np.vstack([dfm['%s_rol20' %station]])
      V = np.vstack([dfm['%s_rol30' %station]])
      X = np.vstack([X,V])
    X = X.T
    if mth in [10,11,12,1]:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 1000, hidden_layer_sizes=(50,), random_state=3)
    elif mth == 2: 
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 100000, hidden_layer_sizes=(1000,), random_state=3)
    elif mth == 3:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 100000, hidden_layer_sizes=(1000,), random_state=3)
    elif mth == 4:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 10000, hidden_layer_sizes=(75,), random_state=3)
    elif mth == 5:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 10000, hidden_layer_sizes=(600,), random_state=3)
    elif mth == 6:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 10000, hidden_layer_sizes=(100,), random_state=3)
    elif mth == 7:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 10000, hidden_layer_sizes=(400,), random_state=3)
    elif mth == 8:
      clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 10000, hidden_layer_sizes=(200,), random_state=3)
    elif mth == 9:
       clf = MLPRegressor(solver='adam', alpha=1e-5,max_iter = 10000, hidden_layer_sizes=(85,), random_state=3)
    ann = clf.fit(X, gains.T)
    pickle.dump(ann, open("orca/data/gains-neural-nets/ann_%s.pkl"%mth, "wb" ))
    month_gains = clf.predict(X)
    dfm['gains_sim'] = month_gains
    dfgs['gains_sim_%s'%mth] = dfm['gains_sim']
    dfgs = dfgs.fillna(0)
    dfgs['gains_sim'] = dfgs['gains_sim'] + dfgs['gains_sim_%s'%mth]
    #to test:
    # if mth == 7:
    #   plt.plot(gains)
    #   plt.plot(month_gains)
    #   plt.show()
  df['gains_sim'] = dfgs.gains_sim.fillna(method = 'bfill') * cfs_tafd #fill in missing beggining values (because of rolling)
  df['netgains'] = df.netgains.fillna(method = 'bfill') * cfs_tafd #fill in missing beggining values (because of rolling)
  # plt.plot(df.netgains)
  # plt.plot(df.gains_sim)  
  # plt.show()

  #if looping gains with historical. probably no longer needed with neural net
  df_g = pd.DataFrame() #gains datafile
  for WYT in ['C','D','BN','AN','W']:
    dfw = df[(df.SR_WYT == WYT)]
    means = dfw.netgains.groupby([dfw.index.strftime('%m-%d')]).mean()
    df_g[WYT] = means 
    days = [dfw.index.strftime('%Y-%m-%d')]
    days = days[0]

    for d in days:
      df.loc[df.index == d,'newgains'] = means[d[5:]]

  df['newgains'] = df.newgains.rolling(5).mean()
  df['newgains']=df.newgains.fillna(method='bfill')
  # df['gains_sim'] = (df['newgains']*0.75+df['gains_sim']*0.25)

  df['OMR'] = df.OMR + df.HRO_pump + df.TRP_pump

  df2 = df[(df.index > '2008-12-01')]#start at 1997 water year
  df['OMR_sim'] = pd.Series()
  df_OM = pd.DataFrame() #for means for OMR loops
  for WYT in ['C','D','BN','AN','W']:
    if WYT =='AN':
      dfWEY = df2[(df2.SR_WYT == 'W')]
      means_wet = dfWEY.OMR.groupby([dfWEY.index.strftime('%m-%d')]).mean()

      dfBN = df2[(df2.SR_WYT == 'BN')]
      means_BN = dfBN.OMR.groupby([dfBN.index.strftime('%m-%d')]).mean()
      means = means_wet.add(means_BN)/2
    else:
      dfw = df2[(df2.SR_WYT == WYT)]
      means = dfw.OMR.groupby([dfw.index.strftime('%m-%d')]).mean()
    df_OM[WYT] = means
    dfh = df[(df.SR_WYT == WYT)]

    days = [dfh.index.strftime('%Y-%m-%d')]
    days = days[0]
    for i,d in enumerate(days):
      if d[5:] != '02-29':
        df.loc[df.index == d,'OMR_sim'] = means[d[5:]]
      elif d[5:] == '02-29':
        df.loc[df.index == d,'OMR_sim'] = means[days[i-1][5:]]
  df['OMR_sim']=df.OMR_sim.fillna(method='ffill')
  return df, df_g, df_OM

def process_maurer(df,df_g,df_OMR,gains_regr,inf_regr,window): #used to process maurer projection data
  SR_pts = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
  SJR_pts = ['NML_fnf', 'TLG_fnf', 'MRC_fnf', 'MIL_fnf']
  df = df[(df.index > '1951-09-30')]
  gains_reg = json.load(open(gains_regr))
  inf_reg = json.load(open(inf_regr))

  df['WY'] = pd.Series([water_year(d) for d in df.index], index=df.index)
  df['DOWY'] = pd.Series([water_year_day(d) for d in df.index], index=df.index)
  #bias corrections for fnfs
  df['SHA_fnf'] = df.SHA_fnf+2800 
  df['BND_fnf'] = df.BND_fnf+2800 
  df['ORO_fnf'] = df.ORO_fnf + 1000
  df['FOL_fnf'] = df.FOL_fnf + 300
  df['YRS_fnf'] = df.YRS_fnf + 400
  df['MIL_fnf'] = df.MIL_fnf + 200
  df['MKM_fnf'] = df.MKM_fnf + 130
  df['TLG_fnf'] = df.TLG_fnf + 400
  df['NML_fnf'] = df.NML_fnf + 500
  
  snow_basins = ['BND_swe','FOL_swe','ORO_swe','YRS_swe']
  for sn in snow_basins:
    df[sn] = df[sn]/25.4
  # convert snow to inches

  #temp conversion and bias correction
  df['SHA_tas'] = (df['SHA_tas'] * 9/5 + 32)
  df['ORO_tas'] = (df['ORO_tas'] * 9/5 + 32)
  df['FOL_tas'] = (df['FOL_tas'] * 9/5 + 32)

  df['SHA_pr'] = (df['SHA_pr'] * 0.0393701)
  df['ORO_pr'] = (df['ORO_pr'] * 0.0393701)
  df['FOL_pr'] = (df['FOL_pr'] * 0.0393701)

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

  ###transfolrmed inflows
  res_ids = ['SHA','ORO','FOL']
  for r in res_ids:
    coeffs = inf_reg["%s_inf_coeffs" %r]
    intercepts = inf_reg["%s_inf_int"%r]
    df['%s_in_tr'%r] = pd.Series()

    for m in range(1,13):
      # for WYT in ['C','D','BN','AN','W']:
        dfm = df[(df.index.month == m)]
        # dfm = dfm[(dfm.SR_WYT == WYT)]
        inter = intercepts[m-1]
        coeff = coeffs[m-1]
        fnf = dfm['%s_fnf'%r].values
        in_tr = inter + fnf*coeff
        dfm['%s_in_tr'%r] = in_tr
        df['%s_in_tr'%r] = df['%s_in_tr'%r].fillna(dfm['%s_in_tr'%r])

########flood incides
  df['SHA_fci'] = rolling_fci(df['SHA_in_tr'], k=0.95, start=100000)
  df.SHA_fci.fillna(method='bfill', inplace=True)

  df['ORO_fci'] = rolling_fci(df['ORO_pr'], k=0.97, start=0)
  df.ORO_fci.fillna(method='bfill', inplace=True)

  df['FOL_fci'] = rolling_fci(df['FOL_pr'], k=0.97, start=0)
  df.ORO_fci.fillna(method='bfill', inplace=True)

  snow_ids = ['GOL_swe','CSL_swe','HYS_swe','SCN_swe','RBB_swe','CAP_swe','RBP_swe',
            'HMB_swe','FOR_swe','RTL_swe','GRZ_swe','SDF_swe','SLT_swe','MED_swe']


  dfs = df[snow_basins] #working with only snow for these calculations
  num = dfs._get_numeric_data()
  num[num < 0 ] = np.NaN
  #num[num > 150 ] = np.NaN#oroville,folsom,shast,new bullards
  num[num > 150 ] = np.NaN
  dfs = dfs.interpolate(method = 'linear')
  dfs = dfs.resample('M').mean()
  df = df.drop(df[snow_basins],axis = 1)
  df = df.join(dfs).fillna(method = 'ffill') #snow stations now cleaned up and back in main datafile 

  df = df[(df.index > '1951-09-30')]#start at 2000 water year

  df = df.drop(df[snow_ids],axis = 1)

  BND = (df['BND_fnf'].to_frame(name='inf'))
  ORO = (df['ORO_fnf'].to_frame(name='inf'))
  YRS = (df['YRS_fnf'].to_frame(name='inf'))
  FOL = (df['FOL_fnf'].to_frame(name='inf'))

  ## delta gains calculations
  dfg = df[['MIL_fnf','NML_fnf','YRS_fnf','TLG_fnf','MRC_fnf','MKM_fnf','NHG_fnf','SR_WYI']] #gains dataframe
  stations = ['SHA','ORO','YRS','FOL','MIL','NML','YRS','TLG','MRC','MKM','NHG']
  stations_WYI = ['SHA','ORO','YRS','FOL','MIL','NML','YRS','TLG','MRC','MKM','NHG','WYI']
  # dfg = df[['MIL_fnf','NML_fnf','YRS_fnf','TLG_fnf','MRC_fnf','MKM_fnf','NHG_fnf','netgains','WYI_sim']] #gains datafile
  # stations = ['MIL','NML','YRS','TLG','MRC','MKM','NHG']

  for station in stations:
    dfg['%s_fnf' %station] = df['%s_fnf' %station].shift(2)
    dfg['%s_rol10' %station] = df['%s_fnf' %station].rolling(10).sum()
    dfg['%s_rol20' %station] = df['%s_fnf' %station].rolling(20).sum()
    dfg['%s_rol30' %station] = df['%s_fnf' %station].rolling(30).sum()
  dfg = dfg.replace([np.inf, -np.inf], np.nan)
  dfg = dfg.dropna()

  # dfg = dfg.dropna()
  dfgs = pd.DataFrame(index = df.index) #gains_sim df
  dfgs['gains_sim'] = pd.Series()
  month_arr = np.arange(1,13)

  for mth in month_arr:
    dfm = dfg[dfg.index.month == mth] #for each month #hist
    WYI = dfm.SR_WYI.values#rolling(3, center = True, min_periods=1).mean().values
    m = dfm.index.month.values
    d = dfm.index.day.values

    X = np.vstack([WYI])

    X = np.vstack([X,[m]])

    for station in stations:
      V = np.vstack([dfm['%s_fnf' %station]])
      V = np.vstack([dfm['%s_rol10' %station]])
      V = np.vstack([dfm['%s_rol20' %station]])
      V = np.vstack([dfm['%s_rol30' %station]])
      X = np.vstack([X,V])
    X = X.T
    ann = pickle.load( open( "orca/data/gains-neural-nets/ann_%s.pkl"%mth, "rb" ))
    month_gains = ann.predict(X)
    if mth == 10:
      month_gains = month_gains - 500
    elif mth == 11:
      month_gains = month_gains - 200
    elif mth == 12:
      month_gains = month_gains *1.5
    elif mth == 1:
      month_gains = month_gains *1.4
    elif mth == 2:
      month_gains = month_gains *0.75
    elif mth == 4:
      month_gains = month_gains *1.5
    elif mth == 7:
      month_gains = month_gains - 20
    elif mth == 8:
      month_gains = (month_gains - 9000)*0.6
    elif mth == 9:
      month_gains = month_gains - 300
    #to test  
    # if mth == 7: 
    #   plt.plot(month_gains)

    dfm['gains_sim'] = month_gains
    dfgs['gains_sim_%s'%mth] = dfm['gains_sim']
    dfgs = dfgs.fillna(0)
    dfgs['gains_sim'] = dfgs['gains_sim'] + dfgs['gains_sim_%s'%mth]
  df['gains_sim'] = dfgs.gains_sim.fillna(method = 'bfill') * cfs_tafd #fill in missing beggining values (because of rolling)
  df['gains_sim'] = df.gains_sim.rolling(5, center = True, min_periods=3).mean()
    # plt.plot(df.gains_sim)  
    # plt.show()


    #get daily gains from regression
    #using coefficients obtained from process()
    
    # df['newgains'] = pd.Series()
    # for WYT in ['C','D','BN','AN','W']:
    #   dfw = df[(df.SR_WYT == WYT)]
    #   means = df_g[WYT] 
    #   days = [dfw.index.strftime('%Y-%m-%d')]
    #   days = days[0]

    #   for d in days:
    #     df.loc[df.index == d,'newgains'] = means[d[5:]]
    # df['newgains'] = df.newgains.rolling(5).mean()
    # df['newgains']=df.newgains.fillna(method='bfill')
    # df['gains_sim'] = (df['newgains']*0.75+df['gains_sim']*0.25)
  df['OMR_sim'] = pd.Series()
  for WYT in ['C','D','BN','AN','W']:
    dfh = df[(df.SR_WYT == WYT)]
    days = [dfh.index.strftime('%Y-%m-%d')]
    days = days[0]
    for i,d in enumerate(days):
      if d[5:] != '02-29':
        df.loc[df.index == d,'OMR_sim'] = df_OMR[WYT][d[5:]]
      elif d[5:] == '02-29':
        df.loc[df.index == d,'OMR_sim'] = df_OMR[WYT][days[i-1][5:]]
  df['OMR_sim'] = df.OMR_sim.rolling(5).mean()
  df['OMR_sim']=df.OMR_sim.fillna(method='bfill')
  return df
