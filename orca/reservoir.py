from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from util import *

0.007963636
class Reservoir():

  def __init__(self, df, key):
    T = len(df)
    self.dayofyear = df.index.dayofyear
    self.month = df.index.month
    self.key = key
    self.wyt = df.SV_WYT# 120 day MA lag
    for k,v in json.load(open('orca/data/%s_properties.json' % key)).items():
      setattr(self,k,v)

    self.Q = df['%s_in_fix'% key].values * cfs_tafd
    self.E = df['%s_evap'% key].values * cfs_tafd
    self.fci = df['%s_fci' % key].values
    self.cum_flow_to_date = df['%s_cum_flow_to_date' % key].values #already in tadf
    self.slope =  df['%s_slope' % key].values
    self.intercept = df['%s_intercept' % key].values
    self.mean = df['%s_mean' % key].values
    self.std = df['%s_std' % key].values  
    self.WYI = df['WYI'].values
    self.S = np.zeros(T)
    self.R = np.zeros(T)
    self.Rtarget = np.zeros(T)
    self.R_to_delta = np.zeros(T)
    self.S[0] = df['%s_storage' % key].iloc[0]
    self.R[0] = 0
    self.storage_bounds = np.zeros(2)
    self.index_bounds = np.zeros(2)
    self.tocs = np.zeros(T)

    #tocs rule variables
    self.tocs_index = []
    for i,v in enumerate(self.tocs_rule['index']):        
        self.tocs_index.append(np.zeros(366))
        for day in range(0, 366):  
            self.tocs_index[i][day] = np.interp(day, self.tocs_rule['dowy'][i], self.tocs_rule['storage'][i])

    self.nodds = np.zeros(367)
    for i in range(0,366):  
        self.nodds[i] = np.interp(i, first_of_month, self.nodd)


  def current_tocs(self, d, ix):
    for i,v in enumerate(self.tocs_rule['index']):
        if ix > v:
            break
    return self.tocs_index[i][d]


  def step(self, t, d, m, wyt, dowy, dmin=0.0, sodd=0.0): #pretty much the same as master, although there are opportunities to speed up this function by using pandas functions elsewhere
    d = self.dayofyear[t]
    dowy = water_day(d)
    m = self.month[t]
    wyt = self.wyt[t]

    envmin = max(self.env_min_flow[wyt][m-1], self.temp_releases[wyt][m-1]) * cfs_tafd
    nodd = np.interp(d, first_of_month, self.nodd)  
    
    sodd *= self.sodd_pct
    ###the variable percentage calculates Folsom & Shasta's contribution to the total releases
    ###based on the current 'available storage' in each reservoir, not a fixed percentage based on the total storage of each reservoir

    self.tocs[t] = self.current_tocs(dowy, self.fci[t])
      
    dout = dmin * self.delta_outflow_pct

    if not self.nodd_meets_envmin:
      envmin += nodd 

    # decide next release
    W = self.S[t-1] + self.Q[t]
    # HB's idea for flood control relesae..
    # fcr = (W-self.tocs[t])*np.exp(4.5*10**-6 * (W-self.capacity))
    fcr = 0.2*(W-self.tocs[t]) 
    self.Rtarget[t] = np.max((fcr, nodd+sodd+dout, envmin))

    # then clip based on constraints
    self.R[t] = min(self.Rtarget[t], W - self.dead_pool)
    self.R[t] = min(self.R[t], self.max_outflow * cfs_tafd)
    self.R[t] +=  max(W - self.R[t] - self.capacity, 0) # spill
    self.S[t] = W - self.R[t] - self.E[t] # mass balance update
    self.R_to_delta[t] = max(self.R[t] - nodd, 0) # delta calcs need this

  # def calc_EOS_storage(self,t):
  #   ##this function is called once per year in the find_available_storage function to calculate the target end-of-september storage
  #   ##at each reservoir which is used to determine how much excess storage is available for delta pumping releases
  #   self.EOS_target = (self.S[t] - self.carryover_target[self.wyt[t]])*self.carryover_excess_use + self.carryover_target[self.wyt[t]] #
  #   #we might wan't to put this function somewhere else. In my opinion is just adding to some messiness with the object-oriented structure. 

  def calc_expected_min_release(self,t):
    ##this function calculates the total expected releases needed to meet environmental minimums used in the find_available_storage function
    ##this is only calculated once per year, at the beginning of the year
    self.cum_min_release[0] = 0.0
    wyt = self.wyt[t]
    ##the cum_min_release is the total expected environmental releases between the current day and the end of september in that water year
    ## (based on the water year type)
    if self.nodd_meets_envmin: #for Shasta and Oroville only 
      for x in range(1,366):
        m = int(self.month[x-1])
        d = int(self.dayofyear[x-1])
        self.cum_min_release[0] += max(self.env_min_flow[wyt][m-1] * cfs_tafd, np.interp(d, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd) #minimum yearly release on first day. Either environmental minimum flow, north of delta demands, or temperature release standards. 
      for x in range(1,365):
        m = int(self.month[x-1])
        self.cum_min_release[x] = self.cum_min_release[x-1] - max(self.env_min_flow[wyt][m-1] * cfs_tafd , np.interp(x-1, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd ) #each day the yearly cumulative minimum release is decreased by that days minimum allowed flow. might be able to re-write this (and def take the interpolate out of the function to save time)
    else: # same idea, but for folsom. env_min_flow and nodd are combined because flow for agricultural users is diverted before the flow reaches the Lower American River (where the env minimunm flows are to be met)
      for x in range(1,366):
        m = int(self.month[x-1])
        d = int(self.dayofyear[x-1])
        self.cum_min_release[0] += max(self.env_min_flow[wyt][m-1] * cfs_tafd + np.interp(d, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd)
      for x in range(1,365):
        m = int(self.month[x-1])
        self.cum_min_release[x] = max(self.cum_min_release[x-1] - self.env_min_flow[wyt][m-1] * cfs_tafd - np.interp(x-1, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd) 


  def find_available_storage(self, t, exceedence_level):
    ##this function uses the linear regression variables calculated in find_release_func (called before simulation loop) to figure out how
    ##much 'excess' storage is available to be released to the delta with the explicit intention of running the pumps.  This function is calculated
    ##each timestep before the reservoirs' individual step function is called
    d = int(self.dayofyear[t-1])
    dowy = water_day(d)
    current_snow = self.SNPK[t-1]
    wyt = self.wyt[t]
    if dowy == 0:
      self.EOS_target = (self.S[t] - self.carryover_target[self.wyt[t]])*self.carryover_excess_use + self.carryover_target[self.wyt[t]] #
      self.calc_expected_min_release(t-1)##what do they expect to need to release for env. requirements through the end of september
      self.obs_flow = 0.0
      self.exceedence_level = (self.WYI[t-1] - 10.0)/3##how conservative are they being about the flow forecasts (ie, 90% exceedence level, 75% exceedence level, etc)
    elif dowy >= 0:
      self.obs_flow += self.Q[t-1]
      self.forecast = max(self.slope[t] * obs_flow[t] + self.intercept[t], 0.0)

    #   apr_jul_forecast = self.regression_ceoffs[dowy][2]*current_snow + self.regression_ceoffs[dowy][3]##prediction based on snowpack

    # if dowy < 182:
    #   oct_mar_forecast = self.regression_ceoffs[dowy][0]*self.oct_mar_obs + self.regression_ceoffs[dowy][1]##prediction based on total flow
    #   self.oct_mar_forecast_adjust[t] = oct_mar_forecast + self.flow_stds[dowy]*self.exceedence_level##correct for how conservative forecasts should be
    #   self.apr_jul_forecast_adjust[t] = apr_jul_forecast + self.snow_stds[dowy]*self.exceedence_level##correct for how conservative forecasts should be
    #   self.oct_mar_forecast_adjust[t] -=  self.oct_mar_obs##remove flows already observed from the forecast (linear regression is for entire period)
    #   if self.oct_mar_forecast_adjust[t] < 0.0:
    #     self.oct_mar_forecast_adjust[t] = 0.0##forcasts cannot be negative (so if observed is greater than forecasts, the extra flow will just show up in the current storage levels)
        
    # else:
    #   oct_mar_forecast = 0.0##no oct-mar forecasts are made after march (already observed) 
    #   self.oct_mar_forecast_adjust[t] = 0.0
    #   self.apr_jul_forecast_adjust[t] = apr_jul_forecast + self.snow_stds[dowy]*z_table_transform[exceedence_level]##apr-jul forecasts keep being corrected after march (but little change b/c most of the information is already baked in by April 1)
    #   self.apr_jul_forecast_adjust[t] -= self.apr_jul_obs##remove flows already observed from the forecast (lineaer regression is for entire period)
    #   if self.apr_jul_forecast_adjust[t] < 0.0:
    #     self.apr_jul_forecast_adjust[t] = 0.0##forecasts cannot be negative
      
    #available storage is storage in reservoir in exceedence of end-of-september target plus forecast for oct-mar (adjusted for already observed flow)
    #plus forecast for apr-jul (adjusted for already observed flow) minus the flow expected to be released for environmental requirements (at the reservoir, not delta)
    # self.available_storage[t] = self.S[t-1] - self.EOS_target + self.apr_jul_forecast_adjust[t] + self.oct_mar_forecast_adjust[t] - self.cum_min_release[dowy]
    self.available_storage[t] = self.S[t-1] - self.EOS_target  + self.forecast[t] - self.cum_min_release[dowy]

  # def step(self, t, d, m, wyt, dowy, dmin=0.0, sodd=0.0):

  #   envmin = self.env_min_flow[wyt][m-1] * cfs_tafd
  #   #nodd = np.interp(d, first_of_month, self.nodd)
  #   nodd = self.nodds[d]
  #   sodd *= self.sodd_pct * self.sodd_curtail_pct[wyt]
  #   self.tocs[t] = self.current_tocs(dowy, self.fci[t])
  #   dout = dmin * self.delta_outflow_pct

  #   if not self.nodd_meets_envmin:
  #     envmin += nodd 

  #   # decide next release
  #   W = self.S[t-1] + self.Q[t]
  #   # HB's idea for flood control relesae..
  #   # fcr = (W-self.tocs[t])*np.exp(4.5*10**-6 * (W-self.capacity))
  #   fcr = 0.2*(W-self.tocs[t])
  #   self.Rtarget[t] = np.max((fcr, nodd+sodd+dout, envmin)) 

  #   # then clip based on constraints
  #   self.R[t] = min(self.Rtarget[t], W - self.dead_pool)
  #   self.R[t] = min(self.R[t], self.max_outflow * cfs_tafd)
  #   self.R[t] +=  max(W - self.R[t] - self.capacity, 0) # spill
  #   self.S[t] = W - self.R[t] - self.E[t] # mass balance update
  #   self.R_to_delta[t] = max(self.R[t] - nodd, 0) # delta calcs need this

  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['storage', 'out', 'target', 'out_to_delta', 'tocs']
    things = [self.S, self.R, self.Rtarget, self.R_to_delta, self.tocs]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df


