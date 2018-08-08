from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from .util import *

class Reservoir():

  def __init__(self, df, dfh, key, projection = False):
    T = len(df)
    self.projection = projection
    self.dayofyear = df.index.dayofyear
    self.month = df.index.month
    self.key = key
    self.wyt = df.WYT_sim# simulated (forecasted)wyi
    for k,v in json.load(open('orca/data/json_files/%s_properties.json' % key)).items():
      setattr(self,k,v)
    self.evap_reg = json.load(open('orca/data/json_files/evap_regression.json'))
    self.evap_coeffs = np.asarray(self.evap_reg['%s_evap_coeffs' % key])
    self.evap_int = self.evap_reg['%s_evap_int' % key]
    # self.sodd_pct_var = self.sodd_pct
    if self.projection:
        self.Q = df['%s_fnf'% key].values * cfs_tafd
        self.E = np.zeros(T)
    if not self.projection:
        self.Q = df['%s_in_tr'%key].values * cfs_tafd
        self.E = df['%s_evap'% key].values * cfs_tafd
    self.fci = df['%s_fci' % key].values
    self.slope =  df['%s_slope' % key].values
    self.intercept = df['%s_intercept' % key].values
    self.rem_flow = df['%s_remaining_flow' % key].values

    self.mean = df['%s_mean' % key].values
    self.std = df['%s_std' % key].values  
    self.tas = df['%s_tas' % key].values
    self.WYI = df['WYI_sim'].values
    self.obs_flow = df['%s_cum_flow_to_date' % key].values
    self.obs_snow = df['%s_snowpack' % key].values
    #initialize time series arrays
    self.S = np.zeros(T)
    self.R = np.zeros(T)
    self.Rtarget = np.zeros(T)
    self.R_to_delta = np.zeros(T)
    self.S[0] = dfh['%s_storage' % key].iloc[0]
    self.R[0] = 0
    self.storage_bounds = np.zeros(2)
    self.index_bounds = np.zeros(2)
    self.tocs = np.zeros(T)
    self.cum_min_release = np.zeros(366)
    self.forecast = np.zeros(T)
    self.available_storage = np.zeros(T)
    self.soddp = np.zeros(T)
    self.spill = np.zeros(T)

    #tocs rule variables
    self.tocs_index = []
    for i,v in enumerate(self.tocs_rule['index']):        
        self.tocs_index.append(np.zeros(366))
        for day in range(0, 366):  
            self.tocs_index[i][day] = np.interp(day, self.tocs_rule['dowy'][i], self.tocs_rule['storage'][i])
    if key == 'ORO':
    	self.tocs_index = np.roll(self.tocs_index, -16, axis = 1)

    self.nodds = np.zeros(367)
    for i in range(0,366):
        self.nodds[i] = np.interp(i, first_of_month, self.nodd)


  def current_tocs(self, d, ix):
    for i,v in enumerate(self.tocs_rule['index']):
        if ix > v:
            break
    return self.tocs_index[i][d]


  def step(self, t, d, m, wyt, dowy, dmin=0.0, sodd=0.0, projection = False): 
    d = self.dayofyear[t]
    dowy = water_day(d)
    m = self.month[t]
    wyt = self.wyt[t]
    envmin = self.env_min_flow[wyt][m-1] * cfs_tafd #minimum allowed environmental flows
    nodd = np.interp(d, first_of_month, self.nodd)  #north of delta demands
    # sodd *= self.sodd_pct_var
    sodd *= self.sodd_pct * self.sodd_curtail_pct[wyt] #south of delta demands
    self.soddp[t] = sodd

    ###the variable percentage calculates Folsom & Shasta's contribution to the total releases
    ###based on the current 'available storage' in each reservoir, not a fixed percentage based on the total storage of each reservoir
    self.tocs[t] = self.current_tocs(dowy, self.fci[t])
    dout = dmin * self.delta_outflow_pct
    if not self.nodd_meets_envmin:
      envmin += nodd 
    # decide next release
    W = self.S[t-1] + self.Q[t]
    fcr = 0.2*(W-self.tocs[t]) #flood control release
    self.Rtarget[t] = np.max((fcr, nodd+sodd+dout, envmin)) #target release
    if self.carryover_rule:
        if m >= 5 and m <= 9: #from may to september, accout for carryover storage targets
            if self.forecast[t] + self.S[t-1] - self.Rtarget[t] * (365-dowy) < self.carryover_target[wyt]: #forecasting rest-of-wateryear inflow
                #how much to curtail releases in attempt to meet carryover targets
                self.carryover_curtail_pct = (self.forecast[t] + self.S[t-1] - self.Rtarget[t] * (365-dowy))/self.carryover_target[wyt]
                self.Rtarget[t] = self.Rtarget[t] * max(self.carryover_curtail_pct,self.carryover_curtail[wyt]) #update target release

    # then clip based on constraints
    self.R[t] = min(self.Rtarget[t], W - self.dead_pool) # dead-pool constraint
    self.R[t] = min(self.R[t], self.max_outflow * cfs_tafd) #max outflow constraint
    self.R[t] +=  max(W - self.R[t] - self.capacity, 0) # spill
    self.spill[t] = max(W - self.R[t] - self.capacity +  self.Q[t],0)

    if self.projection: # if projection mode, calculate evaporation
        X=[]
        storage = self.S[t-1]
        temp = self.tas[t]
        X.append(temp)
        X.append(storage)
        X.append(temp*storage)
        X.append(temp**2)
        X.append(storage**2)
        self.E[t] = max((np.sum(X * self.evap_coeffs) + self.evap_int) * cfs_tafd,0) #evaporation variable
    self.S[t] = W - self.R[t] - self.E[t] # mass balance update
    self.R_to_delta[t] = max(self.R[t] - nodd, 0) # delta calcs need this

  def calc_expected_min_release(self,t):
    ##this function calculates the total expected releases needed to meet environmental minimums used in the find_available_storage function
    ##this is only calculated once per year, at the beginning of the year
    self.cum_min_release[0] = 0.0
    wyt = self.wyt[t]
    ##the cum_min_release is the total expected environmental releases between the current day and the end of september in that water year
    ## (based on the water year type)
    if self.nodd_meets_envmin: #for Shasta and Oroville only 
      for x in range(1,366):
        m = self.month[x-1]
        d = self.dayofyear[x-1]
        self.cum_min_release[0] += max(self.env_min_flow[wyt][m-1] * cfs_tafd, self.nodds[d], self.temp_releases[wyt][m-1] * cfs_tafd) #minimum yearly release on first day. Either environmental minimum flow, north of delta demands, or temperature release standards. 
      for x in range(1,365):
        m = self.month[x-1]
        self.cum_min_release[x] = self.cum_min_release[x-1] - max(self.env_min_flow[wyt][m-1] * cfs_tafd, self.nodds[d], self.temp_releases[wyt][m-1] * cfs_tafd ) #each day the yearly cumulative minimum release is decreased by that days minimum allowed flow. might be able to re-write this (and def take the interpolate out of the function to save time)
    else: # same idea, but for folsom. env_min_flow and nodd are combined because flow for agricultural users is diverted before the flow reaches the Lower American River (where the env minimunm flows are to be met)
      for x in range(1,366):
        m = self.month[x-1]
        d = self.dayofyear[x-1]
        self.cum_min_release[0] += max(self.env_min_flow[wyt][m-1] * cfs_tafd + self.nodds[d], self.temp_releases[wyt][m-1] * cfs_tafd)
      for x in range(1,365):
        m = self.month[x-1]
        self.cum_min_release[x] = max(self.cum_min_release[x-1] - self.env_min_flow[wyt][m-1] * cfs_tafd - self.nodds[d], self.temp_releases[wyt][m-1] * cfs_tafd) 

  def find_available_storage(self, t, d,dowy):#, exceedence_level):
    ##this function uses the linear regression variables calculated in find_release_func (called before simulation loop) to figure out how
    ##much 'excess' storage is available to bešreleased to the delta with the explicit intention of running the pumps.  This function is calculated
    ##each timestep before the reservoirs' individual step function is called
    #also used to obtain inflow forecasts
    self.exceedence_level = -1*min((self.WYI[t-1] - 10.0)*0.8,-2)##how conservative are they being about the flow forecasts (ie, 90% exceedence level, 75% exceedence level, etc)
    self.forecast[t] = max(0,self.slope[t] * self.obs_snow[t] + self.intercept[t] + self.std[t]*z_table_transform[self.exceedence[self.wyt[t]]]) * 1000 #based on forecast regression
    if dowy == 0:
      self.calc_expected_min_release(t)##what do they expect to need to release for env. requirements through the end of september
      self.forecast[t] = max(0,self.slope[t+1] * self.obs_snow[t+1] + self.intercept[t+1]+ self.std[t]*z_table_transform[self.exceedence[self.wyt[t]]]) * 1000 #based on forecast regression
    # if d == 0:
      # self.forecast[t] = max(0,self.slope[t] * self.obs_snow[t] + self.intercept[t]) * 1000 #based on forecast regression
    # self.available_storage[t] = max(0,self.S[t-1] - self.carryover_target[self.wyt[t]] + self.forecast[t] - self.cum_min_release[dowy])
    self.available_storage[t] = max(0,self.S[t-1] - self.carryover_target[self.wyt[t]]/self.exceedence_level + self.forecast[t] - self.cum_min_release[dowy])

  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['storage', 'out', 'target', 'out_to_delta', 'tocs','sodd','spill']
    things = [self.S, self.R, self.Rtarget, self.R_to_delta, self.tocs,self.soddp,self.spill]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df