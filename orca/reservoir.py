from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from .util import *


class Reservoir():

  def __init__(self, df, key):
    T = len(df)
    self.index = df.index
    self.key = key
    self.wyt = df.SR_WYT_rolling # 120 day MA lag
    #self.wyt = self.SRIForecast
    for k,v in json.load(open('orca/data/%s_properties.json' % key)).items():
        setattr(self,k,v)

    self.Q = df['%s_in_fix'% key].values * cfs_tafd
    self.E = df['%s_evap'% key].values * cfs_tafd
    self.fci = df['%s_fci' % key].values
    self.SNPK = df['%s_snowpack' % key].values
    self.S = np.zeros(T)
    self.R = np.zeros(T)
    self.tocs = np.zeros(T)
    self.Rtarget = np.zeros(T)
    self.R_to_delta = np.zeros(T)
    self.available_storage = np.zeros(T)
    self.oct_mar_forecast_adj = np.zeros(T)
    self.apr_jul_forecast_adj = np.zeros(T)
    self.storage_bounds = np.zeros(2)
    self.index_bounds = np.zeros(2)
    self.cum_min_release = np.zeros(366)
    self.historical_storage = df['%s_storage'% key].values
    self.S[0] = df['%s_storage' % key].iloc[0]
    self.R[0] = 0
    self.oct_mar_obs = 0.0
    self.sodd_pct_var = self.sodd_pct
    self.calc_EOS_storage(0)
    self.hist_releases = df['%s_out' % key].values * cfs_tafd

  def current_tocs(self,d,ix):
    for i,v in enumerate(self.tocs_rule['index']):
      if ix > v:
        break
    self.storage_bounds[1] = np.interp(d, self.tocs_rule['dowy'][i-1], self.tocs_rule['storage'][i-1])
    self.storage_bounds[0] = np.interp(d, self.tocs_rule['dowy'][i], self.tocs_rule['storage'][i])
    self.index_bounds[1] = self.tocs_rule['index'][i-1]
    self.index_bounds[0] = self.tocs_rule['index'][i]
    #return np.interp(d, self.tocs_rule['dowy'][i], self.tocs_rule['storage'][i])
    return np.interp(ix, self.index_bounds, self.storage_bounds)
  
  def step(self, t, dmin=0.0, sodd=0.0):
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    wyt = self.wyt[t]

    envmin = max(self.env_min_flow[wyt][m-1], self.temp_releases[wyt][m-1]) * cfs_tafd
    nodd = np.interp(d, first_of_month, self.nodd)  
	
    sodd *= self.sodd_pct_var
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

  def calc_EOS_storage(self,t):
    ##this function is called once per year in the find_available_storage function to calculate the target end-of-september storage
    ##at each reservoir which is used to determine how much excess storage is available for delta pumping releases
    self.EOS_target = (self.S[t] - self.carryover_target[self.wyt[t]])*self.carryover_excess_use + self.carryover_target[self.wyt[t]]

  def calc_expected_min_release(self,t):
    ##this function calculates the total expected releases needed to meet environmental minimums used in the find_available_storage function
    ##this is only calculated once per year, at the beginning of the year
    self.cum_min_release[0] = 0.0
    wyt = self.wyt[t]
    ##the cum_min_release is the total expected environmental releases between the current day and the end of september in that water year
    ## (based on the water year type)
    if self.nodd_meets_envmin:
      for x in range(1,366):
        m = int(self.index.month[x-1])
        d = int(self.index.dayofyear[x-1])
        self.cum_min_release[0] += max(self.env_min_flow[wyt][m-1] * cfs_tafd, np.interp(d, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd)
      for x in range(1,365):
        m = int(self.index.month[x-1])
        self.cum_min_release[x] = self.cum_min_release[x-1] - max(self.env_min_flow[wyt][m-1] * cfs_tafd , np.interp(x-1, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd )
    else:
      for x in range(1,366):
        m = int(self.index.month[x-1])
        d = int(self.index.dayofyear[x-1])
        self.cum_min_release[0] += max(self.env_min_flow[wyt][m-1] * cfs_tafd + np.interp(d, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd)

      for x in range(1,365):
        m = int(self.index.month[x-1])
        self.cum_min_release[x] = max(self.cum_min_release[x-1] - self.env_min_flow[wyt][m-1] * cfs_tafd - np.interp(x-1, first_of_month, self.nodd), self.temp_releases[wyt][m-1] * cfs_tafd) 

  def find_available_storage(self, t, exceedence_level):
    ##this function uses the linear regression variables calculated in find_release_func (called before simulation loop) to figure out how
    ##much 'excess' storage is available to be released to the delta with the explicit intention of running the pumps.  This function is calculated
    ##each timestep before the reservoirs' individual step function is called
    d = int(self.index.dayofyear[t-1])
    dowy = water_day(d)
    current_snow = self.SNPK[t-1]
    wyt = self.wyt[t]
    if dowy == 0:
      self.calc_EOS_storage(t-1)###end-of-september target storage, based on the storage at the beginning of october
      self.calc_expected_min_release(t-1)##what do they expect to need to release for env. requirements through the end of september
      self.oct_mar_obs = 0.0##total observed flows in Oct-Mar, through the current day
      self.apr_jul_obs = 0.0##total observed flows in Apr-Jul, through the current day
      self.exceedence_level = (self.SRIforecast[t-1] - 10.0)/3##how conservative are they being about the flow forecasts (ie, 90% exceedence level, 75% exceedence level, etc)

    if dowy < 182:
      self.oct_mar_obs += self.Q[t-1]##add to the total flow observations (only through March)
    else:
      self.apr_jul_obs += self.Q[t-1]##add to the total flow observations (starting in April)
    
    apr_jul_forecast = self.regression_ceoffs[dowy][2]*current_snow + self.regression_ceoffs[dowy][3]##prediction based on snowpack

    if dowy < 182:
      oct_mar_forecast = self.regression_ceoffs[dowy][0]*self.oct_mar_obs + self.regression_ceoffs[dowy][1]##prediction based on total flow
      self.oct_mar_forecast_adj[t] = oct_mar_forecast + self.flow_stds[dowy]*self.exceedence_level##correct for how conservative forecasts should be
      self.apr_jul_forecast_adj[t] = apr_jul_forecast + self.snow_stds[dowy]*self.exceedence_level##correct for how conservative forecasts should be
      self.oct_mar_forecast_adj[t] -=  self.oct_mar_obs##remove flows already observed from the forecast (linear regression is for entire period)
      if self.oct_mar_forecast_adj[t] < 0.0:
        self.oct_mar_forecast_adj[t] = 0.0##forecasts cannot be negative (so if observed is greater than forecasts, the extra flow will just show up in the current storage levels)
		
    else:
      oct_mar_forecast = 0.0##no oct-mar forecasts are made after march (already observed) 
      self.oct_mar_forecast_adj[t] = 0.0
      self.apr_jul_forecast_adj[t] = apr_jul_forecast + self.snow_stds[dowy]*z_table_transform[exceedence_level]##apr-jul forecasts keep being corrected after march (but little change b/c most of the information is already baked in by April 1)
      self.apr_jul_forecast_adj[t] -= self.apr_jul_obs##remove flows already observed from the forecast (lineaer regression is for entire period)
      if self.apr_jul_forecast_adj[t] < 0.0:
        self.apr_jul_forecast_adj[t] = 0.0##forecasts cannot be negative
	  
    #available storage is storage in reservoir in exceedence of end-of-september target plus forecast for oct-mar (adjusted for already observed flow)
	#plus forecast for apr-jul (adjusted for already observed flow) minus the flow expected to be released for environmental requirements (at the reservoir, not delta)
    self.available_storage[t] = self.S[t-1] - self.EOS_target + self.apr_jul_forecast_adj[t] + self.oct_mar_forecast_adj[t] - self.cum_min_release[dowy]
	
  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['storage', 'out', 'target', 'out_to_delta', 'tocs', 'available_storage', 'apr_for', 'oct_for']
    things = [self.S, self.R, self.Rtarget, self.R_to_delta, self.tocs, self.available_storage, self.apr_jul_forecast_adj, self.oct_mar_forecast_adj,]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df

	
  def find_release_func(self, data): #still will need to work on this for the purposes of speeding things up- that'll be done when merging to master branch
    ##this function is used to make forecasts when calculating available storage for export releases from reservoir
    ##using data from 1996 to 2016 (b/c data is available for all inputs needed), calculate total flows in oct-mar period and apr-jul period
    ##based on linear regression w/snowpack (apr-jul) and w/inflow (oct-mar)
    ##this function is called before simulation loop, and the linear regression coefficient & standard deviation of linear regresion residuals
    ##is used in the find_available_storage function
    ## data used is from release-cdf-data.csv
    self.cum_snow = data['%s_cdf_snow'% self.key].values##cumulative yearly snowpack on each dayin data set
    self.daily_inflow = data['%s_cdf_inf'% self.key].values##cumulative oct-mar inflow (Based on how the data look, I'm inclined to think this is daily inflow- it matches that data in the master branch)
    #time_of_year = 0;##1 is oct-mar, 2 is apr-jul, 3 is aug-sept
    current_year = 0;
    complete_year = 0;##complete year counts all the years that have complete oct-jul data (partial years not used for lin. regression)
    self.oct_mar_cum_inflows = np.zeros(data.index.year[len(data)-1]-data.index.year[0]) #culmulative yearly october- march inflows (to be calculated), one per year
    self.apr_jul_cum_inflows = np.zeros(data.index.year[len(data)-1]-data.index.year[0]) #culmulative yearly april-july inflows (to be calculated), one per year
    self.regression_ceoffs = np.zeros((365,4))##coefficients for linear regression: 2 for oct-mar, 2 for apr-jul
    self.snow_stds = np.zeros(365) #standard deviations for snowpack regressions
    self.flow_stds = np.zeros(182) #standard deviations for oct-march flow regressions
    self.cum_snow_matrix = np.zeros((365,(data.index.year[len(data)-1]-data.index.year[0]))) #will become matrix with cumulative snowpack for each date, structured by year (365 x # of years matrix)
    self.daily_inflow_matrix = np.zeros((365,(data.index.year[len(data)-1]-data.index.year[0])))#will become matrix with daily inflows for each oct-mar date, structured by year (365 x # of years matrix) -only oct-mar inflows are included
	  
    for t in range(1,len(data)): #I'll see if I can make this more consise, although I think it's good for now
      d = int(data.index.dayofyear[t-1])
      dowy = water_day(d)
      m = int(data.index.month[t-1])
      y = int(data.index.year[t-1])
      day = int(data.index.day[t-1])
      if dowy == 1:
        current_year +=1
	  
      if m == 10: #october- beggining of water year
        time_of_year = "oct-mar"
      elif m == 4: # april- start of summer
        time_of_year = "apr-jul"
      elif m == 8 and day == 1: # august- end of apr-july regression
        time_of_year = "aug-sept"
        complete_year += 1##if data exists through jul, counts as a 'complete year' for linear regression purposes
      
      if time_of_year == "oct-mar":
        self.oct_mar_cum_inflows[current_year-1] += self.daily_inflow[t-1] * cfs_tafd##total oct-mar inflow (one value per year - Y vector in lin regression)
      elif time_of_year == "apr-jul":
        self.apr_jul_cum_inflows[current_year-1] += self.daily_inflow[t-1] * cfs_tafd##total apr-jul inflow (one value per year - Y vector in lin regression)
    	
    current_year = 0;
    for t in range(1,len(data)):
      d = int(data.index.dayofyear[t-1])
      dowy = water_day(d)
      m = int(data.index.month[t-1])
      y = int(data.index.year[t-1])
      #cum oct-mar inflows through each day(X vector in oct-mar lin regression - each day has a unique 21 value (year) vector giving us 365 seperate lin regressions)
      if dowy == 1:
        current_year += 1;
        self.daily_inflow_matrix[dowy-1][current_year-1] = self.daily_inflow[t-1] * cfs_tafd
      elif dowy < 182:
        self.daily_inflow_matrix[dowy-1][current_year-1] = self.daily_inflow_matrix[dowy-2][current_year-1] + self.daily_inflow[t-1] * cfs_tafd
      
      ##cum snowpack through each day (X vector in apr-jul lin regression - each day has a unique 21 value (year) vector giving us 365 sepearte lin regressions)	  
      self.cum_snow_matrix[dowy-1][current_year-1] = self.cum_snow[t-1]
      
    for x in range(1,182):
      flow_each_year = self.daily_inflow_matrix[x-1]##this days set of cumulative flow values, one value for each year(X vector) #still not seeing why this is culmalative
      coef = np.polyfit(flow_each_year[0:(complete_year-1)],self.oct_mar_cum_inflows[0:(complete_year-1)],1) #coef is the set of two regression coeffiients for this year's flow regression
      self.regression_ceoffs[x-1][0] = coef[0]
      self.regression_ceoffs[x-1][1] = coef[1]
      pred_dev = np.zeros(complete_year)
      for y in range(1,complete_year):
        pred_dev[y-1] = self.oct_mar_cum_inflows[y-1] - coef[0]*flow_each_year[y-1] - coef[1]##how much was the linear regression off actual observations

      self.flow_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals 
      ##for conservative estimate, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
    for x in range(1,365):
      snow_each_year = self.cum_snow_matrix[x-1]##this days set of cumulative snowpack values (X vector)  one value for each year
      coef = np.polyfit(snow_each_year[0:(complete_year-1)],self.apr_jul_cum_inflows[0:(complete_year-1)],1) #coef is the set of two regression coeffiients for this year's snowpack regression
      self.regression_ceoffs[x-1][2] = coef[0]
      self.regression_ceoffs[x-1][3] = coef[1]
      pred_dev = np.zeros(complete_year)
      for y in range(1,complete_year):
        pred_dev[y-1] = self.apr_jul_cum_inflows[y-1] - coef[0]*snow_each_year[y-1] - coef[1]##how much was the linear regression off actual observations

      self.snow_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals 
      ##for conservative estimate, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
      
