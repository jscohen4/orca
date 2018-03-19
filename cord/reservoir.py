from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from .util import *


class Reservoir():

  def __init__(self, df, key):
    self.T = len(df)
    self.index = df.index
    self.key = key
    self.forecastWYT = "AN"
    for k,v in json.load(open('cord/data/%s_properties.json' % key)).items():
        setattr(self,k,v)
    
	###Reservoir Input Data
    self.Q = df['%s_inf'% key].values * cfs_tafd
    self.E = df['%s_evap'% key].values * cfs_tafd
    self.fci = df['%s_fci' % key].values
    self.SNPK = df['%s_snow' % key].values
    self.historical_storage = df['%s_storage'% key].values
    self.precip = df['%s_precip'% key].values * cfs_tafd
    self.downstream = df['%s_gains'% key].values * cfs_tafd
    self.fnf = df['%s_fnf'% key].values / 1000000.0
    self.hist_releases = df['%s_otf' % key].values * cfs_tafd
	
	##Reservoir Parameters
    self.S = np.zeros(self.T)
    self.R = np.zeros(self.T)
    self.tocs = np.zeros(self.T)
    self.Rtarget = np.zeros(self.T)
    self.R_to_delta = np.zeros(self.T)
    self.S[0] = df['%s_storage' % key].iloc[0] / 1000.0
    self.R[0] = 0
    self.EOS_target = df['%s_storage' % key].iloc[0] / 1000.0
    self.lastYearEOS_target = df['%s_storage' % key].iloc[0] / 1000.0
	
	##Forecasting available storage
    self.wytlist = ['W', 'AN', 'BN', 'D', 'C']
    self.cum_min_release = {}
    self.cum_min_release['W'] = np.zeros(366)
    self.cum_min_release['AN'] = np.zeros(366)
    self.cum_min_release['BN'] = np.zeros(366)
    self.cum_min_release['D'] = np.zeros(366)
    self.cum_min_release['C'] = np.zeros(366)
    self.oct_nov_min_release = {}
    self.oct_nov_min_release['W'] = np.zeros(366)
    self.oct_nov_min_release['AN'] = np.zeros(366)
    self.oct_nov_min_release['BN'] = np.zeros(366)
    self.oct_nov_min_release['D'] = np.zeros(366)
    self.oct_nov_min_release['C'] = np.zeros(366)
    self.aug_sept_min_release = {}
    self.aug_sept_min_release['W'] = np.zeros(366)
    self.aug_sept_min_release['AN'] = np.zeros(366)
    self.aug_sept_min_release['BN'] = np.zeros(366)
    self.aug_sept_min_release['D'] = np.zeros(366)
    self.aug_sept_min_release['C'] = np.zeros(366)
    self.available_storage = np.zeros(self.T)
    self.exceedence_level = 9
    
	##Reservoir "decisions"
    self.din = 0.0
    self.dout = 0.0
    self.envmin = 0.0	
    self.sodd = 0.0
    self.basinuse = 0.0
    self.consumed_releases = 0.0
    
    self.sjrr_release = 0.0
    self.eos_day = 0
	##Vectors for flow projections
    self.rainfnf_stds = np.zeros(365)
    self.snowfnf_stds = np.zeros(365)
    self.raininf_stds = np.zeros(365)
    self.snowinf_stds = np.zeros(365)
    self.rainflood_fnf = np.zeros(self.T)
    self.snowflood_fnf = np.zeros(self.T)
    self.rainflood_inf = np.zeros(self.T)##linear projections (i.e. 50% exceedence)
    self.snowflood_inf = np.zeros(self.T)##linear projections (i.e. 50% exceedence)
    self.rainflood_forecast = np.zeros(self.T)##projections w/ confindence (i.e. 90% exceedence - changes throughout year)
    self.snowflood_forecast = np.zeros(self.T)##projections w/ confindence (i.e. 90% exceedence - changes throughout year)
    self.evap_forecast = 0.0
	
  def find_available_storage(self, t, wyt_exceedence):
    ##this function uses the linear regression variables calculated in find_release_func (called before simulation loop) to figure out how
    ##much 'excess' storage is available to be released to the delta with the explicit intention of running the pumps.  This function is calculated
    ##each timestep before the reservoirs' individual step function is called
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    da = int(self.index.day[t])
    current_snow = self.SNPK[t]
    wyt = self.forecastWYT
	
    daysThroughMonth = [60, 91, 122, 150, 181]
	
	###Find the target end of year storage, and the expected minimum releases, at the beginning of each water year
    if m == 10 and da == 1:
      self.rainflood_flows = 0.0##total observed flows in Oct-Mar, through the current day
      self.snowflood_flows = 0.0##total observed flows in Apr-Jul, through the current day
	  
	  ##Exceedence level for flow forecasts (i.e. 2 is ~90%, 9 is ~50%)
	  ##If year starts in drought conditions, be conservative in Oct-Dec, otherwise, normal operations until January
      self.exceedence_level = 2

      ###Evap. projections are a perfect forecast (not much variation, didn't feel like making a seperate forecast for evap)
      self.evap_forecast = sum(self.E[(t):(t + 364)])
      self.eos_day = t
    if m == 8 and da == 1:
      self.lastYearEOS_target = self.EOS_target

	##Update the target EOS storage as the water year type forecasts change
    self.calc_EOS_storage(t,self.eos_day)###end-of-september target storage, based on the storage at the beginning of october
    ##Update the projected evaporation (its a perfect forecast)
    self.evap_forecast -= self.E[t]
    
	##Forecast exccedence levels set the percentile from which we forecast flows (i.e. 90% exceedence means 90% of years, with same snow conditions, would have higher flow)
	## excedence level 9 is 50% exceedence, level 1 is 90% exceedence.  Be more conservative with forecasts in drier years
    if m < 10:
      self.exceedence_level = min(m+2,9)
	  
    ##YTD observed flows (divided between rainflood and snowflood seasons) 
    if dowy < daysThroughMonth[self.melt_start]:
      self.rainflood_flows += self.Q[t]##add to the total flow observations 
    else:
      self.snowflood_flows += self.Q[t]##add to the total flow observations (
    
	###Rain- and snow-flood forecasts are predictions of future flows to come into the reservoir, for a given confidence interval (i.e. projections minus YTD observations)
    if dowy < daysThroughMonth[self.melt_start]:
      ##Forecasts are adjusted for a given exceedence level (i.e. 90%, 50%, etc)
      self.rainflood_forecast[t] = (self.rainflood_inf[t] + self.raininf_stds[dowy]*z_table_transform[self.exceedence_level]) - self.rainflood_flows
      self.snowflood_forecast[t] = (self.snowflood_inf[t] + self.snowinf_stds[dowy]*z_table_transform[self.exceedence_level])
      if self.rainflood_forecast[t] < 0:
        self.rainflood_forecast[t] = 0.0
      if self.snowflood_forecast[t] < 0:
        self.snowflood_forecast[t] = 0.0	  
    else:
      self.rainflood_forecast[t] = 0.0##no oct-mar forecasts are made after march (already observed) 
      self.snowflood_forecast[t] = (self.snowflood_inf[t] + self.snowinf_stds[dowy]*z_table_transform[self.exceedence_level]) - self.snowflood_flows
      if self.snowflood_forecast[t] < 0:
        self.snowflood_forecast[t] = 0.0	  
  	
	
    #available storage is storage in reservoir in exceedence of end-of-september target plus forecast for oct-mar (adjusted for already observed flow)
	#plus forecast for apr-jul (adjusted for already observed flow) minus the flow expected to be released for environmental requirements (at the reservoir, not delta)
    if self.S[t] < self.EOS_target and dowy > 274:
      self.available_storage[t] = 0.0
    elif m == 8 or m == 9:
      numdaysleft = 365 - dowy + 61
      self.available_storage[t] = self.S[t] - self.lastYearEOS_target + self.Q[t-1]*numdaysleft - self.oct_nov_min_release[wyt][dowy] - self.aug_sept_min_release[wyt][dowy]
      if self.key == "SHA":
        print(t, end = " ")
        print(dowy, end = " ")
        print(wyt, end = " ")
        print("%.2f" %self.available_storage[t], end = "  ")
        print("%.2f" %self.S[t], end = "  ")
        print("%.2f" %self.lastYearEOS_target, end = "  ")

    elif m == 10 or m == 11:
      numdaysleft = 61-dowy
      self.available_storage[t] = self.S[t] - self.lastYearEOS_target + self.Q[t-1]*numdaysleft - self.oct_nov_min_release[wyt][dowy]	
      if self.key == "SHA":
        print(t, end = " ")
        print(dowy, end = " ")
        print(wyt, end = " ")
        print("%.2f" %self.available_storage[t], end = "  ")
        print("%.2f" %self.S[t], end = "  ")
        print("%.2f" %self.lastYearEOS_target, end = "  ")
        print("%.2f" %self.oct_nov_min_release[wyt][dowy])
    else:
      self.available_storage[t] = self.S[t] - self.EOS_target + self.rainflood_forecast[t] + self.snowflood_forecast[t] - self.cum_min_release[wyt][dowy] - self.evap_forecast - self.aug_sept_min_release[wyt][dowy] - self.oct_nov_min_release[wyt][dowy]
	    
  def release_environmental(self,t,basinWYT):
    ###This function calculates how much water will be coming into the delta
    ###based on environmental requirements (and flood control releases).
    ###The additions to the delta contained in self.envmin represent the releases
    ###from the reservoir, minus any calls on water rights that would come from this
    ###reservoir.  This number does not include downstream 'gains' to the delta,
    ###although when those gains can be used to meet demands which would otherwise 'call'
    ###in their water rights, those gains are considered consumed before the delta but
	###no release is required from the reservoir (the reason for this is how water is 
	###accounted at the delta for dividing SWP/CVP pumping)
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    wyt = self.forecastWYT
    	
	####ENVIRONMENTAL FLOWS
	##What releases are needed directly downstream of reservoir
    self.basinuse = np.interp(d, first_of_month, self.nodd)
    self.gains_to_delta += self.basinuse
	
    if self.nodd_meets_envmin:
      reservoir_target_release = max(self.env_min_flow[wyt][m-1]*cfs_tafd - self.basinuse,0.0)
    else:
      reservoir_target_release = self.env_min_flow[wyt][m-1]*cfs_tafd

	###What releases are needed to meet flow requirements further downstream (at a point we can calculate 'gains')
    downstream_target_release = (self.temp_releases[basinWYT][m-1]*cfs_tafd - self.downstream[t])
      
	####FLOOD CONTROL
	##Top of storage pool
    self.tocs[t] = self.current_tocs(dowy, self.fci[t])
    #What size release needs to be made
    W = self.S[t] + self.Q[t]
    fcr = 0.2*(W-self.tocs[t])

	###Based on the above requirements, what flow will make it to the delta?
    self.envmin = max(reservoir_target_release, downstream_target_release,self.sjrr_release, fcr)
    self.envmin = min(self.envmin, W - self.dead_pool)
    self.envmin -= self.consumed_releases
	      
  def step(self, t):	
	###What are the contract obligations north-of-delta (only for Northern Reservoirs)
    self.envmin += (self.basinuse + self.consumed_releases)
	##What is the constraining factor - flood, consumptive demands, or environment?
    self.Rtarget[t] = self.envmin + self.sodd + self.din + self.dout
    # then clip based on constraints
    W = self.S[t] + self.Q[t]
    self.R[t] = min(self.Rtarget[t], W - self.dead_pool)
    self.R[t] = min(self.R[t], self.max_outflow * cfs_tafd)
    self.R[t] +=  max(W - self.R[t] - self.capacity, 0) # spill
    self.S[t+1] = W - self.R[t] - self.E[t] # mass balance update
	  
    self.R_to_delta[t] = max(self.R[t] - self.basinuse - self.consumed_releases, 0) # delta calcs need this

  def current_tocs(self,d,ix):
  ##Interpolates rules from tocs_rule in *_properties.json file to get the top of the conservation
  ##pool in order to determine flood control releases in reservoir.step
    for i,v in enumerate(self.tocs_rule['index']):
      if ix > v:
        break
    storage_bounds = np.zeros(2)
    index_bounds = np.zeros(2)
    storage_bounds[1] = np.interp(d, self.tocs_rule['dowy'][i-1], self.tocs_rule['storage'][i-1])
    storage_bounds[0] = np.interp(d, self.tocs_rule['dowy'][i], self.tocs_rule['storage'][i])
    index_bounds[1] = self.tocs_rule['index'][i-1]
    index_bounds[0] = self.tocs_rule['index'][i]
    return np.interp(ix, index_bounds, storage_bounds)

  def rights_call(self,downstream_flow, reset = 0):
    if reset == 0:
      if downstream_flow < 0.0:
        self.consumed_releases = downstream_flow*-1.0
        self.gains_to_delta = 0.0
      else:
        self.consumed_releases = 0.0
        self.gains_to_delta = downstream_flow
    else:
      if downstream_flow < 0.0:
        self.consumed_releases -= downstream_flow
      else:
        self.gains_to_delta += downstream_flow
	  
  def sj_riv_res_flows(self,t,d):
    ##Interpolates rules from sj_restoration_proj in *_properties.json file (note: only for Lake Millerton)
	##to get the total daily releases needed to be made under the SJ River Restoration Program (note: rules go into effect
	##in 2009 - they are triggered by the model.update_regulations function
    for i,v in enumerate(self.sj_restoration_proj['dowy']):
      if d > v:
        break
    
    ix = self.rainflood_fnf[t] + self.snowflood_fnf[t]
    release_schedule = self.sj_restoration_proj['release'][i]
    return np.interp(ix, self.sj_restoration_proj['index'], release_schedule)*cfs_tafd
	
  def calc_EOS_storage(self,t,eos_day):
    ##calculate the target end-of-september storage at each reservoir which is used to determine how much excess storage is available for delta pumping releases
    if t == 0:
      self.EOS_target = max((self.S[eos_day] - self.carryover_target[self.forecastWYT])*self.carryover_excess_use + self.carryover_target[self.forecastWYT], self.carryover_target[self.forecastWYT])
    else:
      startingStorage = max(self.S[eos_day], self.EOS_target)
      self.EOS_target = max((startingStorage - self.carryover_target[self.forecastWYT])*self.carryover_excess_use + self.carryover_target[self.forecastWYT], self.carryover_target[self.forecastWYT])

  def calc_expected_min_release(self,delta_req,delta_gains):
    ##this function calculates the total expected releases needed to meet environmental minimums used in the find_available_storage function
    ##calclulated as a pre-processing function (w/find_release_func)
    y = 1996
    for wyt in self.wytlist:
      self.cum_min_release[wyt][0] = 0.0

    ##the cum_min_release is a 365x1 vector representing each day of the coming water year.  In each day, the value is equal to 
	##the total expected minimum releases remaining in the water year, so that the 0 index is the sum of all expected releases,
	##with the total value being reduce as the index values go up, until the value is zero at the last index spot
    downstream_release = {}
    downstream_release['W'] = np.zeros(366)
    downstream_release['AN'] = np.zeros(366)
    downstream_release['BN'] = np.zeros(366)
    downstream_release['D'] = np.zeros(366)
    downstream_release['C'] = np.zeros(366)

    if self.has_downstream_target_flow:
      num_obs = np.zeros(366)
      for t in range(1,self.T):
        m = int(self.index.month[t-1])
        d = int(self.index.dayofyear[t-1])
        dowy = water_day(d)
        for wyt in self.wytlist:
          downstream_release[wyt][dowy] += max(self.temp_releases[wyt][m-1]*cfs_tafd - self.downstream[t],0.0)
		  
        num_obs[dowy] += 1
		
      for x in range(1,365):
        for wyt in self.wytlist:
          downstream_release[wyt][x] = delta_req[wyt][x]*self.delta_outflow_pct + downstream_release[wyt][x]/num_obs[x]

    if self.nodd_meets_envmin:
	   ###First, the total environmental minimum flows are summed over the whole year at day 0
      for x in range(1,365):
        for wyt in self.wytlist:
          m = int(self.index.month[x-1])
          reservoir_target_release = self.env_min_flow[wyt][m-1]*cfs_tafd
          downstream_needs = downstream_release[wyt][x]
          if x < 304:
            self.cum_min_release[wyt][0] += max(reservoir_target_release,downstream_needs)
          else:
            self.aug_sept_min_release[wyt][0] += max(reservoir_target_release,downstream_needs)	  
          if x < 61:
            self.oct_nov_min_release[wyt][0] += max(reservoir_target_release,downstream_needs)	  
	    ##THen, we loop through all 365 index spots, removing one day of releases at a time until index 365 = 0.0
      for x in range(1,365):
        for wyt in self.wytlist:
          m = int(self.index.month[x-1])
          reservoir_target_release = self.env_min_flow[wyt][m-1]*cfs_tafd
          downstream_needs = downstream_release[wyt][x]
          if x < 304:
            self.cum_min_release[wyt][x] = self.cum_min_release[wyt][x-1] - max(reservoir_target_release,downstream_needs)
            self.aug_sept_min_release[wyt][x] = self.aug_sept_min_release[wyt][0]
          else:
            self.aug_sept_min_release[wyt][x] = self.aug_sept_min_release[wyt][x-1] - max(reservoir_target_release,downstream_needs)
          if x < 61:
            self.oct_nov_min_release[wyt][x] = self.oct_nov_min_release[wyt][x-1] - max(reservoir_target_release,downstream_needs)
          else:
            self.oct_nov_min_release[wyt][x] = self.oct_nov_min_release[wyt][0]

    else:
	##Same as above, but north of delta demands are included w/ release requirements
      for x in range(1,365):
        for wyt in self.wytlist:
          m = int(self.index.month[x-1])
          reservoir_target_release = self.env_min_flow[wyt][m-1]*cfs_tafd
          downstream_needs = downstream_release[wyt][x] + np.interp(x,first_of_month, self.nodd)
          if x < 304:
            self.cum_min_release[wyt][0] += max(reservoir_target_release,downstream_needs)
          else:
            self.aug_sept_min_release[wyt][0] += max(reservoir_target_release,downstream_needs)	  
          if x < 61:
            self.oct_nov_min_release[wyt][0] += max(reservoir_target_release,downstream_needs)	  
			
      for x in range(1,365):
        for wyt in self.wytlist:
          m = int(self.index.month[x-1])
          reservoir_target_release = self.env_min_flow[wyt][m-1]*cfs_tafd
          downstream_needs = downstream_release[wyt][x] + np.interp(x,first_of_month, self.nodd)
          if x < 304:
            self.cum_min_release[wyt][x] = self.cum_min_release[wyt][x-1] - max(reservoir_target_release,downstream_needs)
            self.aug_sept_min_release[wyt][x] = self.aug_sept_min_release[wyt][0]
          else:
            self.aug_sept_min_release[wyt][x] = self.aug_sept_min_release[wyt][x-1] - max(reservoir_target_release,downstream_needs)
          if x < 61:
            self.oct_nov_min_release[wyt][x] = self.oct_nov_min_release[wyt][x-1] - max(reservoir_target_release,downstream_needs)
          else:
            self.oct_nov_min_release[wyt][x] = self.oct_nov_min_release[wyt][0]


			
  def find_release_func(self):
    ##this function is used to make forecasts when calculating available storage for export releases from reservoir
    ##using data from 1996 to 2016 (b/c data is available for all inputs needed), calculate total flows in oct-mar period and apr-jul period
    ##based on linear regression w/snowpack (apr-jul) and w/inflow (oct-mar)
    ##this function is called before simulation loop, and the linear regression coefficient & standard deviation of linear regresion residuals
    ##is used in the find_available_storage function
    rainfnf = np.zeros(self.index.year[self.T-1]-self.index.year[0] + 1)###total full-natural flow OCT-MAR
    snowfnf = np.zeros(self.index.year[self.T-1]-self.index.year[0] + 1)###total full-natural flow APR-JUL
    fnf_regression = np.zeros((365,4))##constants for linear regression: 2 for oct-mar, 2 for apr-jul
    rainfnf_cumulative = np.zeros((365,(self.index.year[self.T-1]-self.index.year[0] + 1)))###cumulative daily full-natural flow, rainfall runoff
    snowfnf_cumulative = np.zeros((365,(self.index.year[self.T-1]-self.index.year[0] + 1)))###cumulative daily full-natural flow, snowmelt runoff

    raininf = np.zeros(self.index.year[self.T-1]-self.index.year[0])##total reservoir inflow, OCT-Start of snowmelt season
    snowinf = np.zeros(self.index.year[self.T-1]-self.index.year[0])##total reservoir inflow, Start of snowmelt season - SEPT
    inf_regression = np.zeros((365,4))##constants for linear regression: 2 for oct-mar, 2 for apr-jul
    raininf_cumulative = np.zeros((365,(self.index.year[self.T-1]-self.index.year[0] + 1)))##cumulative daily reservoir inflow, rainfall runoff
    snowinf_cumulative = np.zeros((365,(self.index.year[self.T-1]-self.index.year[0] + 1)))##cumulative daily reservoir inflow, snowmelt runoff
    
    snowPattern = np.zeros((365,(self.index.year[self.T-1]-self.index.year[0] + 1)))###daily cumulative snowpack
	
    section = 0;##1 is oct-mar, 2 is apr-jul, 3 is aug-sept
    current_year = 0;
    complete_year = 0;##complete year counts all the years that have complete oct-jul data (partial years not used for lin. regression)
    for t in range(1,self.T):
      ##Get date information
      d = int(self.index.dayofyear[t-1])
      dowy = water_day(d)
      m = int(self.index.month[t-1])
      da = int(self.index.day[t-1])
      
	  #Use date information to determine if its the rainflood season
      if m == 10:
        section_fnf = 1
        section_inf = 1
      elif m == 4:
        section_fnf = 2##SRI & SJI are both divided into Oct-Mar and April-July
      elif m == 8 and da == 1:
        section_fnf = 3
        section_inf = 3
        complete_year += 1##if data exists through end of jul, counts as a 'complete year' for linear regression purposes

      if m == self.melt_start:
        section_inf = 2###for reservoir inflows, section 2 is given a variable start month, runs through end of September
		
	  #find the cumulative full natural flows through each day of the rainflood season - each day has a unique 21 value (year) vector which is the independent variable in the regression to predict total rainflood season flows
      if m == 10 and da == 1:
        current_year +=1
        rainfnf_cumulative[dowy][current_year-1] = self.fnf[t-1]
        snowfnf_cumulative[dowy][current_year-1] = 0
      elif section_fnf == 1:
        rainfnf_cumulative[dowy][current_year-1] = rainfnf_cumulative[dowy-1][current_year-1] + self.fnf[t-1]
        snowfnf_cumulative[dowy][current_year-1] = 0
      elif section_fnf == 2:
        rainfnf_cumulative[dowy][current_year-1] = rainfnf_cumulative[dowy-1][current_year-1] ##no rainflood predictions after the rainflood season ends
        snowfnf_cumulative[dowy][current_year-1] = snowfnf_cumulative[dowy-1][current_year-1] + self.fnf[t-1]
      elif section_fnf == 3:
        rainfnf_cumulative[dowy][current_year-1] = rainfnf_cumulative[dowy-1][current_year-1]
        snowfnf_cumulative[dowy][current_year-1] = snowfnf_cumulative[dowy-1][current_year-1]

	  #find the cumulative reservoir inflows through each day of the rainflood season - each day has a unique 21 value (year) vector which is the independent variable in the regression to predict total rainflood season flows
      if m == 10 and da == 1:
        raininf_cumulative[dowy][current_year-1] = self.Q[t-1]
        snowinf_cumulative[dowy][current_year-1] = 0
      elif section_inf == 1:
        raininf_cumulative[dowy][current_year-1] = raininf_cumulative[dowy-1][current_year-1] + self.Q[t-1]
        snowinf_cumulative[dowy][current_year-1] = 0
      elif section_inf == 2:
        raininf_cumulative[dowy][current_year-1] = raininf_cumulative[dowy-1][current_year-1] ##no rainflood predictions after the rainflood season ends
        snowinf_cumulative[dowy][current_year-1] = snowinf_cumulative[dowy-1][current_year-1] + self.Q[t-1]
      elif section_inf == 3:
        raininf_cumulative[dowy][current_year-1] = raininf_cumulative[dowy-1][current_year-1]
        snowinf_cumulative[dowy][current_year-1] = snowinf_cumulative[dowy-1][current_year-1]
		
	  ##find cumulative snowpack through each day of the year - each day has a unique 21 value (year) vector giving us 365 sepearte lin regressions)	  
      snowPattern[dowy][current_year-1] = self.SNPK[t-1]
 
	  #find the total full natural flows during the rainflood and snowflood seasons
      if section_fnf == 1:
        rainfnf[current_year-1] += self.fnf[t-1] ##total oct-mar inflow (one value per year - Y vector in lin regression)
      elif section_fnf == 2:
        snowfnf[current_year-1] += self.fnf[t-1] ##total apr-jul inflow (one value per year - Y vector in lin regression)

	  #find the total reservoir inflow during the rainflood and snowmelt seasons
      if section_inf == 1:
        raininf[current_year-1] += self.Q[t-1] ##total oct-mar inflow (one value per year - Y vector in lin regression)
      elif section_inf == 2:
        snowinf[current_year-1] += self.Q[t-1] ##total apr-jul inflow (one value per year - Y vector in lin regression)
    print(self.key)	
    for x in range(1,365):
      
	  ########Full natural flow regressions
	  ########################################################################################################################
	  ###rainflood season regression - full natural flow (regress cumulative full natural flow through each day with total full natural flow, Oct-Mar)
      one_year_flow = rainfnf_cumulative[x-1]##this days set of cumulative flow values (X vector)
      if sum(one_year_flow) == 0.0:
        coef[0] = 0.0
        coef[1] = np.mean(rainfnf)
      else:
        coef = np.polyfit(one_year_flow[0:(complete_year-1)],rainfnf[0:(complete_year-1)],1)###regression of cumulative flow through a day of the rain-flood season with total flow in that rain-flood season
      fnf_regression[x-1][0] = coef[0]
      fnf_regression[x-1][1] = coef[1]
      pred_dev = np.zeros(complete_year)
      for y in range(1,complete_year):
        pred_dev[y-1] = rainfnf[y-1] - coef[0]*one_year_flow[y-1] - coef[1]##how much was the linear regression off actual observations
      
      self.rainfnf_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals 
		##use z-score to make estimate at different confidence levels, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
      
	  ###snowflood season regression - full natural flow (regress cumulative snowpack & full natural flow through each day with total full natural flow, April-Jul)
      one_year_snow = snowPattern[x-1]##this days set of cumulative snowpack values (X vector)
      if sum(one_year_snow) == 0.0:
        coef[0] = 0.0
        coef[1] = np.mean(snowfnf)
      else:
        coef = np.polyfit(one_year_snow[0:(complete_year-1)],snowfnf[0:(complete_year-1)],1)
      fnf_regression[x-1][2] = coef[0]
      fnf_regression[x-1][3] = coef[1]
      pred_dev = np.zeros(complete_year)
      for y in range(1,complete_year):
        pred_dev[y-1] = snowfnf[y-1] - coef[0]*one_year_snow[y-1] - coef[1]##how much was the linear regression off actual observations

      self.snowfnf_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals
      ##for conservative estimate, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
	  #########################################################################################################################
	  
	  
	  ################Reservoir Inflow regressions
	  #########################################################################################################################
	  ###rainflood season regression - reservoir inflow (regress cumulative reservoir inflow through each day with total full natural flow, Oct-Start of Snowmelt Season at that reservroi)
      one_year_flow = raininf_cumulative[x-1]##this days set of cumulative flow values (X vector)
      if sum(one_year_flow) == 0.0:
        coef[0] = 0.0
        coef[1] = np.mean(raininf)
      else:
        coef = np.polyfit(one_year_flow[0:(complete_year-1)],raininf[0:(complete_year-1)],1)###regression of cumulative flow through a day of the rain-flood season with total flow in that rain-flood season
      inf_regression[x-1][0] = coef[0]
      inf_regression[x-1][1] = coef[1]
      pred_dev = np.zeros(complete_year)

      for y in range(1,complete_year):
        pred_dev[y-1] = raininf[y-1] - coef[0]*one_year_flow[y-1] - coef[1]##how much was the linear regression off actual observations
      
      self.raininf_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals 
		##use z-score to make estimate at different confidence levels, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
	  ###snowflood season regression - reservoir inflow (regress cumulative snowpack & reservoir inflow through each day with total reservoir inflow, Snowmelta season at the reservoir)
      one_year_snow = snowPattern[x-1]##this days set of cumulative snowpack values (X vector)
      if sum(one_year_snow) == 0.0:
        coef[0] = 0.0
        coef[1] = np.mean(snowinf)
      else:
        coef = np.polyfit(one_year_snow[0:(complete_year-1)],snowinf[0:(complete_year-1)],1)
      inf_regression[x-1][2] = coef[0]
      inf_regression[x-1][3] = coef[1]

      pred_dev = np.zeros(complete_year)
      for y in range(1,complete_year):
        pred_dev[y-1] = snowinf[y-1] - coef[0]*one_year_snow[y-1] - coef[1]##how much was the linear regression off actual observations

      self.snowinf_stds[x-1] = np.std(pred_dev)##standard deviations of linear regression residuals
      ##for conservative estimate, ie 90% exceedence is linear regression plus standard deviation * -1.28, z table in util.py
	  ############################################################################################################################################################
    current_year = 0
    for t in range(1,self.T):
      d = int(self.index.dayofyear[t-1])
      dowy = water_day(d)
      m = int(self.index.month[t-1])
      da = int(self.index.day[t-1])
	  
      if m == 10 and da == 1:
        current_year += 1
      	  
      self.rainflood_fnf[t-1] = fnf_regression[dowy-1][0]*rainfnf_cumulative[dowy][current_year-1] + fnf_regression[dowy][1]
      self.snowflood_fnf[t-1] = fnf_regression[dowy-1][2]*self.SNPK[t-1] + fnf_regression[dowy][3]
      
      self.rainflood_inf[t-1] = inf_regression[dowy-1][0]*raininf_cumulative[dowy][current_year-1] + inf_regression[dowy][1]
      self.snowflood_inf[t-1] = inf_regression[dowy-1][2]*self.SNPK[t-1] + inf_regression[dowy][3]
	  
  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['storage', 'tocs', 'available_storage', 'out']
    things = [self.S, self.tocs, self.available_storage, self.R]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df
