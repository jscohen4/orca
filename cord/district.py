from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from .util import *


class District():

  def __init__(self, key):
    self.key = key
    for k,v in json.load(open('cord/districts/%s_properties.json' % key)).items():
        setattr(self,k,v)

    self.irrdemand = json.load(open('cord/crop/%s_properties.json' % self.zone)).items()
	
  def find_available_storage(self, t):
    ##this function uses the linear regression variables calculated in find_release_func (called before simulation loop) to figure out how
    ##much 'excess' storage is available to be released to the delta with the explicit intention of running the pumps.  This function is calculated
    ##each timestep before the reservoirs' individual step function is called
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    da = int(self.index.day[t])
    current_snow = self.SNPK[t]
    wyt = self.forecastWYT[t]
	
    daysThroughMonth = [60, 91, 122, 150, 181]
	
	###Find the target end of year storage, and the expected minimum releases, at the beginning of each water year
    if m == 10 and da == 1:
      self.calc_expected_min_release(t,dowy)##what do they expect to need to release for env. requirements through the end of september
      self.rainflood_flows = 0.0##total observed flows in Oct-Mar, through the current day
      self.snowflood_flows = 0.0##total observed flows in Apr-Jul, through the current day
      if wyt == 'D' or wyt == 'C':
        self.exceedence_level = 2
      else:
        self.exceedence_level = 9##reset the exccedence level used in forecasts
      
      self.eos_day = t
    
    self.calc_EOS_storage(t,self.eos_day)###end-of-september target storage, based on the storage at the beginning of october	
    if m < 10:
      self.exceedence_level = min(m + 2,9)
    ##Split observed flows (through a given day) into the rainflood and snowflood seasons	
    if dowy < daysThroughMonth[self.melt_start]:
      self.rainflood_flows += self.Q[t]##add to the total flow observations 
    else:
      self.snowflood_flows += self.Q[t]##add to the total flow observations (
    
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
  	
    self.evap_for[t] = sum(self.E[(t):(t - dowy + 364)])
    self.remaining_release[t] = self.cum_min_release[dowy]
    self.storage_target[t] = self.EOS_target
    #available storage is storage in reservoir in exceedence of end-of-september target plus forecast for oct-mar (adjusted for already observed flow)
	#plus forecast for apr-jul (adjusted for already observed flow) minus the flow expected to be released for environmental requirements (at the reservoir, not delta)
    if self.S[t] < self.storage_target[t] and dowy > 274:
      self.available_storage[t] = 0.0
    else:
      self.available_storage[t] = self.S[t] - self.storage_target[t] + self.rainflood_forecast[t] + self.snowflood_forecast[t] - self.remaining_release[t] - self.evap_for[t]
	  
	
