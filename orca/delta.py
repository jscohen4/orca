from __future__ import division
import numpy as np 
import pandas as pd
import json
from .util import *

class Delta():

  def __init__(self, df, key):
    T = len(df)
    self.index = df.index
    self.key = key
    self.wyt = df.SR_WYT
    for k,v in json.load(open('orca/data/Delta_properties.json')).items():
      setattr(self,k,v)

    # what vars to store/save here
    self.gains = np.zeros(T)
    self.gains_d = np.zeros(T)
    self.dmin = np.zeros(T)
    self.sodd_cvp = np.zeros(T)
    self.sodd_swp = np.zeros(T)
    self.TRP_pump = np.zeros(T)
    self.HRO_pump = np.zeros(T)
    self.inflow = np.zeros(T)
    self.outflow = np.zeros(T)
    self.SL_S_stor = np.zeros(T)
    self.SL_F_stor = np.zeros(T)
    self.SL_S_stor[0] = df.SLS_storage.iloc[0]
    self.SL_F_stor[0] = df.SLF_storage.iloc[0]
    self.SL_S_out = df.SLS_out
    self.SL_F_out = df.SLF_out
    
	
    self.hist_OMR = df.OMR * cfs_tafd
    self.hist_TRP_pump = df.TRP_pump * cfs_tafd
    self.hist_HRO_pump = df.HRO_pump * cfs_tafd
    self.OMR = np.zeros(T)


  def calc_flow_bounds(self, t, sum_res_inflows, folsom_nodd, shasta_nodd, oroville_nodd,orovilleAS,folsomAS,shastaAS):
    d = int(self.index.dayofyear[t])
    m = int(self.index.month[t])
    wyt = self.wyt[t]

	#northern project deliveries don't count as in-basin uses for stored/unstored flow
    folsom_n_proj_consump = np.interp(d, first_of_month, folsom_nodd)
    shasta_n_proj_consump = np.interp(d, first_of_month, shasta_nodd)
    oroville_n_proj_consump = np.interp(d, first_of_month, oroville_nodd)
    gains = self.gains_d[t] + folsom_n_proj_consump + shasta_n_proj_consump + oroville_n_proj_consump
    self.gains[t] = gains
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]

    # pump capacities (tracy's reduced to match observations)
    # cvp_max = np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['target']) * cfs_tafd
    # swp_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['target']) * cfs_tafd
    	
    folsomSODDPCT = self.calc_weekly_storage_release(t,gains,min_rule,orovilleAS,folsomAS,shastaAS)
    return folsomSODDPCT
    
  def calc_weekly_storage_release(self,t,gains,min_rule,orovilleAS,folsomAS,shastaAS):
    ##this function takes the total storage available for release at each reservoir and distributes that throughout the year based on the Export/Inflow
    ##ratio 'tax' that exists in that week
    d = self.index.dayofyear[t]
    dowy = water_day(d)
    wyt = self.wyt[t]
    m = int(self.index.month[t])
	
    omrNat = self.hist_OMR[t] + self.hist_TRP_pump[t] + self.hist_HRO_pump[t]
    fish_trigger_adj = np.interp(t, self.omr_reqr['t'], self.omr_reqr['adjustment']) * cfs_tafd
    maxTotPumpInt = omrNat - np.interp(dowy, self.omr_reqr['d'], self.omr_reqr['flow']) * cfs_tafd - fish_trigger_adj
    self.maxTotPump = max(maxTotPumpInt,0.0)
    	
    cvp_max, swp_max = self.find_release(dowy, d, t, wyt, orovilleAS, shastaAS, folsomAS)
    cvp_max, swp_max = self.meet_OMR_requirement(cvp_max, swp_max, t)
    cvp_max, swp_max = self.check_san_luis(cvp_max, swp_max, t, 0)  
   		
    surplus = self.gains[t] - self.cccpumping[t] - self.barkerpump[t] - self.depletions[t] - self.misc[t] - min_rule

    if surplus < 0: # additional flow needed to meet delta reqmnt
      self.dmin[t] = surplus * - 1.0
      self.sodd_cvp[t] = max(cvp_max/self.export_ratio[wyt][m-1] - 0.75 * (self.gains[t] + self.dmin[t]), cvp_max, 0)
      self.sodd_swp[t] = max(swp_max/self.export_ratio[wyt][m-1] - 0.25 * (self.gains[t] + self.dmin[t]), swp_max, 0)
    else: # extra unstored water available
      self.sodd_cvp[t] = max(cvp_max/self.export_ratio[wyt][m-1] - 0.55 * self.gains[t],0)
      self.sodd_swp[t] = max(swp_max/self.export_ratio[wyt][m-1] - 0.45 * self.gains[t],0)
	
    ##calculate folsom/shasta relative release contribution based on ratio of current available storage in each reservoir
    if folsomAS > 0.0 and shastaAS > 0.0:
      folsomSODDPCT = folsomAS/(folsomAS + shastaAS)
    elif folsomAS < 0.0:
      folsomSODDPCT = 0.0
    else:
      folsomSODDPCT = 1.0
    	  
    return folsomSODDPCT

  def step(self, t, cvp_flows, swp_flows):
    d = int(self.index.dayofyear[t])
    m = int(self.index.month[t])
    wyt = self.wyt[t]
    dowy = water_day(d)

    self.inflow[t] = self.gains[t] + cvp_flows + swp_flows
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]
    cvp_max, swp_max = self.find_pumping(d, dowy, t, wyt)

    surplus = self.gains[t] - self.cccpumping[t] - self.barkerpump[t] - self.depletions[t] - self.misc[t] - min_rule
    if surplus > 0:
      self.TRP_pump[t] = max(min((cvp_flows + 0.55 * self.gains[t]) * export_ratio, cvp_flows + 0.55 * surplus, cvp_max),0)
      self.HRO_pump[t] = max(min((swp_flows + 0.45 * self.gains[t]) * export_ratio, swp_flows + 0.45 * surplus, swp_max),0)
    else:
      self.TRP_pump[t] = max(min((cvp_flows + 0.75 * self.gains[t]) * export_ratio, cvp_flows + 0.75 * surplus, cvp_max),0)
      self.HRO_pump[t] = max(min((swp_flows + 0.25 * self.gains[t]) * export_ratio, swp_flows + 0.25 * surplus, swp_max),0)
    if d < 200:
      self.outflow[t] = self.inflow[t] - self.TRP_pump[t] - self.HRO_pump[t]
    	
    self.TRP_pump[t], self.HRO_pump[t] = self.meet_OMR_requirement(self.TRP_pump[t], self.HRO_pump[t], t)
    self.TRP_pump[t], self.HRO_pump[t] = self.check_san_luis(self.TRP_pump[t], self.HRO_pump[t], t, 1)
	
    self.OMR[t] = self.hist_OMR[t] + self.hist_TRP_pump[t] + self.hist_HRO_pump[t] - self.TRP_pump[t] - self.HRO_pump[t]
  
  def find_release(self, dowy, d, t, wyt, orovilleAS, shastaAS, folsomAS):
    swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
    cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd
    san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.sanjoaquin[t] - 1000.0 * cfs_tafd, 0.0)
    if np.interp(t,self.san_joaquin_export_ratio['D1641_dates'],self.san_joaquin_export_ratio['D1641_on_off']) == 1:
      san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['D1641_flow_target'],self.san_joaquin_export_ratio['D1641_export_limit']) * cfs_tafd
    else:
      san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.sanjoaquin[t]
	  
    san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    san_joaquin_ie = san_joaquin_ie_amt * san_joaquin_ie_used
    swp_jas_stor = (self.pump_max['swp']['pmax'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    cvp_jas_stor = (self.pump_max['cvp']['pmax'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    
    if dowy <= 274:
      numdaysSave = 92
    else:
      numdaysSave = 1

    if orovilleAS > numdaysSave*swp_jas_stor:
      swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    else:
      swp_max = 0.0
    if (shastaAS + folsomAS) > numdaysSave*cvp_jas_stor:
      cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)
    else:
      cvp_max = 0.0

    return cvp_max, swp_max
	
  def find_pumping(self, d, dowy, t, wyt):
    swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
    cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd
    san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.sanjoaquin[t] - 1000.0 * cfs_tafd, 0.0)
    san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.sanjoaquin[t]
    san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    san_joaquin_ie = san_joaquin_ie_amt * san_joaquin_ie_used
    swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)

    return cvp_max, swp_max
 

  def check_san_luis(self, Tracy, Harvey, t, update): # updates pumping values based on San Luis storage limits. 
  # called in calc_weekly_storage release and step functions. Tracy Pumping Plant is for CVP, Harvey Pumping Plant is for SWP
    SWP_stor = self.SL_S_stor[t-1] + Harvey - self.SL_S_out[t] #update swp storage in San Luis
    CVP_stor = self.SL_F_stor[t-1] + Tracy - self.SL_F_out[t] # update cvp storage in San Luis
    if SWP_stor > self.sl_cap/2.0: 
      if CVP_stor < self.sl_cap/2.0: 
        if SWP_stor + CVP_stor > self.sl_cap: #updated storages put reservoir over capacity- modifications to pumping from harvey are needed
          Harvey -= SWP_stor + CVP_stor - self.sl_cap #subtract swp storage that is over capacity from harvey pumping 
          SWP_stor = self.sl_cap - CVP_stor #update swp storage because of new harvey pumping value. SL reservoir is now at maximum capacity 

      else: #updated storage values for both projects each take up more than half of the reservoirs capacity This loop alters the pumping 
      #so that each project's storage takes up exactly half of the reservoirs capacity. Loop ends with SL reservoir at maximum capacity
        Harvey -= SWP_stor - self.sl_cap/2.0 #subtract swp storage that is over 1/2 capacity from harvey pumping
        Tracy -= SWP_stor - self.sl_cap/2.0 #subtract cvp storage that is over 1/2 capacity from tracy pumping
    if update == 1: #only in step function, where the pumping values are logged for the data analysis
      self.SL_S_stor[t] = SWP_stor
      self.SL_F_stor[t] = CVP_stor
	  
    Tracy = max(Tracy,0.0)# in case updated pumping values are negative
    Harvey = max(Harvey,0.0) 
    
    return Tracy, Harvey	
  
  
  def meet_OMR_requirement(self, Tracy, Harvey, t): #old and middle river requirements (hence "OMR")
    #cvp_m = cvpm
    #swp_m = swpm
    if Tracy + Harvey > self.maxTotPump: #maxTotPump is calculated in calc_weekly_storage, before this OMR function is called. 
    #current simulated puming is more that the total allowed pumping based on Delta requirements
    #Tracy (CVP) is allocated 55% of available flow for pumping, Harvey (SWP) is allocated 45%. (assuming Delta outflow is greater than it's requirement- I still need to look into where that's determined)
      if Tracy < self.maxTotPump*0.55: #Tracy is pumping less that it's maximum allocated flow. Harvery should pump less flow now. 
        Harvey = self.maxTotPump - Tracy 
      elif Harvey < self.maxTotPump*0.45: #Harvey is pumping less that it's maximum allocated flow. Tracy should pump less flow now. 
        Tracy = self.maxTotPump - Harvey
      else: # in this case, both pumps would be taking their allocated percentage of flow, but the overall flow through the pumps is still greater than the maximum allowed
        Harvey = self.maxTotPump*0.45
        Tracy= self.maxTotPump*0.55
    return Tracy, Harvey
	
  def assign_flow(self, deltadata):
    self.hist_inflows = deltadata['TOT'].values * cfs_tafd
    self.sanjoaquin = deltadata['SJR'].values * cfs_tafd
    self.cccpumping = deltadata['CCC'].values * cfs_tafd
    self.barkerpump = deltadata['NBAQ'].values * cfs_tafd
    self.depletions = deltadata['GCD'].values * cfs_tafd
    self.misc =  deltadata['MISDV'].values * cfs_tafd
	
  def find_gains(self, timestep, folsomreleases, shastareleases, orovillereleases):
      self.gains_d[timestep-1] = self.hist_inflows[timestep-1] - folsomreleases - shastareleases - orovillereleases
      
  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['in','out','TRP_pump','HRO_pump']
    things = [self.inflow, self.outflow, self.TRP_pump, self.HRO_pump]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df