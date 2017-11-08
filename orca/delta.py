from __future__ import division
import numpy as np 
import pandas as pd
import json
from .util import *

class Delta():

  def __init__(self, df, key):
    T = len(df)
    self.dayofyear = df.index.dayofyear
    self.month = df.index.month
    self.key = key
    self.wyt = df.SR_WYT
    self.netgains = df.netgains * cfs_tafd

    for k,v in json.load(open('orca/data/Delta_properties.json')).items():
      setattr(self,k,v)

    # what vars to store/save here
    self.dmin = np.zeros(T)
    self.gains = np.zeros(T)
    self.sodd_cvp = np.zeros(T)
    self.sodd_swp = np.zeros(T)
    self.TRP_pump = np.zeros(T)
    self.HRO_pump = np.zeros(T)
    self.inflow = np.zeros(T)
    self.outflow = np.zeros(T)
    
    self.cvp_target = np.zeros(367)
    self.swp_target = np.zeros(367)
    self.cvp_pmax = np.zeros(367)
    self.swp_pmax = np.zeros(367)

    for i in range(0,365):  
      self.cvp_target[i] = np.interp(i, self.pump_max['cvp']['d'], 
                                self.pump_max['cvp']['target']) * cfs_tafd
      self.swp_target[i] = np.interp(i, self.pump_max['swp']['d'], 
                                self.pump_max['swp']['target']) * cfs_tafd
      self.cvp_pmax[i] = np.interp(i, self.pump_max['cvp']['d'], 
                                self.pump_max['cvp']['pmax']) * cfs_tafd
      self.swp_pmax[i] = np.interp(i, self.pump_max['swp']['d'], 
                                self.pump_max['swp']['pmax']) * cfs_tafd
    
  def calc_flow_bounds(self, t, d, m, wyt, sumnodds): 
    # gains are calculated from (DeltaIn - sum of res. outflow)
    gains = self.netgains[t] + sumnodds 
    self.gains[t] = gains 

    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]

    cvp_max = self.cvp_target[d]
    swp_max = self.swp_target[d]

    # the sodd_* variables tell the reservoirs how much to release
    # for south of delta demands only
    # (dmin is the reservoir release needed to meet delta outflows)
    if gains > min_rule: # extra unstored water available for pumping
      # in this case dmin[t] is 0
      self.sodd_cvp[t] = max((cvp_max - 0.55*(gains - min_rule)) / export_ratio, 0)
      self.sodd_swp[t] = max((swp_max - 0.45*(gains - min_rule)) / export_ratio, 0)
    else: # additional flow needed
      self.dmin[t] = min_rule - gains
      # amount of additional flow from reservoirs that does not need "export tax"
      # because dmin release helps to meet the export ratio requirement
      Q = min_rule*export_ratio/(1-export_ratio) 

      if cvp_max + swp_max < Q:
        self.sodd_cvp[t] = cvp_max
        self.sodd_swp[t] = swp_max
      else:
        self.sodd_cvp[t] = 0.75*Q + (cvp_max - 0.75*Q)/export_ratio
        self.sodd_swp[t] = 0.25*Q + (swp_max - 0.25*Q)/export_ratio

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
  
  def find_pumping(self, d, dowy, t, wyt): #we can do almost all of this interpolating outside of a loop- I'm predicting this entire function will be gone as we restructure
    swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd  #limit to pumping 
    cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd #limit to pumping
    san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.sanjoaquin[t] - 1000.0 * cfs_tafd, 0.0) #not following these variables
    san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.sanjoaquin[t]
    san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    san_joaquin_ie = san_joaquin_ie_amt * san_joaquin_ie_used
    swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)

    return cvp_max, swp_max
 

  def check_san_luis(self, Tracy, Banks, t, update): # updates pumping values based on San Luis storage limits. 
  # called in calc_weekly_storage release and step functions. Tracy Pumping Plant is for CVP, Banks Pumping Plant is for SWP
    SWP_stor = self.SanLuis_SWP_stor[t-1] + Banks - self.SanLuis_SWP_out[t] #update swp storage in San Luis
    CVP_stor = self.SanLuis_CVP_stor[t-1] + Tracy - self.SanLuis_CVP_out[t] # update cvp storage in San Luis
    if SWP_stor > self.sl_cap/2.0: 
      if CVP_stor < self.sl_cap/2.0: 
        if SWP_stor + CVP_stor > self.sl_cap: #updated storages put reservoir over capacity- modifications to pumping from Banks are needed
          Banks -= SWP_stor + CVP_stor - self.sl_cap #subtract swp storage that is over capacity from Banks pumping 
          SWP_stor = self.sl_cap - CVP_stor #update swp storage because of new Banks pumping value. SL reservoir is now at maximum capacity 

      else: #updated storage values for both projects each take up more than half of the reservoirs capacity This loop alters the pumping 
      #so that each project's storage takes up exactly half of the reservoirs capacity. Loop ends with SL reservoir at maximum capacity
        Banks -= SWP_stor - self.sl_cap/2.0 #subtract swp storage that is over 1/2 capacity from Banks pumping
        Tracy -= SWP_stor - self.sl_cap/2.0 #subtract cvp storage that is over 1/2 capacity from tracy pumping
    if update == 1: #only in step function, where the pumping values are logged for the data analysis
      self.SanLuis_SWP_stor[t] = SWP_stor
      self.SanLuis_CVP_stor[t] = CVP_stor
    
    Tracy = max(Tracy,0.0)# in case updated pumping values are negative
    Banks = max(Banks,0.0) 
    
    return Tracy, Banks 
  
  
  def meet_OMR_requirement(self, Tracy, Banks, t): #old and middle river requirements (hence "OMR")
    #cvp_m = cvpm
    #swp_m = swpm
    if Tracy + Banks > self.maxTotPump: #maxTotPump is calculated in calc_weekly_storage, before this OMR function is called. 
    #current simulated puming is more that the total allowed pumping based on Delta requirements
    #Tracy (CVP) is allocated 55% of available flow for pumping, Banks (SWP) is allocated 45%. (assuming Delta outflow is greater than it's requirement- I still need to look into where that's determined)
      if Tracy < self.maxTotPump*0.55: #Tracy is pumping less that it's maximum allocated flow. Harvery should pump less flow now. 
        Banks = self.maxTotPump - Tracy 
      elif Banks < self.maxTotPump*0.45: #Banks is pumping less that it's maximum allocated flow. Tracy should pump less flow now. 
        Tracy = self.maxTotPump - Banks
      else: # in this case, both pumps would be taking their allocated percentage of flow, but the overall flow through the pumps is still greater than the maximum allowed
        Banks = self.maxTotPump*0.45
        Tracy= self.maxTotPump*0.55
    return Tracy, Banks
  def step(self, t, cvp_flows, swp_flows): #going to hold off on this alone until we discuss delta rules
    d = int(self.dayofyear[t])
    m = int(self.month[t])
    wyt = self.wyt[t]
    dowy = water_day(d)

    self.inflow[t] = self.gains[t] + cvp_flows + swp_flows
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]
    cvp_max, swp_max = self.find_pumping(d, dowy, t, wyt)
    self.TRP_pump[t] = cvp_max
    self.HRO_pump[t] = swp_max

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

  # def step(self, t, d, m, wyt, cvp_flows, swp_flows, realinflow):

  #   self.inflow[t] = max(self.gains[t] + cvp_flows + swp_flows, 0) # realinflow * cfs_tafd

  #   min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
  #   export_ratio = self.export_ratio[wyt][m-1]

  #   cvp_max = self.cvp_pmax[d]
  #   swp_max = self.swp_pmax[d]
    
  #   required_outflow = max(min_rule, (1-export_ratio)*self.inflow[t])
  #   surplus = self.gains[t] - required_outflow

  #   if surplus > 0:
  #     # gains cover both the min_rule and the export ratio requirement
  #     # so, pump the full cvp/swp inflows
  #     self.TRP_pump[t] = max(min(cvp_flows + 0.55 * surplus, cvp_max),0)
  #     self.HRO_pump[t] = max(min(swp_flows + 0.45 * surplus, swp_max),0)
  #   else:
  #     # deficit must be made up from cvp/swp flows. Assume 75/25 responsibility for these
  #     # (including meeting the export ratio requirement)
  #     deficit = -surplus
  #     cvp_pump = max(cvp_flows - 0.75 * deficit, 0)
  #     if cvp_pump == 0:
  #       swp_pump = max(swp_flows - (deficit - cvp_flows), 0)
  #     else:
  #       swp_pump = max(swp_flows - 0.25 * deficit, 0)

  #     self.TRP_pump[t] = max(min(cvp_pump, cvp_max),0)
  #     self.HRO_pump[t] = max(min(swp_pump, swp_max),0)

  #   self.outflow[t] = self.inflow[t] - self.TRP_pump[t] - self.HRO_pump[t]

    # if self.outflow[t] < min_rule:
    #   print('\nmin_rule violation: %d' % t)
    #   print('OUTFLOW: %0.2f' % (self.outflow[t] / cfs_tafd))
    #   print('INFLOW: %0.2f' % (self.inflow[t] / cfs_tafd))
    #   print('MIN_RULE: %0.2f' % (min_rule / cfs_tafd))
    #   print('EXPORT_RATIO: %0.2f' % export_ratio)
    #   print('TRACY: %0.2f' % (self.TRP_pump[t] / cfs_tafd))
    #   print('BANKS: %0.2f' % (self.HRO_pump[t] / cfs_tafd))
    #   print('GAINS: %0.2f' % (self.gains[t] / cfs_tafd))
    
    # if self.outflow[t] < (1-export_ratio)*(cvp_flows+swp_flows):
    #   print('\nexport_ratio violation: %d' % t)
      # print('OUTFLOW: %0.2f' % (self.outflow[t] / cfs_tafd))
      # print('INFLOW: %0.2f' % (self.inflow[t] / cfs_tafd))
      # print('MIN_RULE: %0.2f' % (min_rule / cfs_tafd))
      # print('EXPORT_RATIO: %0.2f' % export_ratio)
      # print('TRACY: %0.2f' % (self.TRP_pump[t] / cfs_tafd))
      # print('BANKS: %0.2f' % (self.HRO_pump[t] / cfs_tafd))
      # print('GAINS: %0.2f' % (self.gains[t] / cfs_tafd))
      
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

