from __future__ import division
import numpy as np 
import pandas as pd
import json
from .util import *

class Delta():

  def __init__(self, df, key, sim_gains):
    T = len(df)
    self.dayofyear = df.index.dayofyear
    self.month = df.index.month
    self.key = key
    self.wyt = df.WYT_sim 
    self.sim_gains = sim_gains
    if self.sim_gains:
      self.netgains = df.gains_sim
      self.sanjoaquin = df.gains_sim - df.YRS_fnf - df.NML_fnf

    elif not self.sim_gains:
      self.netgains = df.netgains 
      self.sanjoaquin = df.netgains - df.YRS_fnf - df.NML_fnf
    for k,v in json.load(open('orca/data/json_files/Delta_properties.json')).items():
      setattr(self,k,v)
    # self.assign_flow(df)
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
    self.swp_intake_max = np.zeros(367)
    self.cvp_intake_max = np.zeros(367)
    self.san_joaquin_adj = np.zeros(367)
    self.D1641_on_off = np.zeros(367)
    self.san_joaquin_ie_used = np.zeros(367)
    self.san_joaquin_ie_amt = np.zeros(T)
    for i in range(0,365):  
      self.cvp_target[i] = np.interp(i, self.pump_max['cvp']['d'], #calculate pumping target for day of year (based on target pumping for sodd) 
                                self.pump_max['cvp']['target']) * cfs_tafd
      self.swp_target[i] = np.interp(i, self.pump_max['swp']['d'], 
                                self.pump_max['swp']['target']) * cfs_tafd
      self.cvp_pmax[i] = np.interp(i, self.pump_max['cvp']['d'], 
                                self.pump_max['cvp']['pmax']) * cfs_tafd #calculate pumping targets (based on max allowed pumping) based on time of year 
      self.swp_pmax[i] = np.interp(i, self.pump_max['swp']['d'], 
                                self.pump_max['swp']['pmax']) * cfs_tafd

      self.swp_intake_max[i] = np.interp(i, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
      self.cvp_intake_max[i] = np.interp(i, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd
      self.san_joaquin_adj[i] = np.interp(water_day(i), self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.sanjoaquin[i] - 1000.0 * cfs_tafd, 0.0)
      self.san_joaquin_ie_used[i] = np.interp(water_day(i), self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    
    for i in range(0,T):
      self.san_joaquin_ie_amt[i] = np.interp(self.sanjoaquin[i]*tafd_cfs, self.san_joaquin_export_ratio['D1641_flow_target'],self.san_joaquin_export_ratio['D1641_export_limit']) * cfs_tafd

  def find_release(self, dowy, d, t, wyt, orovilleAS, shastaAS, folsomAS):
    # swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
    # cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd
    
    # san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.sanjoaquin[t] - 1000.0 * cfs_tafd, 0.0)
    # if np.interp(t,self.san_joaquin_export_ratio['D1641_dates'],self.san_joaquin_export_ratio['D1641_on_off']) == 1:
    # san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['D1641_flow_target'],self.san_joaquin_export_ratio['D1641_export_limit']) * cfs_tafd
    # else:
      # san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.sanjoaquin[t]
    # san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])

    san_joaquin_ie = self.san_joaquin_ie_amt[t] * self.san_joaquin_ie_used[dowy]
    swp_jas_stor = (self.pump_max['swp']['pmax'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    cvp_jas_stor = (self.pump_max['cvp']['pmax'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    
    if dowy <= 274:
      numdaysSave = 92
    else:
      numdaysSave = 1
    # swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    # cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)

    if orovilleAS > numdaysSave:#*swp_jas_stor:
      swp_max = min(max(self.swp_intake_max[d] + self.san_joaquin_adj[d], san_joaquin_ie * 0.45), self.swp_pmax[d])
    else:
      swp_max = 0.0
    if (shastaAS + folsomAS) > numdaysSave:#*cvp_jas_stor:
      cvp_max = min(max(self.cvp_intake_max[d], san_joaquin_ie * 0.55), self.cvp_pmax[d])
    else:
      cvp_max = 0.0

    return cvp_max, swp_max

  def calc_flow_bounds(self, t, d, m, wyt, dowy, sumnodds, orovilleAS, shastaAS, folsomAS): 

    # gains are calculated from (DeltaIn - sum of res. outflow)

    gains = self.netgains[t] 
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]

    self.cvp_max = self.cvp_target[d-1]
    self.swp_max = self.swp_target[d-1]
    if d == 366:
      cvp_max = self.cvp_target[d-2]
      swp_max = self.swp_target[d-2]

    # the sodd_* variables tell the reservoirs how much to release
    # for south of delta demands only
    # (dmin is the reservoir release needed to meet delta outflows)
    if gains > min_rule: # extra unstored water available for pumping
      # in this case dmin[t] is 0
      self.sodd_cvp[t] = max((self.cvp_max - 0.55*(gains - min_rule)) / export_ratio, 0) #implementing export ratio "tax"
      self.sodd_swp[t] = max((self.swp_max - 0.45*(gains - min_rule)) / export_ratio, 0)
    else: # additional flow needed
      self.dmin[t] = min_rule - gains
      # amount of additional flow from reservoirs that does not need "export tax"
      # because dmin release helps to meet the export ratio requirement
      Q = min_rule*export_ratio/(1-export_ratio) 

      if self.cvp_max + self.swp_max < Q:
        self.sodd_cvp[t] = self.cvp_max
        self.sodd_swp[t] = self.swp_max
      else:
        self.sodd_cvp[t] = 0.75*Q + (self.cvp_max - 0.75*Q)/export_ratio #implementing export ratio "tax"
        self.sodd_swp[t] = 0.25*Q + (self.swp_max - 0.25*Q)/export_ratio

    #determining percentage of CVP sodd demands from both Shasta and Folsom
    if folsomAS > 0.0 and shastaAS > 0.0:
      self.folsomSODDPCT = folsomAS/(folsomAS + shastaAS)
    elif folsomAS < 0.0:
      self.folsomSODDPCT = 0.0
    else:
      self.folsomSODDPCT = 1.0
    self.shastaSODDPCT = 1.0 - self.folsomSODDPCT
   
  def step(self, t, d, m, wyt, dowy, cvp_flows, swp_flows, orovilleAS, shastaAS, folsomAS,sumnodds):
    self.gains[t] = self.netgains[t] + sumnodds
    self.inflow[t] = max(self.gains[t] + cvp_flows + swp_flows, 0) # realinflow * cfs_tafd
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]
    cvp_max = self.cvp_pmax[d-1] #max pumping allowed 
    swp_max = self.swp_pmax[d-1]
    cvp_max, swp_max = self.find_release(dowy, d, t, wyt, orovilleAS, shastaAS, folsomAS)

    if d == 366:
      cvp_max = self.cvp_pmax[d-2]
      swp_max = self.swp_pmax[d-2]
    required_outflow = max(min_rule, (1-export_ratio)*self.inflow[t])
    surplus = self.gains[t] - required_outflow 

    if surplus > 0:
      # gains cover both the min_rule and the export ratio requirement
      # so, pump the full cvp/swp inflows
      self.TRP_pump[t] = max(min(cvp_flows + 0.55 * surplus, cvp_max),0) #Tracy pumping plant, for CVP exports
      self.HRO_pump[t] = max(min(swp_flows + 0.45 * surplus, swp_max),0) #Harvey 0. Banks pumping plant, for SWP exports
    else:
      # deficit must be made up from cvp/swp flows. Assume 75/25 responsibility for these
      # (including meeting the export ratio requirement)
      deficit = -surplus
      cvp_pump = max(cvp_flows - 0.75 * deficit, 0)
      if cvp_pump == 0:
        swp_pump = max(swp_flows - (deficit - cvp_flows), 0)
      else:
        swp_pump = max(swp_flows - 0.25 * deficit, 0)
      self.TRP_pump[t] = max(min(cvp_pump, cvp_max),0) #overall TRP pumping
      self.HRO_pump[t] = max(min(swp_pump, swp_max),0) #overall HRO pumping

    self.outflow[t] = self.inflow[t] - self.TRP_pump[t] - self.HRO_pump[t]

  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['in','out','TRP_pump','HRO_pump']
    things = [self.inflow, self.outflow, self.TRP_pump, self.HRO_pump]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df