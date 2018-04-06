from __future__ import division
import numpy as np 
import pandas as pd
import json
from .util import *

class Delta():

  def __init__(self, df, key, sim_gains = False):
    T = len(df)
    self.dayofyear = df.index.dayofyear
    self.month = df.index.month
    self.key = key
    self.wyt = df.WYT_sim 
    self.sim_gains = sim_gains
    if self.sim_gains:
      self.netgains = df.gains_sim #.shift(periods = -30, freq = 'D')  * cfs_tafd
    elif not self.sim_gains:
      self.netgains = df.netgains 

    for k,v in json.load(open('orca/data/Delta_properties.json')).items():
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
    # self.hist_OMR = df.OMR * cfs_tafd
    self.hist_TRP_pump = df.TRP_pump * cfs_tafd
    self.hist_HRO_pump = df.HRO_pump * cfs_tafd
    self.SanLuis_SWP_stor = np.zeros(T)
    self.SanLuis_CVP_stor = np.zeros(T)
    self.OMR = np.zeros(T)
    # self.SanLuis_SWP_stor[0] = df.SLS_storage.iloc[0]
    # self.SanLuis_CVP_stor[0] = df.SLF_storage.iloc[0]
    # self.SanLuis_SWP_out = df.SLS_out
    # self.SanLuis_CVP_out = df.SLF_out

    for i in range(0,365):  
      self.cvp_target[i] = np.interp(i, self.pump_max['cvp']['d'], 
                                self.pump_max['cvp']['target']) * cfs_tafd
      self.swp_target[i] = np.interp(i, self.pump_max['swp']['d'], 
                                self.pump_max['swp']['target']) * cfs_tafd
      self.cvp_pmax[i] = np.interp(i, self.pump_max['cvp']['d'], 
                                self.pump_max['cvp']['pmax']) * cfs_tafd
      self.swp_pmax[i] = np.interp(i, self.pump_max['swp']['d'], 
                                self.pump_max['swp']['pmax']) * cfs_tafd
    
  def calc_flow_bounds(self, t, d, m, wyt, dowy, sumnodds, orovilleAS, shastaAS, folsomAS): 

    # gains are calculated from (DeltaIn - sum of res. outflow)
    gains = self.netgains[t] + sumnodds 
    # if (m >= 3) & (m <= 9):
    #   gains = gains * 0.2
    # if (m >= 5) & (m <= 8):
    #   gains = gains - 5
    self.gains[t] = gains 

    # if m ==3:
    #   gains = gains/4

    # if m >=6 and m<=9:
    #   if gains <= 2:
    #     gains = abs(gains)*-200
    # if gains >= 0 and gains <4000:
    #   gains = gains *1.9
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]

    self.cvp_max = self.cvp_target[d]
    self.swp_max = self.swp_target[d]
    # omrNat = self.hist_OMR[t] + self.hist_TRP_pump[t] + self.hist_HRO_pump[t]
    # fish_trigger_adj = np.interp(t, self.omr_reqr['t'], self.omr_reqr['adjustment']) * cfs_tafd
    # maxTotPumpInt = omrNat - np.interp(dowy, self.omr_reqr['d'], self.omr_reqr['flow']) * cfs_tafd - fish_trigger_adj
    # self.maxTotPump = max(maxTotPumpInt,0.0)

    # the sodd_* variables tell the reservoirs how much to release
    # for south of delta demands only
    # (dmin is the reservoir release needed to meet delta outflows)
    if gains > min_rule: # extra unstored water available for pumping
      # in this case dmin[t] is 0
      self.sodd_cvp[t] = max((self.cvp_max - 0.55*(gains - min_rule)) / export_ratio, 0)
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
        self.sodd_cvp[t] = 0.75*Q + (self.cvp_max - 0.75*Q)/export_ratio
        self.sodd_swp[t] = 0.25*Q + (self.swp_max - 0.25*Q)/export_ratio

    if folsomAS > 0.0 and shastaAS > 0.0:
      self.folsomSODDPCT = folsomAS/(folsomAS + shastaAS)
    elif folsomAS < 0.0:
      self.folsomSODDPCT = 0.0
    else:
      self.folsomSODDPCT = 1.0
    self.shastaSODDPCT = 1.0 - self.folsomSODDPCT



  def find_pumping(self, d, dowy, t, wyt): #we can do almost all of this interpolating outside of a loop- I'm predicting this entire function will be gone as we restructure
    swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd  #limit to pumping 
    cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd #limit to pumping
    swp_max = swp_intake_max
    cvp_max = cvp_intake_max
    # san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.sanjoaquin[t] - 1000.0 * cfs_tafd, 0.0) #not following these variables
    san_joaquin_ie_amt = np.interp(self.sanjoaquin[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.sanjoaquin[t]
    san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    san_joaquin_ie = san_joaquin_ie_amt * san_joaquin_ie_used
    swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)

    return cvp_max, swp_max
   
  def step(self, t, d, m, wyt, dowy, cvp_flows, swp_flows, realinflow):

    self.inflow[t] = max(self.gains[t] + cvp_flows + swp_flows, 0) # realinflow * cfs_tafd

    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]
    # self.calc_weekly_storage_release(t,gains,min_rule,orovilleAS,folsomAS,shastaAS)

    cvp_max = self.cvp_pmax[d]
    swp_max = self.swp_pmax[d]
    # cvp_max, swp_max = self.find_pumping(d, dowy, t, wyt)

    required_outflow = max(min_rule, (1-export_ratio)*self.inflow[t])
    surplus = self.gains[t] - required_outflow

    if surplus > 0:
      # gains cover both the min_rule and the export ratio requirement
      # so, pump the full cvp/swp inflows
      self.TRP_pump[t] = max(min(cvp_flows + 0.55 * surplus, cvp_max),0)
      self.HRO_pump[t] = max(min(swp_flows + 0.45 * surplus, swp_max),0)
    else:
      # deficit must be made up from cvp/swp flows. Assume 75/25 responsibility for these
      # (including meeting the export ratio requirement)
      deficit = -surplus
      cvp_pump = max(cvp_flows - 0.75 * deficit, 0)
      if cvp_pump == 0:
        swp_pump = max(swp_flows - (deficit - cvp_flows), 0)
      else:
        swp_pump = max(swp_flows - 0.25 * deficit, 0)

      self.TRP_pump[t] = max(min(cvp_pump, cvp_max),0)
      self.HRO_pump[t] = max(min(swp_pump, swp_max),0)

    self.outflow[t] = self.inflow[t] - self.TRP_pump[t] - self.HRO_pump[t]
      
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

