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
    self.wyt = df.SR_WYT_rolling
    self.gains_sac = df.SAC_gains
    self.gains_d = df.Delta_gains

    for k,v in json.load(open('cord/data/Delta_properties.json')).items():
      setattr(self,k,v)

    # what vars to store/save here
    self.gains = np.zeros(T)
    self.dmin = np.zeros(T)
    self.sodd_cvp = np.zeros(T)
    self.sodd_swp = np.zeros(T)
    self.TRP_pump = np.zeros(T)
    self.HRO_pump = np.zeros(T)
    self.inflow = np.zeros(T)
    self.outflow = np.zeros(T)

  def calc_flow_bounds(self, t, sum_res_inflows, folsom_nodd, shasta_nodd, oroville_nodd):
    d = self.index.dayofyear[t]
    m = self.index.month[t]
    wyt = self.wyt[t]

    # gains = sum_res_inflows * self.gains_factor[m-1]
    # self.gains[t] = gains
	
	#northern project deliveries don't count as in-basin uses for stored/unstored flow
    folsom_n_proj_consump = np.interp(d, first_of_month, folsom_nodd)
    shasta_n_proj_consump = np.interp(d, first_of_month, shasta_nodd)
    oroville_n_proj_consump = np.interp(d, first_of_month, oroville_nodd)
    gains = self.gains_d[t] * cfs_tafd + folsom_n_proj_consump + shasta_n_proj_consump + oroville_n_proj_consump
    self.gains[t] = gains
    #min_rule = np.interp(d, first_of_month, self.min_outflow[wyt]) * cfs_tafd
    #export_ratio = np.interp(d, first_of_month, self.export_ratio[wyt])
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]

    # pump capacities (tracy's reduced to match observations)
    cvp_max = np.interp(d, self.pump_max['cvp']['d'], 
                           self.pump_max['cvp']['target']) * cfs_tafd
    swp_max = np.interp(d, self.pump_max['swp']['d'], 
                           self.pump_max['swp']['target']) * cfs_tafd
	
	#I think you've got the export/outflow relationship wrong
	#changed it to only consider exports made from unstored flows
	#consider export of stored flows separately, later
    if min_rule > gains: # additional flow needed to meet delta reqmnt
      self.dmin[t] = min_rule - gains
      #self.sodd_cvp[t] = cvp_max / export_ratio
      #self.sodd_swp[t] = swp_max / export_ratio
      self.sodd_cvp[t] = max(cvp_max/export_ratio - 0.75 * min_rule, 0)
      self.sodd_swp[t] = max(swp_max/export_ratio - 0.25 * min_rule, 0)
	  
    else: # extra unstored water available
      #self.sodd_cvp[t] = max((cvp_max - 0.55*(gains - min_rule)) / export_ratio, 0)
      #self.sodd_swp[t] = max((swp_max - 0.45*(gains - min_rule)) / export_ratio, 0)
      self.sodd_cvp[t] = max(cvp_max/export_ratio - 0.55 * gains,0)
      self.sodd_swp[t] = max(swp_max/export_ratio - 0.45 * gains,0)
    

  def step(self, t, cvp_flows, swp_flows):
    d = self.index.dayofyear[t]
    m = self.index.month[t]
    wyt = self.wyt[t]

    self.inflow[t] = self.gains[t] + cvp_flows + swp_flows

    #min_rule = np.interp(d, first_of_month, self.min_outflow[wyt]) * cfs_tafd
    #export_ratio = np.interp(d, first_of_month, self.export_ratio[wyt])
    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]
    
    cvp_max = np.interp(d, self.pump_max['cvp']['d'], 
                           self.pump_max['cvp']['pmax']) * cfs_tafd
    swp_max = np.interp(d, self.pump_max['swp']['d'], 
                           self.pump_max['swp']['pmax']) * cfs_tafd

    surplus = self.gains[t] - min_rule
    
    if surplus > 0:
      self.TRP_pump[t] = max(min((cvp_flows + 0.55 * self.gains[t]) * export_ratio, cvp_max),0)
      self.HRO_pump[t] = max(min((swp_flows + 0.45 * self.gains[t]) * export_ratio, swp_max),0)
    else:
      self.TRP_pump[t] = max(min((cvp_flows + 0.75 * self.gains[t]) * export_ratio, cvp_max),0)
      self.HRO_pump[t] = max(min((swp_flows + 0.25 * self.gains[t]) * export_ratio, swp_max),0)

    self.outflow[t] = self.inflow[t] - self.TRP_pump[t] - self.HRO_pump[t]
	
    # if outflow < max(min_rule, (1-export_ratio)*inflow):
    #   print('DELTA PROBLEM: %d' % t)
    #   print('OUTFLOW: %f, MIN: %f, PUMP: %f' % (outflow / cfs_tafd, min_rule / cfs_tafd, self.TRP_pump[t] / cfs_tafd))



  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['in','out','TRP_pump','HRO_pump']
    things = [self.inflow, self.outflow, self.TRP_pump, self.HRO_pump]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df