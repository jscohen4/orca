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
    self.wyt = df.WYT
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

  def step(self, t, d, m, wyt, cvp_flows, swp_flows, realinflow):

    self.inflow[t] = max(self.gains[t] + cvp_flows + swp_flows, 0) # realinflow * cfs_tafd

    min_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    export_ratio = self.export_ratio[wyt][m-1]

    cvp_max = self.cvp_pmax[d]
    swp_max = self.swp_pmax[d]
    
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


  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['in','out','TRP_pump','HRO_pump']
    things = [self.inflow, self.outflow, self.TRP_pump, self.HRO_pump]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df

