from __future__ import division
import numpy as np 
import pandas as pd
import json
from util import *

class Delta():

  def __init__(self, df, key):
    T = len(df)
    self.index = df.index
    self.key = key
    self.wyt = df.SR_WYT

    for k,v in json.load(open('cord/data/Delta_properties.json')).items():
      setattr(self,k,v)

    # what vars to store/save here
    self.dmin = np.zeros(T)
    self.sodd_cvp = np.zeros(T)
    self.TRP_pump = np.zeros(T)
    self.HRO_pump = np.zeros(T)

  def calc_flow_bounds(self, t, sum_res_inflows):
    d = self.index.dayofyear[t]
    m = self.index.month[t]
    wyt = self.wyt[t]

    gains = sum_res_inflows * self.gains_factor[m-1]
    min_rule = np.interp(d, first_of_month, self.min_outflow[wyt]) * cfs_tafd
    export_ratio = np.interp(d, first_of_month, self.export_ratio[wyt])
    outflow_ratio = 1 - export_ratio
    # ratio_rule = gains * 
    

    # pump capacities
    cvp_max = np.interp(d, self.pump_max['cvp']['d'], 
                           self.pump_max['cvp']['pmax']) * cfs_tafd

    if min_rule > gains: # additional flow needed
      self.dmin[t] = min_rule - gains
      self.sodd_cvp[t] = cvp_max / export_ratio
    else: # extra unstored water available
      self.sodd_cvp[t] = max((cvp_max - 0.55*(gains - min_rule)) / export_ratio, 0)
      # self.TRP_pump[t] = 0.55*(gains - min_rule)

    # cvp_pmax = 
    # swp_pmax = 


