from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import json
from util import *


class Reservoir():

  def __init__(self, df, key):
    T = len(df)
    self.index = df.index
    self.key = key
    self.wyt = df.SR_WYT
    for k,v in json.load(open('cord/data/%s_properties.json' % key)).items():
        setattr(self,k,v)

    self.Q = df['%s_in'% key].values * cfs_tafd
    self.fci = df['%s_fci' % key].values
    self.S = np.zeros(T)
    self.R = np.zeros(T)
    self.Rtarget = np.zeros(T)
    self.Rexport = np.zeros(T)
    self.S[0] = self.capacity / 1.5 # ?? assumption
    self.R[0] = 0

  def current_tocs(self,d,ix):
    for i,v in enumerate(self.tocs['index']):
      if ix > v:
        break
    return np.interp(d, self.tocs['dowy'][i], self.tocs['storage'][i])

  def step(self, t, dmin=0.0, sodd=0.0):

    d = self.index.dayofyear[t]
    dowy = water_day(d)
    m = self.index.month[t]
    wyt = self.wyt[t]

    envmin = self.env_min_flow[wyt][m-1] * cfs_tafd
    nodd = np.interp(d, first_of_month, self.nodd)
    sodd *= self.sodd_pct * self.sodd_curtail_pct[wyt]
    toc = self.current_tocs(dowy, self.fci[t])
    dout = dmin * self.delta_outflow_pct

    # decide next release
    W = self.S[t-1] + self.Q[t]
    self.Rtarget[t] = np.max((0.2*(W - toc), nodd+sodd+dout, envmin))
    self.Rexport[t] = sodd

    # then clip based on constraints
    self.R[t] = min(self.Rtarget[t], W - self.dead_pool)
    self.R[t] = min(self.R[t], self.max_outflow * cfs_tafd)
    self.R[t] +=  max(W - self.R[t] - self.capacity, 0) # spill
    self.S[t] = W - self.R[t] # mass balance update


  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['storage', 'out', 'target', 'rexport']
    things = [self.S, self.R, self.Rtarget, self.Rexport]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df



