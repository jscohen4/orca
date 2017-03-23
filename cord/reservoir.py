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
    self.properties = json.load(open('cord/data/%s_properties.json' % key))
    self.K = self.properties['capacity']
    self.Q = df['%s_in'% key].values * cfs_tafd
    self.fci = df['%s_fci' % key].values
    self.S = np.zeros(T)
    self.R = np.zeros(T)
    self.Rtarget = np.zeros(T)
    self.S[0] = self.K / 1.5 # ?? assumption
    self.R[0] = 0

  def tocs(self,d,ix):
    t = self.properties['tocs']
    for i,v in enumerate(t['index']):
      if ix > v:
        break
    return np.interp(d, t['dowy'][i], t['storage'][i])

  def step(self, t):

    d = water_day(self.index.dayofyear[t])
    m = self.index.month[t]
    wyt = self.wyt[t]

    envmin = self.properties['env_min_flow'][wyt][m-1] * cfs_tafd
    toc = self.tocs(d, self.fci[t])


    # decide next release
    W = self.S[t-1] + self.Q[t]
    if W > toc:
      self.Rtarget[t] = 0.2*(W - toc)
    else:
      self.Rtarget[t] = envmin

    # then clip based on constraints
    self.R[t] = min(self.Rtarget[t], W)
    self.R[t] = min(self.R[t], self.properties['max_outflow'] * cfs_tafd)
    self.R[t] +=  max(W - self.R[t] - self.K, 0) # spill
    self.S[t] = W - self.R[t] # mass balance update


  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['storage', 'out', 'target']
    things = [self.S, self.R, self.Rtarget]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df



