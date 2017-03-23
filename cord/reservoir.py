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
    self.wyt = df.SR_WYT
    self.properties = json.load(open('cord/data/%s_properties.json' % key))
    self.K = self.properties['capacity']
    self.S = np.zeros(T)
    self.R = np.zeros(T)
    self.Q = df['%s_in'% key].values * cfs_tafd
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
    toc = self.tocs(d, self.Q[t] + self.Q[t]*0.95)

    # envmin = 

    # mass balance
    self.S[t] = min(self.S[t-1] + self.Q[t-1] - self.R[t-1], self.K)

    # decide next release
    if self.S[t] + self.Q[t] > toc:
      self.R[t] = 0.2*(self.S[t] + self.Q[t] - toc)
    else:
      self.R[t] = envmin+10

