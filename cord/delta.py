from __future__ import division
import numpy as np 
import pandas as pd
import json

tafw_cfsd = 1000 / 86400 * 43560 / 7
cfsd_tafw = 2.29568411*10**-5 * 86400 / 1000 * 7

class Model():

  def __init__(self, datafile, sd, ed):
    # should not technically be using any of this delta data
    self.df = pd.read_csv(datafile, index_col=0, parse_dates=True)[sd:ed]
    self.T = len(self.df.index)

  def simulate():
    delta = Delta('delta_rules.json')
    
    
    for t in range(1,T):
      v = delta.get_required_outflow(t, wyt)

class Delta():

  def __init__(self, datafile):
    self.rules = json.load(open(datafile))

  # datetime and wyt is all that's needed
  def get_required_outflow(t, wyt):
    pass
    # all constraints here


# blah blah blah
model = Model('data-clean.csv', sd='10-01-1999', ed='09-30-2015')
print model.df
print json.load(open('delta_rules.json'))