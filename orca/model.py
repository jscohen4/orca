import numpy as np
import pandas as pd
from .reservoir import Reservoir
from .delta import Delta
from .util import *
import matplotlib.pyplot as plt


class Model(): 

  def __init__(self, datafile, hist_datafile, sd='10-01-1999',projection = False, sim_gains = False):
    self.df = pd.read_csv(datafile, index_col=0, parse_dates=True)
    self.dfh = pd.read_csv(hist_datafile, index_col=0, parse_dates=True)
    self.sim_gains = sim_gains
    self.projection = projection
    self.T = len(self.df)
    self.shasta = Reservoir(self.df, self.dfh, 'SHA', self.projection)
    self.folsom = Reservoir(self.df, self.dfh, 'FOL', self.projection)
    self.oroville = Reservoir(self.df, self.dfh, 'ORO', self.projection)
    self.reservoirs = [self.shasta, self.folsom, self.oroville]
    self.delta = Delta(self.df, 'DEL', self.sim_gains)
    self.dayofyear = self.df.index.dayofyear
    self.month = self.df.index.month    
    self.wyt = self.df.WYT_sim

  def simulate(self):
    self.sumnodds = np.zeros(367)
    for d in range(0,365):
      for r in self.reservoirs:
        self.sumnodds[d] = sum([np.interp(d, first_of_month, r.nodd) for r in self.reservoirs])

    for t in range(1,self.T):
      d = self.dayofyear[t]
      dowy = water_day(d)
      m = self.month[t]
      wyt = self.wyt[t]
      self.oroville.find_available_storage(t,dowy,self.oroville.exceedence[wyt])
      self.folsom.find_available_storage(t,dowy,self.folsom.exceedence[wyt])
      self.shasta.find_available_storage(t,dowy,self.shasta.exceedence[wyt])   
      self.delta.calc_flow_bounds(t, d, m, wyt, dowy, self.sumnodds[d], self.oroville.available_storage[t], self.shasta.available_storage[t], self.folsom.available_storage[t])
      # self.shasta.sodd_pct = self.delta.shastaSODDPCT
      # self.folsom.sodd_pct = self.delta.folsomSODDPCT
      self.shasta.step(t, d, m, wyt, dowy, self.delta.dmin[t], self.delta.sodd_cvp[t], self.projection)
      self.folsom.step(t, d, m, wyt, dowy, self.delta.dmin[t], self.delta.sodd_cvp[t], self.projection)
      self.oroville.step(t, d, m, wyt, dowy, self.delta.dmin[t], self.delta.sodd_swp[t], self.projection)
      self.delta.step(t, d, m, wyt, dowy,
                      self.shasta.R_to_delta[t] + self.folsom.R_to_delta[t],
                      self.oroville.R_to_delta[t])
    return self.results_as_df()

  def results_as_df(self):
    df = pd.DataFrame(index=self.df.index)
    for x in [self.shasta, self.folsom, self.oroville, self.delta]:
      df = pd.concat([df, x.results_as_df(df.index)], axis=1)
    return df