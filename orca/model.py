import numpy as np
import pandas as pd
from .reservoir import Reservoir
from .delta import Delta
from .util import *
import matplotlib.pyplot as plt


class Model(): 

  def __init__(self, datafile, hist_datafile, SHA_shift, ORO_shift, FOL_shift, sd='10-01-1999',projection = False, sim_gains = False):
    self.df = pd.read_csv(datafile, index_col=0, parse_dates=True)
    self.dfh = pd.read_csv(hist_datafile, index_col=0, parse_dates=True)
    self.sim_gains = sim_gains
    self.projection = projection
    self.T = len(self.df)
    self.shasta = Reservoir(self.df, self.dfh, 'SHA', SHA_shift, self.projection)
    self.folsom = Reservoir(self.df, self.dfh, 'FOL', FOL_shift, self.projection)
    self.oroville = Reservoir(self.df, self.dfh, 'ORO', ORO_shift, self.projection)
    self.reservoirs = [self.shasta, self.folsom, self.oroville]
    self.delta = Delta(self.df, 'DEL', self.sim_gains)
    self.dayofyear = self.df.index.dayofyear
    self.month = self.df.index.month    
    self.wyt = self.df.WYT_sim

  def simulate(self):
    self.sumnodds = np.zeros(367)
    for d in range(0,366):
      for r in self.reservoirs:
        self.sumnodds[d] = sum([np.interp(d, first_of_month, r.nodd) for r in self.reservoirs])
    leap_year = False
    for t in range(1,self.T):
      d = self.dayofyear[t]
      if d == 0:
        leap_year = False
      dowy = water_day(d)
      if leap_year == True: 
        dowy +=1
      if dowy == 92 and doy == 92:
        leap_year = True
      doy = dowy
      m = self.month[t]
      wyt = self.wyt[t]
      self.oroville.find_available_storage(t,d,dowy)#,self.oroville.exceedance[wyt])
      self.folsom.find_available_storage(t,d,dowy)#,self.folsom.exceedance[wyt])
      self.shasta.find_available_storage(t,d,dowy)#,self.shasta.exceedance[wyt])   
      self.delta.calc_flow_bounds(t, d, m, wyt, dowy, self.sumnodds[d], self.oroville.available_storage[t], self.shasta.available_storage[t], self.folsom.available_storage[t])
      self.shasta.sodd_pct = self.delta.shastaSODDPCT
      self.folsom.sodd_pct = self.delta.folsomSODDPCT
      self.shasta.step(t, d, m, wyt, dowy, self.delta.dmin[t], self.delta.sodd_cvp[t], self.projection)
      self.folsom.step(t, d, m, wyt, dowy, self.delta.dmin[t], self.delta.sodd_cvp[t], self.projection)
      self.oroville.step(t, d, m, wyt, dowy, self.delta.dmin[t], self.delta.sodd_swp[t], self.projection)
      self.delta.step(t, d, m, wyt, dowy,
                      self.shasta.R_to_delta[t] + self.folsom.R_to_delta[t],
                      self.oroville.R_to_delta[t],self.oroville.available_storage[t], self.shasta.available_storage[t], self.folsom.available_storage[t],self.sumnodds[d])
    return self.results_as_df()

  def results_as_df(self):
    df = pd.DataFrame(index=self.df.index)
    for x in [self.shasta, self.folsom, self.oroville, self.delta]:
      df = pd.concat([df, x.results_as_df(df.index)], axis=1)
    return df