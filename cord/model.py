import numpy as np
import pandas as pd
from .reservoir import Reservoir
from .delta import Delta
from .util import *


class Model():

  def __init__(self, datafile, sd='10-01-1999'):
    self.df = pd.read_csv(datafile, index_col=0, parse_dates=True)
    self.T = len(self.df)

    self.shasta = Reservoir(self.df, 'SHA')
    self.folsom = Reservoir(self.df, 'FOL')
    self.oroville = Reservoir(self.df, 'ORO')
    self.delta = Delta(self.df, 'DEL')


  def simulate(self):

    for t in range(1,self.T):
      sumin = (self.shasta.Q[t] 
             + self.df.FOL_in[t] * cfs_tafd 
             + self.df.ORO_in[t] * cfs_tafd)
	  
      d = self.delta.index.dayofyear[t]
      dowy = water_day(d)
      wyt = self.delta.wyt[t]
      ##find the available storage in each resevoir - needs to be before delta.calc_flow_bounds (where actual releases are calculated)
      self.oroville.find_available_storage(t,self.oroville.exceedence[wyt])
      self.folsom.find_available_storage(t,self.folsom.exceedence[wyt])
      self.shasta.find_available_storage(t,self.shasta.exceedence[wyt])	  
      self.folsom.sodd_pct_var = self.delta.calc_flow_bounds(t, sumin, self.shasta.nodd, self.folsom.nodd, self.oroville.nodd, self.oroville.available_storage[t], self.folsom.available_storage[t], self.shasta.available_storage[t])
      self.shasta.sodd_pct_var = 1.0 - self.folsom.sodd_pct_var##set shasta/folsom contribution % (calculated together in delta class function)
      self.shasta.step(t, self.delta.dmin[t], self.delta.sodd_cvp[t])
      self.folsom.step(t, self.delta.dmin[t], self.delta.sodd_cvp[t])
      self.oroville.step(t, self.delta.dmin[t], self.delta.sodd_swp[t])
      self.delta.step(t, 
                      self.shasta.R_to_delta[t] + self.folsom.R_to_delta[t],
                      self.oroville.R_to_delta[t])

    return self.results_as_df()


  def results_as_df(self):
    df = pd.DataFrame(index=self.df.index)
    for x in [self.shasta, self.folsom, self.oroville, self.delta]:
      df = pd.concat([df, x.results_as_df(df.index)], axis=1)
    return df
	
  def create_flow_cdf(self, datafile, sd='10-01-1996'):
    self.pdf_data = pd.read_csv(datafile, index_col=0, parse_dates=True)
    self.shasta.find_release_func(self.pdf_data)
    self.oroville.find_release_func(self.pdf_data)
    self.folsom.find_release_func(self.pdf_data)
	
	