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
    self.reservoirs = [self.shasta, self.folsom, self.oroville]
    self.delta = Delta(self.df, 'DEL')


  def simulate(self):

    for t in range(1,self.T):
      self.delta.calc_flow_bounds(t, [r.nodd for r in self.reservoirs])
      self.shasta.step(t, self.delta.dmin[t], self.delta.sodd_cvp[t])
      self.folsom.step(t, self.delta.dmin[t], self.delta.sodd_cvp[t])
      self.oroville.step(t, self.delta.dmin[t], self.delta.sodd_swp[t])
      self.delta.step(t, 
                      self.shasta.R_to_delta[t] + self.folsom.R_to_delta[t],
                      self.oroville.R_to_delta[t], self.df.DeltaIn[t])

    return self.results_as_df()


  def results_as_df(self):
    df = pd.DataFrame(index=self.df.index)
    for x in [self.shasta, self.folsom, self.oroville, self.delta]:
      df = pd.concat([df, x.results_as_df(df.index)], axis=1)
    return df
