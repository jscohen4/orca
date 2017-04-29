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

    for t in range(1,self.deltaT):
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
	
  def read_delta_inflow(self, datafile, sd = '10-01-1999'):
    self.deltadf = pd.read_csv(datafile, index_col = 0, parse_dates = True)
    self.delta.assign_flow(self.deltadf)
    self.deltaT = len(self.deltadf)
    for t in range(1,len(self.deltadf)):
      self.delta.find_gains(t,self.folsom.hist_releases[t-1], self.shasta.hist_releases[t-1], self.oroville.hist_releases[t-1])

  def find_running_WYI(self,zscore1, zscore2):
    self.fnf_SAC = self.df.BND_fnf.values * cfsd_mafd
    self.fnf_FEA = self.df.ORO_fnf.values * cfsd_mafd
    self.fnf_YUB = self.df.YRS_fnf.values * cfsd_mafd
    self.fnf_AME = self.df.FOL_fnf.values * cfsd_mafd
    self.snow_SAC = self.df.SHA_snowpack.values
    self.snow_FEA = self.df.ORO_snowpack.values
    self.snow_AME = self.df.FOL_snowpack.values
    apr_jul_fnf = np.zeros(self.df.index.year[len(self.df)-1] - self.df.index.year[0])
    oct_mar_fnf = np.zeros(self.df.index.year[len(self.df)-1] - self.df.index.year[0])
    remaining_flow = np.zeros((365,(self.df.index.year[len(self.df)-1]-self.df.index.year[0])))
    total_snowpack = np.zeros((365,(self.df.index.year[len(self.df)-1]-self.df.index.year[0])))
    self.forecastWYI_expt = np.zeros(len(self.df))
    self.forecastWYI_envr = np.zeros(len(self.df))

    current_year = 0
    for t in range(1,len(self.df)):
      d = int(self.df.index.dayofyear[t-1])	  
      dowy = water_day(d)
      if dowy == 0:
        current_year = current_year + 1
        oct_mar_fnf[current_year-1] = 0.0
        apr_jul_fnf[current_year-1] = 0.0
      if dowy < 182:
        oct_mar_fnf[current_year - 1] += self.fnf_SAC[t-1] + self.fnf_FEA[t-1] + self.fnf_YUB[t-1] + self.fnf_AME[t-1]
        total_snowpack[dowy][current_year-1] = self.snow_SAC[t-1] + self.snow_FEA[t-1] + self.snow_AME[t-1]
      elif dowy < 305:
        apr_jul_fnf[current_year - 1] += self.fnf_SAC[t-1] + self.fnf_FEA[t-1] + self.fnf_YUB[t-1] + self.fnf_AME[t-1]
        total_snowpack[dowy][current_year-1] = self.snow_SAC[t-1] + self.snow_FEA[t-1] + self.snow_AME[t-1]
      else:
        total_snowpack[dowy][current_year-1] = self.snow_SAC[t-1] + self.snow_FEA[t-1] + self.snow_AME[t-1]
    
    current_year = 0
    complete_year = 0	
    for t in range(1,len(self.df)):
      d = int(self.df.index.dayofyear[t-1])	  
      dowy = water_day(d)
      if dowy == 0:
        current_year = current_year + 1
        remaining_flow[dowy][current_year-1] = oct_mar_fnf[current_year-1]
      elif dowy < 182:
        remaining_flow[dowy][current_year-1] = remaining_flow[dowy-1][current_year-1] - (self.fnf_SAC[t-1] + self.fnf_FEA[t-1] + self.fnf_YUB[t-1] + self.fnf_AME[t-1])
      else:
        remaining_flow[dowy][current_year-1] = 0.0
      
      if dowy == 304:
        complete_year = complete_year + 1	  

    snowfall_std = np.zeros(365)
    earlyflow_std = np.zeros(365)
    earlyflow_mean = np.zeros(365)
    snowfallCoef = np.zeros((365,2))	
    pred_devi = np.zeros(complete_year)
    for x in range(1,365):
      one_year_flow = remaining_flow[x-1]
      one_year_snow = total_snowpack[x-1]
      coef = np.polyfit(one_year_snow[0:(complete_year-1)], apr_jul_fnf[0:(complete_year-1)],1)
      snowfallCoef[x-1][0] = coef[0]
      snowfallCoef[x-1][1] = coef[1]
      for y in range(1,complete_year):
        pred_devi[y-1] = apr_jul_fnf[y-1] - coef[0]*one_year_snow[y-1] - coef[1]
	  
      snowfall_std[x-1] = np.std(pred_devi)
      if x < 183:
        earlyflow_std[x-1] = np.std(one_year_flow)
      else:
        earlyflow_std[x-1] = 0.0

      earlyflow_mean[x-1] = np.mean(one_year_flow)
	
    prevValue = 10.0
    for t in range(1,len(self.df)):
      d = int(self.df.index.dayofyear[t-1])	  
      dowy = water_day(d)
      if dowy == 0:
        current_mar_oct = 0.0
        current_apr_jul = 0.0
        if t > 1:
          prevValue = min(self.forecastWYI_envr[t-2],10.0)
		  
      if dowy < 182:
        current_mar_oct += (self.fnf_SAC[t-1] + self.fnf_FEA[t-1] + self.fnf_YUB[t-1] + self.fnf_AME[t-1])
      elif dowy < 305:
        current_apr_jul += (self.fnf_SAC[t-1] + self.fnf_FEA[t-1] + self.fnf_YUB[t-1] + self.fnf_AME[t-1])
      
      snowpackEstimate = snowfallCoef[dowy][0]*(self.snow_SAC[t-1] + self.snow_FEA[t-1] + self.snow_AME[t-1]) + snowfallCoef[dowy][1]
      self.forecastWYI_envr[t-1] = 0.3 * prevValue + 0.3 * (current_mar_oct + earlyflow_std[dowy]*zscore1 + earlyflow_mean[dowy]) + 0.4 *( current_apr_jul + max(snowpackEstimate + snowfall_std[dowy]*zscore1 - current_apr_jul,0.0) )
      self.forecastWYI_expt[t-1] = 0.3 * prevValue + 0.3 * (current_mar_oct + earlyflow_std[dowy]*zscore2 + earlyflow_mean[dowy]) + 0.4 *( current_apr_jul + max(snowpackEstimate + snowfall_std[dowy]*zscore2 - current_apr_jul,0.0) )
    self.shasta.SRIforecast = self.forecastWYI_expt
    self.folsom.SRIforecast = self.forecastWYI_expt
    self.oroville.SRIforecast = self.forecastWYI_expt
    self.forecastWYT = []
    for x in range(1,len(self.df)):
      if self.forecastWYI_envr[x-1] <= 5.4:
        self.forecastWYT.append("C")
      elif self.forecastWYI_envr[x-1] <= 6.5:
        self.forecastWYT.append("D")
      elif self.forecastWYI_envr[x-1] <= 7.8:
        self.forecastWYT.append("BN")
      elif self.forecastWYI_envr[x-1] <= 9.2:
        self.forecastWYT.append("AN")
      else:
        self.forecastWYT.append("W")
    
    self.shasta.wyt = self.forecastWYT	
    self.oroville.wyt = self.forecastWYT
    self.folsom.wyt = self.forecastWYT
	