import numpy as np
import pandas as pd
from .reservoir import Reservoir
from .delta import Delta
from .district import District
from .util import *


class Model():

  def __init__(self, datafile, sd='10-01-1996'):
    self.df = pd.read_csv(datafile, index_col=0, parse_dates=True)
    self.index = self.df.index
    self.T = len(self.df)
	###Reservoir Initialization
    self.shasta = Reservoir(self.df, 'SHA')
    self.folsom = Reservoir(self.df, 'FOL')
    self.oroville = Reservoir(self.df, 'ORO')
    self.yuba = Reservoir(self.df,'YRS')
    self.newmelones = Reservoir(self.df,'NML')
    self.donpedro = Reservoir(self.df,'DNP')
    self.exchequer = Reservoir(self.df,'EXC')
    self.millerton = Reservoir(self.df,'MIL')
	##Delta Initialization
    self.delta = Delta(self.df, 'DEL')

  def simulate(self):
    
	###Daily Operations###
	##Step forward environmental parameters (snow & flow)
    ##Set Delta operating rules
    ##Water Balance on each reservoir
    ##Decisions - release water for delta export, flood control
    NMI = self.calc_wytypes(0,1)
    for t in range(1,self.T):
	  ##TIMESTEP
	  ##Step forward time, set date (for operating rules)
      d = self.index.dayofyear[t-1]
      m = int(self.index.month[t])
      dowy = water_day(d)
	  
      ##WATER YEAR TYPE CLASSIFICATION (for operating rules)
	  ##WYT uses flow forecasts - gets set every day, may want to decrease frequency (i.e. every month, season)
      NMI = self.calc_wytypes(t-1,dowy)
	  
	  ##REAL-WORLD RULE ADJUSTMENTS
	  ##Updates to reflect SJRR & Yuba Accords occuring during time period (1996-2016)
      self.update_regulations(t-1,dowy)
	  
	  ####NON-PROJECT USES
      ##Find out if reservoir releases need to be made for in-stream uses	  
      self.shasta.rights_call(self.shasta.downstream[t-1])
      self.oroville.rights_call(self.oroville.downstream[t-1])
      self.yuba.rights_call(self.yuba.downstream[t-1])
      self.folsom.rights_call(self.folsom.downstream[t-1])
      self.newmelones.rights_call(self.newmelones.downstream[t-1])
      self.donpedro.rights_call(self.donpedro.downstream[t-1])
      self.exchequer.rights_call(self.exchequer.downstream[t-1])
	  ##any additional losses before the delta inflow (.downstream member only accounts to downstream trib gauge)
	  ##must be made up by releases from Shasta and New Melones, respectively
      #self.newmelones.rights_call(self.delta.gains_sj[t-1],1)

      ##FIND MINIMUM ENVIRONMENTAL RELEASES
      self.exchequer.release_environmental(t-1,self.delta.forecastSJWYT)
      self.donpedro.release_environmental(t-1,self.delta.forecastSJWYT)
      self.newmelones.release_environmental(t-1,self.delta.forecastSJWYT)
      self.yuba.release_environmental(t-1,self.delta.forecastSCWYT)
      self.folsom.release_environmental(t-1,self.delta.forecastSCWYT)
      self.oroville.release_environmental(t-1,self.delta.forecastSCWYT)
      self.shasta.release_environmental(t-1,self.delta.forecastSCWYT)
	  
	  ##MINIMUM FLOW AT VERNALIS GAUGE(SAN JOAQUIN DELTA INFLOW)
      self.delta.vernalis_gains = self.newmelones.gains_to_delta + self.donpedro.gains_to_delta + self.exchequer.gains_to_delta + self.delta.gains_sj[t-1]
      self.delta.vernalis_gains += self.newmelones.envmin + self.donpedro.envmin + self.exchequer.envmin 
      self.exchequer.din, self.donpedro.din, self.newmelones.din  = self.delta.calc_vernalis_rule(t-1, NMI)
	 
	  ##MINIMUM FLOW AT RIO VIST GAUGE (SACRAMENTO DELTA INFLOW)
      self.delta.rio_gains = self.shasta.gains_to_delta + self.oroville.gains_to_delta + self.yuba.gains_to_delta + self.folsom.gains_to_delta + self.delta.gains_sac[t-1]
      self.delta.rio_gains += self.shasta.envmin + self.oroville.envmin + self.yuba.envmin + self.folsom.envmin
      self.shasta.din = self.delta.calc_rio_vista_rule(t-1)
	  
      ##MINIMUM DELTA OUTFLOW REQUIREMENTS
      self.delta.total_inflow = self.delta.eastside_streams[t-1] + self.delta.rio_gains + self.delta.vernalis_gains
      self.delta.cvp_stored_release = self.shasta.envmin + self.folsom.envmin + self.shasta.din
      self.delta.swp_stored_release = self.oroville.envmin + self.yuba.envmin
      cvp_dout, swp_dout = self.delta.calc_outflow_release(t-1)
	  	  
	  ###FLOW FORECASTING AT RESERVOIR w/ SNOWPACK
      ##future flow projections based on snow, current storage, and projected environmental releases
      self.oroville.find_available_storage(t-1,self.delta.forecastSCWYT)
      self.folsom.find_available_storage(t-1,self.delta.forecastSCWYT)
      self.shasta.find_available_storage(t-1,self.delta.forecastSCWYT)
      self.yuba.find_available_storage(t-1,self.delta.forecastSCWYT)
      shastaDelta, orovilleDelta, folsomDelta, yubaDelta = self.delta.assign_releases(self.shasta.available_storage[t-1],self.folsom.available_storage[t-1], self.oroville.available_storage[t-1], self.yuba.available_storage[t-1],self.shasta.S[t-1]-self.shasta.dead_pool, self.folsom.S[t-1]-self.folsom.dead_pool,self.oroville.S[t-1]-self.oroville.dead_pool,self.yuba.S[t-1]-self.yuba.dead_pool)	  
      self.shasta.dout = (cvp_dout + swp_dout)*shastaDelta
      self.folsom.dout = (cvp_dout + swp_dout)*folsomDelta
      self.oroville.dout = (cvp_dout + swp_dout)*orovilleDelta
      self.yuba.dout = (cvp_dout + swp_dout)*yubaDelta
      
	  ###DETERMINE RELEASES REQUIRED FOR DESIRED PUMPING
      ###Uses gains and environmental releases to determine additional releases required for
	  ###pumping (if desired), given inflow/export requirements, pump constraints, and CVP/SWP sharing of unstored flows
      cvp_available_storage = max(self.folsom.available_storage[t-1],0.0) + max(self.shasta.available_storage[t-1],0.0)
      swp_available_storage = max(self.oroville.available_storage[t-1],0.0) + max(self.yuba.available_storage[t-1],0.0)
      cvp_max, swp_max = self.delta.time_pumping(t-1,cvp_available_storage,swp_available_storage)
      self.delta.calc_flow_bounds(t-1, cvp_max, swp_max,cvp_available_storage,swp_available_storage)
      if cvp_available_storage > 0.0:
        self.shasta.sodd = self.delta.sodd_cvp[t-1]*max(self.shasta.available_storage[t-1],0.0)/cvp_available_storage
        self.folsom.sodd = self.delta.sodd_cvp[t-1]*max(self.folsom.available_storage[t-1],0.0)/cvp_available_storage
      else:
        self.shasta.sodd = 0.0
        self.folsom.sodd = 0.0
      if swp_available_storage > 0.0:
        self.oroville.sodd = self.delta.sodd_swp[t-1]*max(self.oroville.available_storage[t-1],0.0)/swp_available_storage
        self.yuba.sodd = self.delta.sodd_swp[t-1]*max(self.yuba.available_storage[t-1],0.0)/swp_available_storage
      else:
        self.oroville.sodd = 0.0
        self.yuba.sodd = 0.0
	        
	  ##SAN JOAQUIN RESERVOIR OPERATIONS
	  ##lower SJ basins - no 'release for exports' but used to meet delta targets @ vernalis
      ##Water Balance	  
      self.exchequer.step(t-1)
      self.donpedro.step(t-1)
      self.newmelones.step(t-1)	  

	  #SACRAMENTO RESERVOIR OPERATIONS
	  ##Water balance at each Northern Reservoir
      self.shasta.rights_call(self.delta.ccc[t-1]*-1.0,1)
      self.oroville.rights_call(self.delta.barkerslough[t-1]*-1.0,1)
      self.shasta.step(t-1)
      self.folsom.step(t-1)
      self.oroville.step(t-1)
      self.yuba.step(t-1)
	  
	  ###DELTA OPERATIONS
	  ##Given delta inflows (from gains and reservoir releases), find pumping
      cvp_stored_flow = self.shasta.R_to_delta[t-1] + self.folsom.R_to_delta[t-1]
      swp_stored_flow = self.oroville.R_to_delta[t-1] + self.yuba.R_to_delta[t-1]
      unstored_flow = self.delta.vernalis_gains + self.delta.eastside_streams[t-1] + self.shasta.gains_to_delta + self.oroville.gains_to_delta + self.yuba.gains_to_delta + self.folsom.gains_to_delta + self.delta.gains_sac[t-1]

      self.delta.step(t-1, cvp_stored_flow, swp_stored_flow, unstored_flow)
	  
    return self.results_as_df()
	
  def update_regulations(self,t,dowy):
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    ##San Joaquin River Restoration Project, started in October of 2009 (WY 2009)
	##Additional Releases from Millerton Lake depending on WYT
    if y == 2009 and m >= 10:
      self.millerton.sjrr_release = self.millerton.sj_riv_res_flows(t, dowy)
    elif y > 2009:
      self.millerton.sjrr_release = self.millerton.sj_riv_res_flows(t, dowy)
    

	##Yuba River Accord, started in Jan of 2006 (repaces minimum flow requirements)
    if y >= 2006:
      self.yuba.env_min_flow = self.yuba.env_min_flow_ya
      self.yuba.temp_releases = self.yuba.temp_releases_ya
      self.yuba.carryover_target['W'] = 650
      self.yuba.carryover_target['AN'] = 650
      self.yuba.carryover_target['BN'] = 650
      self.yuba.carryover_target['D'] = 650
      self.yuba.carryover_target['C'] = 650
	
  

  def results_as_df(self):
    df = pd.DataFrame(index=self.df.index)
    for x in [self.shasta, self.folsom, self.oroville, self.yuba, self.newmelones, self.donpedro, self.exchequer, self.delta]:
      df = pd.concat([df, x.results_as_df(df.index)], axis=1)
    return df
	
  def create_flow_cdf(self):
  
    expected_outflow_req = self.delta.calc_expected_delta_outflow(self.shasta.downstream,self.oroville.downstream,self.yuba.downstream,self.folsom.downstream, self.shasta.temp_releases, self.oroville.temp_releases, self.yuba.temp_releases, self.folsom.temp_releases)
    self.shasta.find_release_func()
    self.oroville.find_release_func()
    self.folsom.find_release_func()
    self.yuba.find_release_func()
    self.newmelones.find_release_func()
    self.donpedro.find_release_func()
    self.exchequer.find_release_func()
    self.millerton.find_release_func()
		
    self.shasta.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)##what do they expect to need to release for env. requirements through the end of september
    self.oroville.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)
    self.folsom.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)
    self.yuba.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)
    self.newmelones.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)
    self.donpedro.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)
    self.exchequer.calc_expected_min_release(expected_outflow_req,self.delta.gains_sac)
      
  def find_running_WYI(self):
    ###Pre-processing function
	##Finds the 8 River, Sacramento, and San Joaquin indicies based on flow projections
    lastYearSRI = 10.26 # WY 1996
    lastYearSJI = 4.12 # WY 1996
    current_year = 0
    startMonth = int(self.index.month[0])
    startYear = int(self.index.year[0])
    rainflood_sac_obs = 0.0
    snowflood_sac_obs = 0.0
    rainflood_sj_obs = 0.0
    snowflood_sj_obs = 0.0
    index_exceedence = 2

    for t in range(1,self.T):
      d = int(self.index.dayofyear[t-1])
      dowy = water_day(d)
      m = int(self.index.month[t-1])
      y = int(self.index.year[t-1])
      da = int(self.index.day[t-1])
      index_exceedence_sac = 9
      index_exceedence_sj = 5
	  ##8 River Index
      self.delta.eri[m-startMonth + (y - startYear)*12] += (self.shasta.fnf[t] + self.folsom.fnf[t-1] + self.oroville.fnf[t-1] + self.yuba.fnf[t-1] + self.newmelones.fnf[t-1] + self.donpedro.fnf[t-1] + self.exchequer.fnf[t-1] + self.millerton.fnf[t-1])*1000.0
	  ####################Sacramento Index#############################################################################################
	  ##Individual Rainflood Forecast - either the 90% exceedence level prediction, or the observed WYTD fnf value
      if m >=10:
        self.delta.forecastSJI[t-1] = lastYearSJI
        self.delta.forecastSRI[t-1] = lastYearSRI
      else:
        shasta_rain_forecast = self.shasta.rainflood_fnf[t] + self.shasta.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
        oroville_rain_forecast = self.oroville.rainflood_fnf[t] + self.oroville.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
        folsom_rain_forecast = self.folsom.rainflood_fnf[t] + self.folsom.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
        yuba_rain_forecast = self.yuba.rainflood_fnf[t] + self.yuba.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
	  ##
	  ##SAC TOTAL RAIN
        sac_rain = max(rainflood_sac_obs, shasta_rain_forecast + oroville_rain_forecast + folsom_rain_forecast + yuba_rain_forecast)
	  ##
	  ##Individual Snowflood Forecast - either the 90% exceedence level prediction, or the observed WYTD fnf value
        shasta_snow_forecast = self.shasta.snowflood_fnf[t] + self.shasta.snowfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
        oroville_snow_forecast = self.oroville.snowflood_fnf[t] + self.oroville.snowfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
        folsom_snow_forecast = self.folsom.snowflood_fnf[t] + self.folsom.snowfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
        yuba_snow_forecast = self.yuba.snowflood_fnf[t] + self.yuba.snowfnf_stds[dowy]*z_table_transform[index_exceedence_sac]
      ##
	  ##SAC TOTAL SNOW
        sac_snow = max(snowflood_sac_obs, shasta_snow_forecast + oroville_snow_forecast + folsom_snow_forecast + yuba_snow_forecast)
      ##
	  #######################################################################################################################################
	  #####################San Joaquin Index################################################################################################
      ##Individual Rainflood Forecast - either the 90% exceedence level prediction, or the observed WYTD fnf value
        nm_rain_forecast = self.newmelones.rainflood_fnf[t] + self.newmelones.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sj]
        dp_rain_forecast = self.donpedro.rainflood_fnf[t] + self.donpedro.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sj]
        exc_rain_forecast = self.exchequer.rainflood_fnf[t] + self.exchequer.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sj]
        mill_rain_forecast = self.millerton.rainflood_fnf[t] + self.millerton.rainfnf_stds[dowy]*z_table_transform[index_exceedence_sj]
	  ##
	  ##SJ TOTAL RAIN
        sj_rain = max(rainflood_sj_obs, nm_rain_forecast + dp_rain_forecast + exc_rain_forecast + mill_rain_forecast)
	  ##
	  ##Individual Snowflood Forecast - either the 90% exceedence level prediction, or the observed WYTD fnf value
        nm_snow_forecast = self.newmelones.snowflood_fnf[t] + self.newmelones.snowfnf_stds[dowy]*z_table_transform[index_exceedence]
        dp_snow_forecast = self.donpedro.snowflood_fnf[t] + self.donpedro.snowfnf_stds[dowy]*z_table_transform[index_exceedence]
        exc_snow_forecast = self.exchequer.snowflood_fnf[t] + self.exchequer.snowfnf_stds[dowy]*z_table_transform[index_exceedence]
        mill_snow_forecast = self.millerton.snowflood_fnf[t] + self.millerton.snowfnf_stds[dowy]*z_table_transform[index_exceedence]
      ##
	  ##SJ TOTAL SNOW
        sj_snow = max(snowflood_sj_obs, nm_snow_forecast + dp_snow_forecast + exc_snow_forecast + mill_snow_forecast)
      ###INDEX FORECASTS########################################################################################################################
        self.delta.forecastSJI[t-1] = min(lastYearSJI,4.5)*0.2 + sj_rain*0.2 + sj_snow*0.6
        self.delta.forecastSRI[t-1] = min(lastYearSRI,10)*0.3 + sac_rain*0.3 + sac_snow*0.4

      ##REAL-TIME OBSERVATIONS
      if m >= 10 or m <= 3:
        rainflood_sac_obs += self.shasta.fnf[t] + self.oroville.fnf[t] + self.folsom.fnf[t] + self.yuba.fnf[t]
        rainflood_sj_obs += self.newmelones.fnf[t] + self.donpedro.fnf[t] + self.exchequer.fnf[t] + self.millerton.fnf[t]
      else:
        snowflood_sac_obs += self.shasta.fnf[t] + self.oroville.fnf[t] + self.folsom.fnf[t] + self.yuba.fnf[t]
        snowflood_sj_obs += self.newmelones.fnf[t] + self.donpedro.fnf[t] + self.exchequer.fnf[t] + self.millerton.fnf[t]
		
	##SAVE INDEX FROM EACH YEAR (FOR USE IN NEXT YEAR'S FORECAST	
      if m == 9 and da == 30:
        lastYearSRI = 0.3*min(lastYearSRI,10) + 0.3*rainflood_sac_obs + 0.4*snowflood_sac_obs
        lastYearSJI = 0.2*min(lastYearSJI,4.5) + 0.2*rainflood_sj_obs + 0.6*snowflood_sj_obs
        rainflood_sac_obs = 0.0
        snowflood_sac_obs = 0.0
        rainflood_sj_obs = 0.0
        snowflood_sj_obs = 0.0
 
  def calc_wytypes(self,t,dowy):
  
####NOTE:  Full natural flow data is in MAF, inflow data is in TAF  
##Index for Shasta Min Flows
############################
    if self.delta.forecastSRI[t] <= 5.4:
      self.shasta.forecastWYT = "C"
      self.delta.forecastSCWYT = "C"
    elif self.delta.forecastSRI[t] <= 6.6:
      self.shasta.forecastWYT = "D"
      self.delta.forecastSCWYT = "D" 
    elif self.delta.forecastSRI[t] <= 7.8:
      self.shasta.forecastWYT = "BN"
      self.delta.forecastSCWYT = "BN"
    elif self.delta.forecastSRI[t] <= 9.2:
      self.shasta.forecastWYT = "AN"
      self.delta.forecastSCWYT = "AN" 
    else:
      self.shasta.forecastWYT = "W"
      self.delta.forecastSCWYT = "W"

##Index for Oroville Min Flows
############################	  
    if self.oroville.snowflood_fnf[t] < 0.55*1.942:
      self.oroville.forecastWYT = "D"
    else:
      self.oroville.forecastWYT = "W"
    
    if self.delta.forecastSRI[t] <= 5.4:
      self.oroville.forecastWYT = "C"
  
##Index for Yuba Min Flows
############################	
    eos_date = t - dowy
    if eos_date < 0:
      eos_date = 0
	  
    yubaIndex = (self.yuba.rainflood_fnf[t] + self.yuba.snowflood_fnf[t])*1000 + self.yuba.S[eos_date] - 234.0
    if yubaIndex >= 1400:
      self.yuba.forecastWYT = "W" 
    elif yubaIndex >= 1040:
      self.yuba.forecastWYT = "AN"
    elif yubaIndex >= 920:
      self.yuba.forecastWYT = "BN"
    elif yubaIndex >= 820:
      self.yuba.forecastWYT = "D"
    elif yubaIndex >= 693:
      self.yuba.forecastWYT = "C"
    else:
      self.yuba.forecastWYT = "EC"
  
##Index for Folsom Min Flows
############################
##Folsom has the most ridiculous operating rules, and combines a bunch of different 'indicies' throughout the year to determine min flows	
    if dowy < 91:
      folsomIndex = self.folsom.S[eos_date] + (361.701 - self.folsom.fci[eos_date])
      if folsomIndex >= 848:
        self.folsom.forecastWYT = "W"
      elif folsomIndex >= 746:
        self.folsom.forecastWYT = "AN" 
      elif folsomIndex >= 600:
        self.folsom.forecastWYT = "BN" 
      elif folsomIndex >= 300:
        self.folsom.forecastWYT = "D"
      else:
        self.folsom.forecastWYT = "C"
    elif dowy < 150:
      folsomIndex = self.folsom.S[eos_date] + (361.701 - self.folsom.fci[eos_date])
      if self.delta.forecastSRI[t] <= 5.4 and folsomIndex < 600:
        self.folsom.forecastWYT = "C"
      elif self.delta.forecastSRI[t] <= 5.4 and folsomIndex < 746:
        self.folsom.forecastWYT = "D"
      elif self.delta.forecastSRI[t] <= 5.4 and folsomIndex < 848:
        self.folsom.forecastWYT = "AN"
      elif self.delta.forecastSRI[t] < 7.8 and folsomIndex < 600:
        self.folsom.forecastWYT = "D"
      elif self.delta.forecastSRI[t] < 7.8 and folsomIndex < 746:
        self.folsom.forecastWYT = "BN"
      else:
        self.folsom.forecastWYT = "W"
    else:
      folsomIndex = (self.folsom.snowflood_fnf[t] - sum(self.folsom.fnf[(t-dowy+181):(t-dowy+211)]) + sum(self.folsom.fnf[(t-dowy+304):(t-dowy+364)]))*1000
      if folsomIndex < 250:
        self.folsom.forecastWYT = "C"
      elif folsomIndex < 375:
        self.folsom.forecastWYT = "D"
      elif folsomIndex < 460:
        self.folsom.forecastWYT = "BN"
      elif folsomIndex < 550:
        self.folsom.forecastWYT = "AN"
      else:
        self.folsom.forecastWYT = "W"
  
##Index for New Melones Min Flows
############################
    if dowy <= 150:
      eof_storage = t - dowy - 215
      if eof_storage < 0:
        eof_storage == 0
      newmelonesIndex = self.newmelones.S[eof_storage] + sum(self.newmelones.fnf[(eof_storage+1):(t-dowy)])*1000
    else:
      eof_storage = t - dowy + 149
      newmelonesIndex = self.newmelones.S[eof_storage] + (sum(self.newmelones.fnf[(eof_storage+1):(t-dowy+181)]) + self.newmelones.snowflood_fnf[t] + sum(self.newmelones.fnf[(t-dowy+304):(t-dowy+365)]))*1000
	
    if newmelonesIndex < 1400:
      self.newmelones.forecastWYT = "C"
    elif newmelonesIndex < 2000:
      self.newmelones.forecastWYT = "D"
    elif newmelonesIndex < 2500:
      self.newmelones.forecastWYT = "BN"
    elif newmelonesIndex < 3000:
      self.newmelones.forecastWYT = "AN"
    else:
      self.newmelones.forecastWYT = "W"
  
##Index for Don Pedro Min Flows
############################
    if self.delta.forecastSJI[t] <= 2.1:
      self.donpedro.forecastWYT = "C"
      self.delta.forecastSJWYT = "C"
    elif self.delta.forecastSJI[t] <= 2.5:
      self.donpedro.forecastWYT = "D"
      self.delta.forecastSJWYT = "D"
    elif self.delta.forecastSJI[t] <= 3.1:
      self.donpedro.forecastWYT = "BN"
      self.delta.forecastSJWYT = "BN"
    elif self.delta.forecastSJI[t] <= 3.8:
      self.donpedro.forecastWYT = "AN"
      self.delta.forecastSJWYT = "AN"
    else:
      self.donpedro.forecastWYT = "W"
      self.delta.forecastSJWYT = "W"
    
  
##Index for Exchequer Min Flows
############################	  
    if self.exchequer.snowflood_fnf[t] < .45:
      self.exchequer.forecastWYT = "D"
    else:
      self.exchequer.forecastWYT = "AN"
  	
    return newmelonesIndex
