from __future__ import division
import numpy as np 
import pandas as pd
import json
from .util import *

class Delta():

  def __init__(self, df, key):
    self.T = len(df)
    self.index = df.index
    self.key = key
    self.forecastSCWYT = "AN"
    self.forecastSJWYT = "AN"
    self.last_year_vamp = 5.0

    for k,v in json.load(open('cord/data/Delta_properties.json')).items():
      setattr(self,k,v)

    # Vectors for delta Inflows
    self.gains = np.zeros(self.T)
    self.gains_sac = df.SAC_gains * cfs_tafd
    self.gains_sj = df.SJ_gains * cfs_tafd
    self.depletions = df.delta_depletions * cfs_tafd
    self.ccc = df.CCC_pump * cfs_tafd
    self.barkerslough = df.BRK_pump *cfs_tafd
    self.vernalis_flow = np.zeros(self.T)
    self.eastside_streams = df.EAST_gains * cfs_tafd
    self.inflow = np.zeros(self.T)

	##Vectors for delta outflows/exports
    self.dmin = np.zeros(self.T)
    self.sodd_cvp = np.zeros(self.T)
    self.sodd_swp = np.zeros(self.T)
    self.TRP_pump = np.zeros(self.T)
    self.HRO_pump = np.zeros(self.T)
    self.outflow = np.zeros(self.T)
    self.x2 = np.zeros(self.T)
    self.x2[0] = 82.0
    self.x2constraint = {}
    self.x2constraint['W'] = np.zeros(366)
    self.x2constraint['AN'] = np.zeros(366)
    self.x2constraint['BN'] = np.zeros(366)
    self.x2constraint['D'] = np.zeros(366)
    self.x2constraint['C'] = np.zeros(366)

	##River Indicies
    self.eri = np.zeros(self.T)	
    self.forecastSRI = np.zeros(self.T)
    self.forecastSJI = np.zeros(self.T)

	##Old/Middle River Calculations
    self.hist_OMR = df.OMR * cfs_tafd
    self.hist_TRP_pump = df.TRP_pump * cfs_tafd
    self.hist_HRO_pump = df.HRO_pump * cfs_tafd
    self.OMR = np.zeros(self.T)
	
	##Variables for determining releases for export (initialize)
    self.cvp_aval_stor = 0.5
    self.swp_aval_stor = 0.5
    self.cvp_delta_outflow_pct = 0.75
    self.swp_delta_outflow_pct = 0.25
	
    self.swp_allocation = np.zeros(self.T)
    self.cvp_allocation = np.zeros(self.T)
    self.annual_HRO_pump = np.zeros(int(self.index.year[self.T-1]) - int(self.index.year[0]))
    self.annual_TRP_pump = np.zeros(int(self.index.year[self.T-1]) - int(self.index.year[0]))


	
  def calc_expected_delta_outflow(self,shastaD,orovilleD,yubaD,folsomD,shastaMIN,orovilleMIN,yubaMIN,folsomMIN):
    expected_outflow_releases = {}
    expected_outflow_releases['W'] = np.zeros(366)
    expected_outflow_releases['AN'] = np.zeros(366)
    expected_outflow_releases['BN'] = np.zeros(366)
    expected_outflow_releases['D'] = np.zeros(366)
    expected_outflow_releases['C'] = np.zeros(366)

	
    num_obs = np.zeros(366)
    for t in range(1,self.T):
      d = int(self.index.dayofyear[t-1])
      dowy = water_day(d)
      m = int(self.index.month[t])
      zone = int(np.interp(dowy, self.san_joaquin_min_flow['d'], self.san_joaquin_min_flow['zone']))
      for wyt in ['W', 'AN', 'BN', 'D', 'C']:
        outflow_rule = self.min_outflow[wyt][m-1] * cfs_tafd
        vernalis_flows = max(self.gains_sj[t],self.san_joaquin_min_flow[wyt][zone-1]* cfs_tafd)
        trib_flow = max(shastaMIN[wyt][m-1]*cfs_tafd,shastaD[t]) + max(orovilleMIN[wyt][m-1]*cfs_tafd,orovilleD[t]) + max(yubaMIN[wyt][m-1]*cfs_tafd,yubaD[t]) + max(folsomMIN[wyt][m-1]*cfs_tafd,folsomD[t])
        expected_outflow_releases[wyt][dowy] += max(outflow_rule - self.gains_sac[t] - trib_flow - self.eastside_streams[t] - vernalis_flows - self.depletions[t], 0.0)
	  
      num_obs[dowy] += 1.0
    	
    for x in range(1,365):
      self.x2constraint['C'][x] = 90.0
      for wyt in ['W', 'AN', 'BN', 'D']:
        if x > 180 and x < 274:
          self.x2constraint[wyt][x] = 74.0 + 5.0*(x-181)/94
        elif x > 274 and x < 318:
          self.x2constraint[wyt][x] = 79.0 + 6.0*(x-275)/44
        else:
          self.x2constraint[wyt][x] = 90.0
		
      for wyt in ['W', 'AN', 'BN', 'D', 'C']:
        expected_outflow_releases[wyt][x] = expected_outflow_releases[wyt][x]/num_obs[x]
		
    return expected_outflow_releases
  
  def calc_rio_vista_rule(self,t):
    m = int(self.index.month[t])
    wyt = self.forecastSCWYT
    shasta_contr = max(self.rio_vista_min[wyt][m-1]*cfs_tafd - self.rio_gains,0.0)
    self.rio_gains += shasta_contr
    return shasta_contr
	
  def calc_vernalis_rule(self,t,NMI):
    ##Calculates delta rules at Vernalis (San Joaquin/delta confluence)
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    wyt = self.forecastSJWYT
    
	##D_1641 RULES
    ##zone refers to one of four 'min flow' groupings throughout the year, based on WYT
	##Note: I forget why I wrote the Delta_property files like this, with the 'zones' but it works
    zone = int(np.interp(dowy, self.san_joaquin_min_flow['d'], self.san_joaquin_min_flow['zone']))
    d_1641_min = self.san_joaquin_min_flow[wyt][zone-1] * cfs_tafd	
    d_1641_flow = max(d_1641_min - self.vernalis_gains, 0.0)

    ##Exchequer and Don Pedro only make contributions during a given time of the year
    if m == 9 and d == 1:
      self.last_year_vamp = self.get_vamp_no(self.forecastSJWYT)
        
    if zone == 4:
      if y < 2009:
        vamp_pulse = max(self.vamp_rule(self.vernalis_gains) - self.vernalis_gains,0.0)
        merced_contr = vamp_pulse*0.5
        tuolumne_contr = vamp_pulse*0.25
        stanislaus_contr = vamp_pulse*0.25
      else:
        vamp_pulse = max(self.new_vamp_rule[wyt]*cfs_tafd - self.vernalis_gains,0.0)
        merced_contr = vamp_pulse*0.5
        tuolumne_contr = vamp_pulse*0.25
        stanislaus_contr = vamp_pulse*0.25     
    else:
      vamp_pulse = 0.0
      stanislaus_contr = d_1641_flow
      merced_contr = 0.0
      tuolumne_contr = 0.0
    
	##BIOPS RULES for Vernalis start in 2009, before then only d_1641 rule is used
    if y > 2009:
      #biops_min = np.interp(dowy,self.san_joaquin_min['biops_d'],self.san_joaquin_min['biops_on_off'])*np.interp(NMI,self.san_joaquin_min['biops_NMI'],self.san_joaquin_min['biops_flow']) * cfs_tafd
      biops_min = 0.0
    else:
      biops_min = 0.0
	
    ##BIOPS Releases are only made from New Melones (CVP reservoir)
    if biops_min > d_1641_min:
      stanislaus_contr += max(biops_min-self.vernalis_gains,0.0) - max(max(d_1641_min,vamp_pulse)-self.vernalis_gains,0.0)
	
	###Releases for vernalis are added to gains at the delta
    self.vernalis_gains += merced_contr + tuolumne_contr + stanislaus_contr

    return merced_contr, tuolumne_contr, stanislaus_contr
  
  def vamp_rule(self,ix):
  ##Interpolates rules from tocs_rule in *_properties.json file to get the top of the conservation
  ##pool in order to determine flood control releases in reservoir.step
    flow = ix*tafd_cfs
    for i,v in enumerate(self.vamp_flows['forecast']):
      if flow < v:
        break
    ##Double-step or dry-year reduction: if its been wet, flow targets are
	##taken from the next-highest 'bracket' (double-stepping the brackets), if
	##its been dry, no VAMP releases are required
    this_year_no = self.get_vamp_no(self.forecastSJWYT)
    step_number = self.last_year_vamp + this_year_no

    if step_number > 6.9:
      vamp_target = self.vamp_flows['target'][i+1]*cfs_tafd##double step
    elif step_number < 4.1:
      vamp_target = 0.0##dry-year reduction
    else:
      vamp_target = self.vamp_flows['target'][i]*cfs_tafd#normal conditions
		
    return vamp_target
    
  def get_vamp_no(self,wyt):
    if wyt == 'W':
      vamp_no = 5.0
    elif wyt == 'AN':
      vamp_no = 4.0
    elif wyt == 'BN':
      vamp_no = 3.0
    elif wyt == 'D':
      vamp_no = 2.0
    elif wyt == 'C':
      vamp_no = 1.0
    return vamp_no
	
  def calc_outflow_release(self,t):
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    wyt = self.forecastSCWYT
    outflow_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    if dowy > 180 and dowy < 318:
      if self.x2[t-1] > self.x2constraint[wyt][dowy]:
        x2outflow = 10**((self.x2constraint[wyt][dowy] - 10.16 - 0.945*self.x2[t])/(-1.487))
      else:
        x2outflow = 0.0
    else:
      x2outflow = 0.0
		
    salinity_rule = min(x2outflow*cfs_tafd,12000.0*cfs_tafd)
	  
    min_rule = max(outflow_rule, salinity_rule)
    dout = max(min_rule - self.depletions[t] - (self.total_inflow - self.cvp_stored_release - self.swp_stored_release),0.0)
    
    if dout > (self.cvp_stored_release + self.swp_stored_release):
      swp_dout = 0.25*(dout - self.cvp_stored_release - self.swp_stored_release)
      cvp_dout = 0.75*(dout - self.cvp_stored_release - self.swp_stored_release)
      self.cvp_outflow = self.cvp_stored_release + cvp_dout
      self.swp_outflow = self.swp_stored_release + swp_dout
    else:
      cvp_dout = 0.0
      swp_dout = 0.0
      if self.cvp_stored_release > 0.75*dout:
        self.cvp_outflow = dout - self.swp_stored_release
        self.swp_outflow = self.swp_stored_release
      else:
        self.cvp_outflow = self.cvp_stored_release
        self.swp_outflow = dout - self.cvp_stored_release
	  
    self.cvp_stored_release += cvp_dout
    self.swp_stored_release += swp_dout
    self.total_inflow += (cvp_dout + swp_dout)
    return cvp_dout, swp_dout   

  def assign_releases(self,shastaAS,folsomAS,orovilleAS,yubaAS,shastaODP,folsomODP,orovilleODP,yubaODP):
    cvp_AS = max(shastaAS,0.0) + max(folsomAS,0.0)
    swp_AS = max(yubaAS,0.0) + max(orovilleAS,0.0)
    cvp_ODP = max(shastaODP,0.0) + max(folsomODP,0.0)
    swp_ODP = max(yubaODP,0.0) + max(orovilleODP,0.0)
    total_ODP = cvp_ODP + swp_ODP
    if cvp_ODP and swp_ODP > 0.0:
      if swp_AS > 0.0:
        orovilleFrac = 0.25*max(orovilleODP,0.0)/swp_ODP
        yubaFrac = 0.25*max(yubaODP,0.0)/swp_ODP
      else:
        orovilleFrac = 0.25*max(orovilleODP,0.0)/swp_ODP
        yubaFrac = 0.25*max(yubaODP,0.0)/swp_ODP
 
      if cvp_AS > 0.0:
        shastaFrac = 0.75*max(shastaODP,0.0)/cvp_ODP
        folsomFrac = 0.75*max(folsomODP,0.0)/cvp_ODP
      else:
        shastaFrac = 0.75*max(shastaODP,0.0)/cvp_ODP
        folsomFrac = 0.75*max(folsomODP,0.0)/cvp_ODP
      
    else:
      if total_ODP > 0.0:
        shastaFrac = max(shastaODP,0.0)/total_ODP
        folsomFrac = max(folsomODP,0.0)/total_ODP
        orovilleFrac = max(orovilleODP,0.0)/total_ODP
        yubaFrac = max(yubaODP,0.0)/total_ODP
      else:
        shastaFrac = 0.0
        folsomFrac = 0.0
        orovilleFrac = 0.0
        yubaFrac = 0.0
 
    return shastaFrac, folsomFrac, orovilleFrac, yubaFrac
	
  def time_pumping(self, t, cvpAS,swpAS):
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    wyt = self.forecastSCWYT
	
	####Old-Middle River Rule################################################################################################
	#The delta pumps make the old-middle river run backwards - there are restrictions (starting in 2008) about how big the negative flows can
	#become.  To model this, we use either a liner adjustment of flow at vernalis, or the observed flow (adding back in the pumping, to estimate
	#the 'natural' flow on the Old-Middle Rivers).  Once we have a model of the Old-Middle River, we assume that this flow is reduced (can become negative)
	# by 90% of the total pumping.  So, if the OMR limit is -5000CFS, & the natural flow is 1000CFS, pumping can be, maximum, 6000CFS/0.9
    if t < 4441:
      omrNat = self.vernalis_flow[t-1]*.462 + 120*cfs_tafd#vernalis flow linear adjustment
    else:
      omrNat = self.hist_OMR[t] + (self.hist_TRP_pump[t] + self.hist_HRO_pump[t])*.9
	  
    ## OMR max negative flows are only binding Jan-June, the rest of the year flow can be anything
    omr_condition = np.interp(dowy, self.omr_reqr['d'], self.omr_reqr['flow']) * cfs_tafd
	##The OMR constraints are governed by fish-kill conditions at the pumps.  negative 5000 CFS is normal, but can also be reduced
	##here we use the actual declaration - any forward-in-time projection runs will have to be either simulated probabilistically or just run under
	##the normal condition (5000 cfs)
	###Note: OMR Flow adjustments come from here: http://www.water.ca.gov/swp/operationscontrol/calfed/calfedwomt.cfm
    fish_trigger_adj = np.interp(t, self.omr_reqr['t'], self.omr_reqr['adjustment']) * cfs_tafd
    if t < 4441:
      omrRequirement = fish_trigger_adj
    else:
      omrRequirement = max(omr_condition, fish_trigger_adj)
    
    maxTotPumpInt = (omrNat - omrRequirement)/.9
    self.maxTotPump = max(maxTotPumpInt,0.0)    
    ##Apply OMR rules and other constraints on pumping
    cvp_max, swp_max = self.find_max_pumping(dowy, d, t, wyt)
    cvp_max, swp_max = self.meet_OMR_requirement(cvp_max, swp_max, t)
	##Check to see if there is enough projected storage (reservoir & snowpack) to meet these requirements (if not, cvp & swp_max are equal to zero)
    cvp_max, swp_max = self.find_release(t, dowy, cvp_max, swp_max, cvpAS, swpAS)
	
    return cvp_max, swp_max

  def calc_flow_bounds(self, t, cvp_max, swp_max, cvp_AS, swp_AS):
    ### stored and unstored flow agreement between SWP & CVP
	### project releases for pumping must consider the I/E 'tax'
	### releases already made (environmental flows, delta outflows) and downstream gains can be credited against this tax, but
	### we first have to determine how much of each is assigned to the delta outflow.
    
	##the first step is applying all 'unstored' flows to the delta outflows
	##note: unstored_flows should be strictly positive, but available_unstored can be negative (need outflow/environmental releases to meet delta outflow requirements)
    m = int(self.index.month[t])
    wyt = self.forecastSCWYT
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    outflow_rule = self.min_outflow[wyt][m-1] * cfs_tafd
	##Same salinity rule as in calc_flow_bounds
    if dowy > 180 and dowy < 318:
      if self.x2[t-1] > self.x2constraint[wyt][dowy]:
        x2outflow = 10**((self.x2constraint[wyt][dowy] - 10.16 - 0.945*self.x2[t])/(-1.487))
      else:
        x2outflow = 0.0
    else:
      x2outflow = 0.0
		
    salinity_rule = min(x2outflow*cfs_tafd,12000.0*cfs_tafd)

    min_rule = max(outflow_rule, salinity_rule)
    unstored_flows = self.total_inflow - self.cvp_stored_release - self.swp_stored_release
    available_unstored = unstored_flows + self.depletions[t] - min_rule	
	
    if available_unstored < 0:
      ##if unstored flows cannot meet delta outflow requirements, 75% of the deficit comes
	  ##from releases made by CVP (25% from SWP releases)
      cvp_frac = self.cvp_outflow/(self.cvp_outflow + self.swp_outflow)
      swp_frac = self.swp_outflow/(self.cvp_outflow + self.swp_outflow)
    else:
	  ##if unstored flows remain after meeting delta outflow requirements, the additional
	  ##flows are split 55/45 between CVP/SWP
      cvp_frac = 0.55
      swp_frac = 0.45
	
	##how much do we need to release for pumping diversions? (i.e. water balance)
	##we can use stored releases and any available unstored, as divided above (55/45)
	##if available unstored is less than zero, it comes out of stored release, as divided above (75/25)
    cvp_releases_for_pumping = cvp_max - self.cvp_stored_release - cvp_frac*available_unstored
    swp_releases_for_pumping = swp_max - self.swp_stored_release - swp_frac*available_unstored
	##how much (if any) do we need to release to meet the I/E tax? (i.e. environmental requirements)
	##unstored_flows is always positive, so we split them 55/45 to make the 'tax' on each project pumping
    cvp_releases_for_tax = cvp_max/self.export_ratio[wyt][m-1] - 0.55*unstored_flows - self.cvp_stored_release
    swp_releases_for_tax = swp_max/self.export_ratio[wyt][m-1] - 0.45*unstored_flows - self.swp_stored_release
    if cvp_AS > 0.0:
      cvp_extra_pump = max(self.export_ratio[wyt][m-1]*(0.75*(min_rule - self.depletions[t]))/(1-self.export_ratio[wyt][m-1]) - self.cvp_stored_release - 0.55*unstored_flows,0.0)
    else:
      cvp_extra_pump = 0.0

    if swp_AS > 0.0:
      swp_extra_pump = max(self.export_ratio[wyt][m-1]*(0.25*(min_rule - self.depletions[t]))/(1-self.export_ratio[wyt][m-1]) - self.swp_stored_release - 0.55*unstored_flows,0.0)
    else:
      swp_extra_pump = 0.0
    
	##need to release the larger of the two requirements
    self.sodd_cvp[t] = max(cvp_releases_for_pumping, cvp_releases_for_tax,cvp_extra_pump)
    self.sodd_swp[t] = max(swp_releases_for_pumping, swp_releases_for_tax,swp_extra_pump)
    ##unused flows from one project can be used to meet requirements for other project
    ##(i.e., if large stored releases (for environmental requirements) made from one project,
	##they may be greater than the total pumping/tax requirement, meaning the other project can
	##pump some of that water or use it to meet the 'tax'
    
    if self.sodd_cvp[t] < 0.0:
      extraCVP = self.sodd_cvp[t]
      self.sodd_cvp[t] = 0.0
    else:
      extraCVP = 0.0 
    if self.sodd_swp[t] < 0.0:
      extraSWP = self.sodd_cvp[t]
      self.sodd_swp[t] = 0.0
    else:
      extraSWP = 0.0
	##extraSWP/CVP variable will be negative  
    self.sodd_cvp[t] = max(self.sodd_cvp[t] + extraSWP, 0.0)
    self.sodd_swp[t] = max(self.sodd_swp[t] + extraCVP, 0.0)
 

		  
  def find_max_pumping(self, d, dowy, t, wyt):
    ##Find max pumping uses the delta pumping rules (from delta_properties) to find the maximum the pumps can
	##be run on a given day, from the D1641 and the BIOPS rules
	
    swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
    cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd
    san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.vernalis_flow[t], 0.0)
    
    if np.interp(t,self.san_joaquin_export_ratio['D1641_dates'],self.san_joaquin_export_ratio['D1641_on_off']) == 1:
      san_joaquin_ie_amt = np.interp(self.vernalis_flow[t]*tafd_cfs, self.san_joaquin_export_ratio['D1641_flow_target'],self.san_joaquin_export_ratio['D1641_export_limit']) * cfs_tafd
    else:
      san_joaquin_ie_amt = np.interp(self.vernalis_flow[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.vernalis_flow[t]
    
    san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    san_joaquin_ie = san_joaquin_ie_amt * san_joaquin_ie_used
    swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)

    return cvp_max, swp_max
   
  def meet_OMR_requirement(self, cvpm, swpm, t):
    ##Delta pumps reverse flow on the Old and Middle River, and they must be
	##operated so that the flow on those rivers does not exceed a certain negative (reverse) value
	
    cvp_m = cvpm
    swp_m = swpm
    if cvp_m + swp_m > self.maxTotPump:
      if cvp_m < self.maxTotPump*0.55:
        swp_m = self.maxTotPump - cvp_m
      elif swp_m < self.maxTotPump*0.45:
        cvp_m = self.maxTotPump - swp_m
      else:
        swp_m = self.maxTotPump*0.45
        cvp_m = self.maxTotPump*0.55
    return cvp_m, swp_m
  
  def find_release(self, t, dowy, cvp_max, swp_max, cvpAS, swpAS):
  ###This function looks at how much water is available in storage & snowpack to be exported,
  ##then determines when the best time to release that water is (water is saved for part of year
  ##when the inflow/export 'tax' is lowest)
    wyt = self.forecastSCWYT
    y = int(self.index.year[t])
    p1_stor_swp = (self.pump_max['swp']['intake_limit'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    p1_stor_cvp = (self.pump_max['cvp']['intake_limit'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    p2_stor_swp = (self.pump_max['swp']['intake_limit'][0] * cfs_tafd)/self.export_ratio[wyt][2]
    p2_stor_cvp = (self.pump_max['cvp']['intake_limit'][0] * cfs_tafd)/self.export_ratio[wyt][2]
    p3_stor_swp = (self.pump_max['swp']['intake_limit'][2] * cfs_tafd)/self.export_ratio[wyt][4]
    p3_stor_cvp = (self.pump_max['cvp']['intake_limit'][2] * cfs_tafd)/self.export_ratio[wyt][4]
  
    if dowy <=61:
      numdaysSave = 1
      numdaysSave2 = 0	  
    elif dowy <= 123:
      numdaysSave = 92
      numdaysSave2 = 0
    elif dowy <= 197:
      numdaysSave = 92
      numdaysSave2 = 196-dowy + 45
    elif dowy <= 227:
      numdaysSave = 92
      numdaysSave2 = 45  
    elif dowy <=273:
      numdaysSave = 92
      numdaysSave2 = 273-dowy
    else:
      numdaysSave = 1
      numdaysSave2 = 0
   
    if swpAS < (numdaysSave*p1_stor_swp + numdaysSave2*p2_stor_swp):
      swp_max = 0.0
    if cvpAS < (numdaysSave*p1_stor_cvp + numdaysSave2*p1_stor_cvp):
      cvp_max = 0.0
	  
    numdaysP1 = 92 - max(dowy-274, 0.0)##July-September
    numdaysP1a = max(123 - dowy, 0.0)##October - January
    numdaysP2 = max(59 - max(dowy-124, 0.0), 0.0)## February, March
    numdaysP3 = max(61 - max(dowy-180, 0.0), 0.0)## April, May
    numdaysP2a = max(30 - max(dowy - 241, 0.0), 0.0)## June

    totalP1 = min(swpAS, p1_stor_swp*numdaysP1)*self.export_ratio[wyt][8]
    totalP1a = min(swpAS - p1_stor_swp*numdaysP1, p1_stor_swp*numdaysP2)*self.export_ratio[wyt][8]
    totalP2 = min(swpAS - p1_stor_swp*(numdaysP1 + numdaysP1a), p2_stor_swp*numdaysP2)*self.export_ratio[wyt][2]
    totalP3 = min(swpAS - p1_stor_swp*(numdaysP1 + numdaysP1a) - p2_stor_swp*numdaysP2, p3_stor_swp*numdaysP3)*self.export_ratio[wyt][4]
    totalP2a = min(swpAS - p1_stor_swp*(numdaysP1 + numdaysP1a) - p2_stor_swp*numdaysP2 - p3_stor_swp*numdaysP3, p2_stor_swp*numdaysP2a)*self.export_ratio[wyt][2]
    self.forecastSWPPUMP = max(totalP1,0.0) + max(totalP1a,0.0) + max(totalP2,0.0) + max(totalP2a,0.0) + max(totalP3,0.0)

    totalP1 = min(cvpAS, p1_stor_cvp*numdaysP1)*self.export_ratio[wyt][8]
    totalP1a = min(cvpAS - p1_stor_cvp*numdaysP1, p1_stor_cvp*numdaysP2)*self.export_ratio[wyt][8]
    totalP2 = min(cvpAS - p1_stor_cvp*(numdaysP1 + numdaysP1a), p2_stor_cvp*numdaysP2)*self.export_ratio[wyt][2]
    totalP3 = min(cvpAS - p1_stor_cvp*(numdaysP1 + numdaysP1a) - p2_stor_cvp*numdaysP2, p3_stor_cvp*numdaysP3)*self.export_ratio[wyt][4]
    totalP2a = min(cvpAS - p1_stor_cvp*(numdaysP1 + numdaysP1a) - p2_stor_cvp*numdaysP2 - p3_stor_cvp*numdaysP3, p2_stor_cvp*numdaysP2a)*self.export_ratio[wyt][2]
    self.forecastCVPPUMP = max(totalP1,0.0) + max(totalP1a,0.0) + max(totalP2,0.0) + max(totalP2a,0.0) + max(totalP3,0.0)
        
    return cvp_max, swp_max


  def step(self, t, cvp_flows, swp_flows, unstored_flows):
    ##Takes releases (cvp_flows & swp_flows) and gains, and divides water between delta outflows, CVP exports, and SWP exports (and delta depletions)
	##Basically runs through the operations of calc_flow_bounds, only w/actual reservoir releases
    d = int(self.index.dayofyear[t])
    m = int(self.index.month[t])
    wyt = self.forecastSCWYT
    y = int(self.index.year[t])
    dowy = water_day(d)

	##Same outflow rule as in calc_flow_bounds
    outflow_rule = self.min_outflow[wyt][m-1] * cfs_tafd
	##Same salinity rule as in calc_flow_bounds
    if dowy > 180 and dowy < 318:
      if self.x2[t-1] > self.x2constraint[wyt][dowy]:
        x2outflow = 10**((self.x2constraint[wyt][dowy] - 10.16 - 0.945*self.x2[t])/(-1.487))
      else:
        x2outflow = 0.0
    else:
      x2outflow = 0.0
		
    salinity_rule = min(x2outflow*cfs_tafd,12000.0*cfs_tafd)

    min_rule = max(outflow_rule, salinity_rule)
    surplus = unstored_flows + self.depletions[t] - min_rule
    	
	#Same export ratio as in calc_weekly_storage_release
    export_ratio = self.export_ratio[wyt][m-1]
	#Same max pumping rules as in calc_weekly_storage release
    cvp_max, swp_max = self.find_max_pumping(d, dowy, t, wyt)

    if surplus < 0:
      ##if unstored flows cannot meet delta outflow requirements, 75% of the deficit comes
	  ##from releases made by CVP (25% from SWP releases)		
      cvp_frac = self.cvp_outflow/(self.cvp_outflow+self.swp_outflow)
      swp_frac = self.swp_outflow/(self.cvp_outflow+self.swp_outflow)
    else:
	  ##if unstored flows remain after meeting delta outflow requirements, the additional
	  ##flows are split 55/45 between CVP/SWP
      cvp_frac = 0.55
      swp_frac = 0.45

    self.TRP_pump[t] = max(min((cvp_flows + 0.55 * unstored_flows) * export_ratio, cvp_flows + cvp_frac * surplus),0.0)
    self.HRO_pump[t] = max(min((swp_flows + 0.45 * unstored_flows) * export_ratio, swp_flows + swp_frac * surplus),0.0)
		
    if self.TRP_pump[t] > cvp_max:
      self.TRP_pump[t] = cvp_max
      self.HRO_pump[t] = max(min(cvp_flows + swp_flows + surplus - self.TRP_pump[t], (cvp_flows + swp_flows + unstored_flows)*export_ratio - self.TRP_pump[t],swp_max),0.0)
  
    if self.HRO_pump[t] > swp_max:
      self.HRO_pump[t] = swp_max
      self.TRP_pump[t] = max(min(cvp_flows + swp_flows + surplus - self.HRO_pump[t], (cvp_flows + swp_flows + unstored_flows)*export_ratio - self.HRO_pump[t],cvp_max),0.0)
    
    if self.TRP_pump[t] < 0.0:
      self.HRO_pump[t] = max(self.HRO_pump[t] + self.TRP_pump[t],0.0)
      self.TRP_pump[t] = 0.0
    elif self.HRO_pump[t] < 0.0:
      self.TRP_pump[t] = max(self.TRP_pump[t] + self.HRO_pump[t],0.0)
      self.HRO_pump[t] = 0.0


	##Same as in calc_weekly_storage_release
    self.TRP_pump[t], self.HRO_pump[t] = self.meet_OMR_requirement(self.TRP_pump[t], self.HRO_pump[t], t)

    self.outflow[t] = cvp_flows + swp_flows + unstored_flows - self.TRP_pump[t] - self.HRO_pump[t] + self.depletions[t]

	##Calculate X2 values, for salinity rules - note, if outflow is negative (may happen under extreme conditions in which not enough storage to meet negative gains)
	##X2 is calculated w/an outflow of 50 cfs b/c log calcs require positive flow
    if self.outflow[t] > 0.0:
      self.x2[t+1] = 10.16 + 0.945*self.x2[t] - 1.487*np.log10(self.outflow[t]*tafd_cfs)
    else:
      self.x2[t+1] = 10.16 + 0.945*self.x2[t] - 1.487*np.log10(50.0)
	  
    self.OMR[t] = self.hist_OMR[t] + self.hist_TRP_pump[t] + self.hist_HRO_pump[t] - self.TRP_pump[t] - self.HRO_pump[t]
    self.calc_project_allocation(t)
  	
  def calc_project_allocation(self,t):
    y = int(self.index.year[t])
    m = int(self.index.month[t])
    startYear = int(self.index.year[0])
    if m < 10:
      wateryear = y -(startYear + 1)
    else:
      wateryear = y - startYear       
	
    self.swp_allocation[t] = self.forecastSWPPUMP + self.annual_HRO_pump[wateryear-1]
    self.cvp_allocation[t] = self.forecastCVPPUMP + self.annual_TRP_pump[wateryear-1]
    self.annual_HRO_pump[wateryear-1] += self.HRO_pump[t]
    self.annual_TRP_pump[wateryear-1] += self.TRP_pump[t]
	  
  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['TRP_pump','HRO_pump', 'total_outflow','SWP_allocation', 'CVP_allocation', 'X2']
    things = [self.TRP_pump, self.HRO_pump, self.outflow, self.swp_allocation, self.cvp_allocation, self.x2]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df