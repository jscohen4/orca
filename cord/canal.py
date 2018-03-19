from __future__ import division
import numpy as np 
import pandas as pd
import json
from .util import *

class Canal():

  def __init__(self, df, key):
    self.T = len(df)
    self.index = df.index
    self.key = key

    for k,v in json.load(open('cord/data/%s_properties.json' % key)).items():
        setattr(self,k,v)
	
    self.cvp_aval_stor = 0.5
    self.swp_aval_stor = 0.5
    self.cvp_delta_outflow_pct = 0.0
    self.swp_delta_outflow_pct = 0.0
    self.swp_allocation = np.zeros(self.T)
    self.cvp_allocation = np.zeros(self.T)
    self.fk_class1 = np.zeros(self.T)
    self.fk_class2 = np.zeros(self.T)
    self.annual_fk1 = np.zeros((self.index.year[self.T-1]-self.index.year[0]))
    self.annual_fk2 = np.zeros((self.index.year[self.T-1]-self.index.year[0]))
    self.annual_deliveries = np.zeros((self.index.year[self.T-1]-self.index.year[0]))
    self.total_diversions = np.zeros(self.T)
    self.annual_swp_deliveries = np.zeros((self.index.year[self.T-1]-self.index.year[0]))
    self.annual_cvp_deliveries = np.zeros((self.index.year[self.T-1]-self.index.year[0]))
    self.total_swp_diversions = np.zeros(self.T)
    self.total_cvp_diversions = np.zeros(self.T)
    self.swp_spill = np.zeros(self.T)
    self.cvp_spill = np.zeros(self.T)
    self.forecastSWPPUMP = 0.0
    self.forecastCVPPUMP = 0.0
    self.storage215 = 0.0
    self.cumulativeInflow = 0.0
	

  def calc_friant_kern(self, t, millerton_available, current_storage, summer_flows):
	
	tot_deliverySchedule = [65.175, 107.8, 103.125, 204.875, 341.825, 499.1250, 544.225, 402.05, 226.6, 165.275, 59.95, 30.25]
    fk_deliverySchedule = [59.84, 100.32, 91.08, 166.1, 268.4, 364.32, 400.4, 318.78, 198.44, 147.84, 55.66, 28.82]
    totalC1A = 800#total 'Class 1' Irrigation water allocations from Millerton Lake (thousand acre feet)
    totalC2A = 1400#total 'Class 2' irrigation water allocations from Millerton Lake (thousand acre feet)
    
    todayFK1 = self.fk_class1[t]*(8/22)*fk_deliverySchedule[m-1]/days_in_month[m-1]
    todayFK2 = self.fk_class2[t]*(14/22)*fk_deliverySchedule[m-1]/days_in_month[m-1]
    todayFKMAD = (self.fk_class1[t]*(8/22) + self.fk_class2[t]*(14/22))*tot_deliverySchedule[m-1]/days_in_month[m-1]
    if todayFK1 > millerton_available*0.8:
      todayFK1 = millerton_available*0.8
      todayFK2 = 0.0
    elif (todayFK1+todayFK2) > millerton_available*0.8:
      todayFK2 = millerton_available*0.8 - todayFK1
	  
 
    if todayFKMAD > millerton_available:
      todayFKMAD = millerton_available
    if m == 4 and da == 1:
      total_room = millertonCapacity - current_storage + (tot_deliverySchedule[4]+tot_deliverySchedule[5])*(self.fk_class1[t]*(8/22)+self.fk_class2[t]*(14/22))
      total_flows = summer_flows/2.0
      if total_flows > total_room:
        self.storage215 = total_flows - total_room
      else:
        self.storage215 = 0.0
	      
    fkcap = 5000*cfs_tafd
    spill215 = max(min(fkcap-(todayFK1+todayFK2),self.storage215),0.0)
    self.storage215 -= spill215
    self.annual_fk1[wateryear] += (todayFK1)
    self.annual_fk2[wateryear] += (todayFK2)
    self.annual_deliveries[wateryear] += (todayFKMAD+spill215)
    self.total_diversions[t] = (todayFKMAD+spill215)
	
  def calc_ca_aqueduct(self, t):
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    da = int(self.index.day[t])
    startYear = int(self.index.year[0])
    if m < 10:
      wateryear = y -(startYear + 1)
    else:
      wateryear = y - startYear
    if m == 10 and da == 1:
      self.max_swp_forecast = 4100
      self.max_cvp_forecast = 2600

    self.max_swp_forecast -= np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
    self.max_cvp_forecast -= np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['intake_limit']) * cfs_tafd
	
    self.swp_allocation[t] = self.SL_S_stor[t]-42.0 + min(self.forecastSWPPUMP,self.max_swp_forecast) + self.annual_swp_deliveries[wateryear-1]
    self.cvp_allocation[t] = self.SL_F_stor[t]-42.0 + min(self.forecastCVPPUMP,self.max_cvp_forecast) + self.annual_cvp_deliveries[wateryear-1]
	
    swp_deliverySchedule = [.0418, .0575, .0654, .0854, .1131, .1564, .1422, .0971, .0812, .0579, .0571, .0450]
    cvp_deliverySchedule = [.0560, .0626, .0684, .1084, .1555, .1669, .1273, .0822, .0678, .0396, .0313, .0340]

    self.annual_swp_deliveries[wateryear-1] += (self.swp_allocation[t])*swp_deliverySchedule[m-1]/days_in_month[m-1]
    self.annual_cvp_deliveries[wateryear-1] += (self.cvp_allocation[t])*cvp_deliverySchedule[m-1]/days_in_month[m-1]
    self.total_swp_diversions[t] = (self.swp_allocation[t])*swp_deliverySchedule[m-1]/days_in_month[m-1]
    self.total_cvp_diversions[t] = (self.cvp_allocation[t])*cvp_deliverySchedule[m-1]/days_in_month[m-1]
    self.SL_S_stor[t+1] -= (self.swp_allocation[t])*swp_deliverySchedule[m-1]/days_in_month[m-1]
    self.SL_F_stor[t+1] -= (self.cvp_allocation[t])*cvp_deliverySchedule[m-1]/days_in_month[m-1]
    if self.SL_S_stor[t+1] < 0.0:
      self.SL_S_stor[t+1] = 0.0
    if self.SL_F_stor[t+1] < 0.0:
      self.SL_F_stor[t+1] = 0.0

  def calc_vernalis_rule(self,t, tuolumne_flow, merced_flow, stanislaus_flow, NMI):
    d = int(self.index.dayofyear[t])
    dowy = water_day(d)
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    wyt = self.forecastSJWYT[t]
    
    if y > 2009:
      biops_min = np.interp(dowy,self.san_joaquin_min['biops_d'],self.san_joaquin_min['biops_on_off'])*np.interp(NMI,self.san_joaquin_min['biops_NMI'],self.san_joaquin_min['biops_flow']) * cfs_tafd
    else:
      biops_min = 0.0
    
    zone = int(np.interp(dowy, self.san_joaquin_min_flow['d'], self.san_joaquin_min_flow['zone']))
    d_1641_min = self.san_joaquin_min_flow[wyt][zone-1] * cfs_tafd
    self.vernalis_gains = self.gains_sj[t] + tuolumne_flow + merced_flow + stanislaus_flow
    SJGRA_flow = max(d_1641_min - self.vernalis_gains, 0.0)
    if zone == 4:
      merced_contr = SJGRA_flow*0.5
      tuolumne_contr = SJGRA_flow*0.25
      stanislaus_contr = SJGRA_flow*.025
    else:
      stanislaus_contr = SJGRA_flow
      merced_contr = 0.0
      tuolumne_contr = 0.0


    if biops_min > d_1641_min:
      stanislaus_contr += max(biops_min-self.vernalis_gains,0.0) - max(d_1641_min-self.vernalis_gains,0.0)

    return merced_contr, tuolumne_contr, stanislaus_contr
  
  def calc_flow_bounds(self, t, folsom_nodd, shasta_nodd, oroville_nodd,orovilleAS,folsomAS,shastaAS,yubaAS):
    d = int(self.index.dayofyear[t])
    m = int(self.index.month[t])
    y = int(self.index.year[t])
    
    wyt = self.forecastSCWYT[t]
	#northern project deliveries don't count as in-basin uses for stored/unstored flow
    folsom_n_proj_consump = np.interp(d, first_of_month, folsom_nodd)
    shasta_n_proj_consump = np.interp(d, first_of_month, shasta_nodd)
    oroville_n_proj_consump = np.interp(d, first_of_month, oroville_nodd)
    self.gains[t] = self.vernalis_flow[t] + self.eastside_streams[t] + self.gains_sac[t] + folsom_n_proj_consump + shasta_n_proj_consump + oroville_n_proj_consump
    outflow_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    monthindex = int(m - int(self.index.month[0]) + (y - int(self.index.year[0]))*12) 
    if self.x2[t] > self.salinity_target[wyt][m-1]:	
      salinity_rule = 11400*cfs_tafd
    else:
      salinity_rule = 0.0
	  
    min_rule = max(outflow_rule, salinity_rule)	 
    self.calc_weekly_storage_release(t,min_rule,orovilleAS,folsomAS,shastaAS,yubaAS)
  
  def calc_weekly_storage_release(self,t,min_rule,orovilleAS,folsomAS,shastaAS, yubaAS):
    ##this function takes the total storage available for release at each reservoir and distributes that throughout the year based on the Export/Inflow
    ##ratio 'tax' that exists in that week
    d = self.index.dayofyear[t]
    dowy = water_day(d)
    wyt = self.forecastSCWYT[t]
    m = int(self.index.month[t])
    if t < 4441:
      omrNat = self.vernalis_flow[t-1]*.462 + 120*cfs_tafd
    else:
      omrNat = self.hist_OMR[t] + (self.hist_TRP_pump[t] + self.hist_HRO_pump[t])*.9
	###Note: OMR Flow adjustments come from here: http://www.water.ca.gov/swp/operationscontrol/calfed/calfedwomt.cfm
    omr_condition = np.interp(dowy, self.omr_reqr['d'], self.omr_reqr['flow']) * cfs_tafd
    fish_trigger_adj = np.interp(t, self.omr_reqr['t'], self.omr_reqr['adjustment']) * cfs_tafd
    omrRequirement = max(omr_condition, fish_trigger_adj)
    maxTotPumpInt = (omrNat - omrRequirement)/.9
    self.maxTotPump = max(maxTotPumpInt,0.0)
    	
    cvp_max, swp_max = self.find_max_pumping(dowy, d, t, wyt)

    cvp_max, swp_max = self.meet_OMR_requirement(cvp_max, swp_max, t)

    cvp_max, swp_max = self.check_san_luis(cvp_max, swp_max, t, 0)  
	
    surplus = self.gains[t] + self.depletions[t] - min_rule
	
    cvp_max, swp_max = self.find_release(t, dowy, cvp_max, swp_max, orovilleAS, shastaAS, folsomAS, yubaAS)

    if surplus < 0: # additional flow needed to meet delta reqmnt
      self.dmin[t] = surplus * - 1.0
      self.sodd_cvp[t] = max(cvp_max/self.export_ratio[wyt][m-1] - 0.75 * self.gains[t]  - self.cvp_delta_outflow_pct * self.dmin[t], cvp_max, 0.0)
      self.sodd_swp[t] = max(swp_max/self.export_ratio[wyt][m-1] - 0.25 * self.gains[t] - self.swp_delta_outflow_pct * self.dmin[t], swp_max, 0.0)
    else: # extra unstored water available
      self.sodd_cvp[t] = max(cvp_max/self.export_ratio[wyt][m-1] - 0.45 * self.gains[t],0)
      self.sodd_swp[t] = max(swp_max/self.export_ratio[wyt][m-1] - 0.55 * self.gains[t],0)
    
    surplus_sac = self.sodd_cvp[t] + self.sodd_swp[t] + self.gains_sac[t] + self.dmin[t] - self.sac_min_flow[wyt][m-1]*cfs_tafd
 	
    if surplus_sac < 0.0:
      self.dmin[t] += surplus_sac * -1.0

	    	  		
  def check_san_luis(self, tp, hp, t, update):
    sls_stor = self.SL_S_stor[t] + hp
    slf_stor = self.SL_F_stor[t] + tp
    if sls_stor > self.sl_cap/2.0:
      if slf_stor < self.sl_cap/2.0:
        if sls_stor > (self.sl_cap - slf_stor):
          self.swp_spill[t] = sls_stor - (self.sl_cap - slf_stor)
          self.cvp_spill[t] = 0.0
          sls_stor = self.sl_cap - slf_stor
        else:
          self.swp_spill[t] = 0.0
          self.cvp_spill[t] = 0.0    
      else:
        self.swp_spill[t] = (sls_stor - self.sl_cap/2.0)
        self.cvp_spill[t] = (slf_stor - self.sl_cap/2.0)
        sls_stor = self.sl_cap/2.0
        slf_stor = self.sl_cap/2.0
    elif slf_stor > self.sl_cap - sls_stor:
      self.cvp_spill[t] = slf_stor - (self.sl_cap - sls_stor)
      self.swp_spill[t] = 0.0
      slf_stor = self.sl_cap - sls_stor
    else:
      self.swp_spill[t] = 0.0
      self.cvp_spill[t] = 0.0

    
    if update == 1:
      self.SL_S_stor[t+1] = sls_stor
      self.SL_F_stor[t+1] = slf_stor
	  
    tp = max(tp,0.0)
    hp = max(hp,0.0)
    return tp, hp	
  
  
  def meet_OMR_requirement(self, cvpm, swpm, t):
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
  
  def find_release(self, t, dowy, cvp_max, swp_max, orovilleAS, shastaAS, folsomAS, yubaAS):
    wyt = self.forecastSCWYT[t]
    y = int(self.index.year[t])
    swp_jas_stor = (self.pump_max['swp']['pmax'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    cvp_jas_stor = (self.pump_max['cvp']['pmax'][5] * cfs_tafd)/self.export_ratio[wyt][8]
    swp_daj_stor = (self.pump_max['swp']['pmax'][3] * cfs_tafd)/self.export_ratio[wyt][3]
    cvp_daj_stor = (self.pump_max['cvp']['pmax'][3] * cfs_tafd)/self.export_ratio[wyt][3]
  
    shastaAS = max(shastaAS,0.0)
    orovilleAS = max(orovilleAS, 0.0)
    folsomAS = max(folsomAS,0.0)
    yubaAS = max(yubaAS,0.0)
	
    self.cvp_aval_stor = shastaAS + folsomAS
    self.swp_aval_stor = orovilleAS + yubaAS

    if dowy <= 123:
      numdaysSave = 92
      numdaysSave2 = 0
    elif dowy <=274:
      numdaysSave = 92
      numdaysSave2 = 274-dowy
    else:
      numdaysSave = 1
      numdaysSave2 = 0
   
    if self.swp_aval_stor < (numdaysSave*swp_jas_stor + numdaysSave2*swp_daj_stor):
      swp_max = 0.0
    if self.cvp_aval_stor < (numdaysSave*cvp_jas_stor + numdaysSave2*cvp_daj_stor):
      cvp_max = 0.0
    
    numdaysP1 = 92 - max(dowy-274, 0)
    numdaysP2 = max(123 - dowy,0.0)
    numdaysP3 = max(150 - max(dowy-124, 0.0),0.0)
    self.forecastSWPPUMP = max(min(self.swp_aval_stor,swp_jas_stor*numdaysP1),0.0)*self.export_ratio[wyt][8] + max(min(self.swp_aval_stor-swp_jas_stor*numdaysP1,swp_jas_stor*numdaysP2),0.0)*self.export_ratio[wyt][8] + max(min(self.swp_aval_stor-swp_jas_stor*(numdaysP1+numdaysP2),cvp_daj_stor*numdaysP3),0.0)*self.export_ratio[wyt][3]
    self.forecastCVPPUMP = max(min(self.cvp_aval_stor,cvp_jas_stor*numdaysP1),0.0)*self.export_ratio[wyt][8] + max(min(self.cvp_aval_stor-cvp_jas_stor*numdaysP1,cvp_jas_stor*numdaysP2),0.0)*self.export_ratio[wyt][8] + max(min(self.cvp_aval_stor-cvp_jas_stor*(numdaysP1+numdaysP2),cvp_daj_stor*numdaysP3),0.0)*self.export_ratio[wyt][3]
    
    return cvp_max, swp_max


  def step(self, t, cvp_flows, swp_flows):
    d = int(self.index.dayofyear[t])
    m = int(self.index.month[t])
    wyt = self.forecastSCWYT[t]
    y = int(self.index.year[t])
    dowy = water_day(d)

    self.inflow[t] = self.gains[t] + cvp_flows + swp_flows
    outflow_rule = self.min_outflow[wyt][m-1] * cfs_tafd
    monthindex = int(m - int(self.index.month[0]) + (y - int(self.index.year[0]))*12)
    if self.x2[t] > self.salinity_target[wyt][m-1]:	
      salinity_rule = 11400.0*cfs_tafd
    else:
      salinity_rule = 0.0
	  
    min_rule = max(outflow_rule, salinity_rule)	 
    export_ratio = self.export_ratio[wyt][m-1]
    cvp_max, swp_max = self.find_max_pumping(d, dowy, t, wyt)

    surplus = self.gains[t] + self.depletions[t] - min_rule
    if surplus > 0:
      self.TRP_pump[t] = max(min((cvp_flows + 0.55 * self.gains[t]) * export_ratio, cvp_flows + 0.55 * surplus, cvp_max),0)
      self.HRO_pump[t] = max(min((swp_flows + 0.45 * self.gains[t]) * export_ratio, swp_flows + 0.45 * surplus, swp_max),0)
    else:
      self.TRP_pump[t] = min((cvp_flows + 0.75 * self.gains[t]) * export_ratio, cvp_flows + 0.75 * surplus, cvp_max)
      self.HRO_pump[t] = min((swp_flows + 0.25 * self.gains[t]) * export_ratio, swp_flows + 0.25 * surplus, swp_max)
      if self.TRP_pump[t] < 0.0:
        self.HRO_pump[t] += self.TRP_pump[t]
        self.TRP_pump[t] = 0.0
      elif self.HRO_pump[t] < 0.0:
        self.TRP_pump[t] += self.HRO_pump[t]
        self.HRO_pump[t] = 0.0
	
    self.TRP_pump[t], self.HRO_pump[t] = self.meet_OMR_requirement(self.TRP_pump[t], self.HRO_pump[t], t)
    self.TRP_pump[t], self.HRO_pump[t] = self.check_san_luis(self.TRP_pump[t], self.HRO_pump[t], t, 1)
	
    self.outflow[t] = self.inflow[t] - self.TRP_pump[t] - self.HRO_pump[t] + self.depletions[t]
    
    if self.outflow[t] < 0.0:
      print(t, end = " ")
      print("%.2f" % self.gains[t], end = " ")
      print("%.2f" % surplus, end = " ")
      print("%.2f" % self.depletions[t], end = " ")
      print("%.2f" % cvp_flows, end = " ")
      print("%.2f" % swp_flows, end = " ")
      print("%.2f" % min_rule, end = " ")
      print("%.2f" % export_ratio, end = " ")
      print("%.2f" % self.inflow[t], end = " ")
      print("%.2f" % self.TRP_pump[t], end = " ")
      print("%.2f" % self.HRO_pump[t])

    if self.outflow[t] > 0.0:
      self.x2[t+1] = 10.16 + 0.945*self.x2[t] - 1.487*np.log10(self.outflow[t]*tafd_cfs)
    else:
      self.x2[t+1] = 10.16 + 0.945*self.x2[t] - 1.487*np.log10(50.0)
	  
    self.OMR[t] = self.hist_OMR[t] + self.hist_TRP_pump[t] + self.hist_HRO_pump[t] - self.TRP_pump[t] - self.HRO_pump[t]
    self.cumulativeInflow += self.inflow[t]
    #if t > 360:
      #self.cumulativeInflow -= self.inflow[t-360]
      #if self.cumulativeInflow < 10500.0:
        #totalPumping = self.TRP_pump[t] + self.HRO_pump[t]
        #magnitudeDecrease = max(10.0* (10500.0 - self.cumulativeInflow)/1500.0,10.0)
        #totalReduction = min(magnitudeDecrease, max(totalPumping - 1.5,0.0))
        #hroAdjustment = totalReduction*self.HRO_pump[t]/totalPumping
        #trpAdjustment = totalReduction*self.TRP_pump[t]/totalPumping
        #self.HRO_pump[t] = max(self.HRO_pump[t] - hroAdjustment,0.0)
        #self.TRP_pump[t] = max(self.TRP_pump[t] - trpAdjustment, 0.0)		
  	
	
  def find_max_pumping(self, d, dowy, t, wyt):
    swp_intake_max = np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['intake_limit']) * cfs_tafd
    cvp_intake_max = np.interp(d, self.pump_max['cvp']['d'],self.pump_max['cvp']['intake_limit']) * cfs_tafd
    san_joaquin_adj = np.interp(dowy, self.san_joaquin_add['d'], self.san_joaquin_add['mult']) * max(self.vernalis_flow[t], 0.0)
    
    if np.interp(t,self.san_joaquin_export_ratio['D1641_dates'],self.san_joaquin_export_ratio['D1641_on_off']) == 1:
      san_joaquin_ie_amt = np.interp(self.vernalis_flow[t]*tafd_cfs, self.san_joaquin_export_ratio['D1641_flow_target'],self.san_joaquin_export_ratio['D1641_export_limit']) * cfs_tafd
    else:
      san_joaquin_ie_amt = np.interp(self.vernalis_flow[t]*tafd_cfs, self.san_joaquin_export_ratio['flow'], self.san_joaquin_export_ratio['ratio']) * self.vernalis_flow[t]
    
    san_joaquin_ie_used = np.interp(dowy, self.san_joaquin_export_ratio['d'], self.san_joaquin_export_ratio['on_off'])
    san_joaquin_ie = san_joaquin_ie_amt * san_joaquin_ie_used
    swp_max = min(max(swp_intake_max + san_joaquin_adj, san_joaquin_ie * 0.55), np.interp(d, self.pump_max['swp']['d'], self.pump_max['swp']['pmax']) * cfs_tafd)
    cvp_max = min(max(cvp_intake_max, san_joaquin_ie * 0.45), np.interp(d, self.pump_max['cvp']['d'], self.pump_max['cvp']['pmax']) * cfs_tafd)

    return cvp_max, swp_max
		
  #def assign_flow(self, deltadata):
    #self.hist_inflows = deltadata['TOT'].values * cfs_tafd
    #self.sanjoaquin = deltadata['SJR'].values * cfs_tafd
    #self.cccpumping = deltadata['CCC'].values * cfs_tafd
    #self.barkerpump = deltadata['NBAQ'].values * cfs_tafd
    #self.depletions = deltadata['GCD'].values * cfs_tafd
    #self.misc =  deltadata['MISDV'].values * cfs_tafd
	
  #def find_gains(self, timestep, folsomreleases, shastareleases, orovillereleases):2
      #self.gains_d[timestep-1] = self.hist_inflows[timestep-1] - folsomreleases - shastareleases - orovillereleases
      
  def results_as_df(self, index):
    df = pd.DataFrame()
    names = ['TRP_pump','HRO_pump','CVP_allocation','SWP_allocation','Friant1_allocation','Friant2_allocation','swp_delivery','cvp_delivery','fk_delivery','tot_inflow','tot_outflow','vernalis_inflow','SCI','SJI','SanLuisS','SanLuisF','SWPForecast','CVPForecast','SWPStorage','CVPStorage', 'X2']
    things = [self.TRP_pump, self.HRO_pump, self.cvp_allocation, self.swp_allocation, self.fk_class1, self.fk_class2,self.total_swp_diversions,self.total_cvp_diversions,self.total_diversions,self.inflow,self.outflow, self.vernalis_flow, self.forecastSRI, self.forecastSJI,self.SL_S_stor,self.SL_F_stor,self.forecastSWPPUMP,self.forecastCVPPUMP, self.swp_aval_stor,self.cvp_aval_stor, self.x2]
    for n,t in zip(names,things):
      df['%s_%s' % (self.key,n)] = pd.Series(t, index=index)
    return df