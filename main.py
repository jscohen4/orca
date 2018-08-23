import numpy as np
np.warnings.filterwarnings('ignore') #to not display numpy warnings... be careful
import pandas as pd
import matplotlib.pyplot as plt
from orca import *
from orca.data import *
from subprocess import call
from datetime import datetime
now = datetime.now().strftime('Last modified %Y-%m-%d %H:%M:%S')
#Each of these booleans determines the actions that will be run by the model 

projection = False #True if running a single climate projection
calc_R2s = True #True if calculating R2s (only relevant for historical scenario)
plot = True #True if plotting outputs, need calc_R2s to also be true if plotting historical results!!!!
change_inflow_exeedance = False

#######Define a few parameters
SHA_shift = 0
ORO_shift = 0
FOL_shift = 0
index_exceedance_sac = 8
window_type = 'rolling'
window_length = 50
SHA_exceedance = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
ORO_exceedance = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
FOL_exceedance = {"W": 10, "AN": 10, "BN": 5, "D": 2, "C": 1}

#########these are the parameters originally used for historical runs
# SHA_shift = 0
# ORO_shift = 0
# FOL_shift = 0
# index_exceedance_sac = 8
# window_type = 'historical'
# window_length = 50
SHA_exceedance_hist = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
ORO_exceedance_hist = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
FOL_exceedance_hist = {"W": 10, "AN": 10, "BN": 5, "D": 2, "C": 1}


process_hist_data = False#True if changing any historical data inputs, or downloading updated data from cdec
###Only relevant if processing historical data
cdec = False # True if downloading up-to-date cdec data
hist_indices = True #True if running calc_indices scriptwater_day
hist_forcast = True #True if running updated forecast

sc = 'access1-0_rcp45_r1i1p1' #cmip5 climate scenario to use, if projection = True
process_climate_data = False #only mark True if running climate projection and/or processing projection input data
####### only relevant if processing projection data
climate_indices = True
climate_forecasts = True
#Nothing below here should be changed!
###############################################
###############################################
###############################################
if process_hist_data or not projection: 
  text_file = open("orca/data/historical_runs_data/datetime.txt", "w")
  text_file.write("%s" %now)
  text_file.close()

if process_climate_data or projection:
  call(['mkdir', 'orca/data/individual_projection_runs/%s'%sc])
  text_file = open("orca/data/individual_projection_runs/%s/datetime.txt"%sc, "w")
  text_file.write("%s" %now)
  text_file.close()
####################### below this line data processing and model runs are executed

if process_hist_data: 
  from orca.data import *
  if cdec:
    cdec_df = cdec_scraper.scrape_cdec()
    cdec_df.to_csv('orca/data/historical_runs_data/cdec-data.csv')
  if hist_indices:
    if not cdec: 
      cdec_df = pd.read_csv('orca/data/historical_runs_data/cdec-data.csv', index_col=0, parse_dates=True)
    ind_df, gains_df, OMR_df= process(cdec_df,'orca/data/json_files/evap_regression.json','orca/data/json_files/gains_regression.json','orca/data/json_files/inf_regression.json')  
    ind_df.to_csv('orca/data/historical_runs_data/orca-data-processed.csv')
    gains_df.to_csv('orca/data/historical_runs_data/gains_loops.csv')
    OMR_df.to_csv('orca/data/historical_runs_data/OMR_loops.csv')
  if hist_forcast:
    if not hist_indices:
      ind_df = pd.read_csv('orca/data/historical_runs_data/orca-data-processed.csv', index_col=0, parse_dates=True)
    forc_df, stats_df, WYI_stats= forecast(ind_df,8)
    forc_df.to_csv('orca/data/historical_runs_data/orca-data-forecasted.csv')
    stats_df.to_csv('orca/data/forecast_regressions/carryover_regression_statistics.csv')
    WYI_stats.to_csv('orca/data/forecast_regressions/WYI_forcasting_regression_stats.csv')

if change_inflow_exeedance:
  write_json.modify('orca/data/json_files/SHA_properties.json', 'exceedance', SHA_exceedance)
  write_json.modify('orca/data/json_files/ORO_properties.json', 'exceedance', ORO_exceedance)
  write_json.modify('orca/data/json_files/FOL_properties.json', 'exceedance', FOL_exceedance)

if not projection:
  model = Model('orca/data/historical_runs_data/orca-data-forecasted.csv', 'orca/data/historical_runs_data/orca-data-forecasted.csv',SHA_shift, ORO_shift, FOL_shift,sd='10-01-1999',projection = False, sim_gains = False) #beacuse of rolling calc in gains, we start on 10th day of
  results = model.simulate() # takes a while... save results
  results.to_csv('orca/data/historical_runs_data/results.csv')
  if calc_R2s:
    results['Combined_pump'] = results['DEL_HRO_pump'] + results['DEL_TRP_pump']
    sim = [results['DEL_HRO_pump'] / cfs_tafd,
         results['DEL_TRP_pump'] / cfs_tafd, 
         # (results['DEL_HRO_pump'] + results['DEL_TRP_pump']) / cfs_tafd,
         results['Combined_pump']/ cfs_tafd,
         results['SHA_storage'], 
         results['SHA_out'] / cfs_tafd,
         results['FOL_storage'],
         results['FOL_out'] / cfs_tafd,
         results['ORO_storage'],
         results['ORO_out'] / cfs_tafd,
         results['DEL_in'] / cfs_tafd,
         results['DEL_out'] / cfs_tafd, 
         results['DEL_X2']]
    obs = [model.df['HRO_pump'],
           model.df['TRP_pump'],
           (model.df['HRO_pump'] + model.df['TRP_pump']),
           model.df['SHA_storage'],
           model.df['SHA_out'],
           model.df['FOL_storage'],
           model.df['FOL_out'],
           model.df['ORO_storage'],
           model.df['ORO_out'],
           model.df['DeltaIn'],
           model.df['DeltaOut'], 
           model.df['X2']]
    plotter.Rsquares(sim,obs,'orca/data/historical_runs_data/Rsquares.txt')
    if plot:
      calibr_pts = ['HRO_pump','TRP_pump','Combined_pump','SHA_storage','SHA_out','FOL_storage','FOL_out','ORO_storage','ORO_out','DeltaIn','DeltaOut','X2']
      text_file = open("orca/figs/historical/datetime.txt", "w")
      text_file.write("%s" %now)
      text_file.close()
      for f in ['D','W','M','AS-OCT']:
        for s,o,c in zip(sim,obs,calibr_pts):
          plotter.compare(s, o, freq=f)
          plt.savefig('orca/figs/historical/%s_%s.pdf' % (f,c), dpi=150)
          plt.close()  

if process_climate_data:
  from orca.data import *
  call(['cp','orca/data/input_climate_files/%s_input_data.csv'%sc,'orca/data/individual_projection_runs/%s/%s_input_data.csv'%(sc,sc)]) 
  if climate_indices:
    input_df = pd.read_csv('orca/data/individual_projection_runs/%s/%s_input_data.csv'%(sc,sc), index_col = 0, parse_dates = True)
    gains_loop_df = pd.read_csv('orca/data/historical_runs_data/gains_loops.csv', index_col = 0, parse_dates = True)
    OMR_loop_df = pd.read_csv('orca/data/historical_runs_data/OMR_loops.csv', index_col = 0, parse_dates = True)
    proj_ind_df = process_projection(input_df,gains_loop_df,OMR_loop_df,'orca/data/json_files/gains_regression.json','orca/data/json_files/inf_regression.json',window = window_type)  
    proj_ind_df.to_csv('orca/data/individual_projection_runs/%s/orca-data-processed-%s.csv'%(sc,sc))
  if climate_forecasts:
    if not climate_indices:
      proj_ind_df = pd.read_csv('orca/data/individual_projection_runs/%s/orca-data-processed-%s.csv'%(sc,sc), index_col = 0, parse_dates = True)
    WYI_stats_file = pd.read_csv('orca/data/forecast_regressions/WYI_forcasting_regression_stats.csv', index_col = 0, parse_dates = True)
    carryover_stats_file = pd.read_csv('orca/data/forecast_regressions/carryover_regression_statistics.csv', index_col = 0, parse_dates = True)
    forc_df= projection_forecast(proj_ind_df,WYI_stats_file,carryover_stats_file,window_type,window_length, index_exceedance_sac)
    forc_df.to_csv('orca/data/individual_projection_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc))

if projection:
  model = Model('orca/data/individual_projection_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc), 'orca/data/historical_runs_data/results.csv',SHA_shift, ORO_shift, FOL_shift,sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
  results = model.simulate() # takes a while... save results
  results.to_csv('orca/data/individual_projection_runs/%s/%s-results.csv'%(sc,sc))
# calibration points (lists of pandas series)
# results = pd.read_csv('orca/data/results.csv', index_col=0, parse_dates=True)
  results['Combined_pump'] = results['DEL_HRO_pump'] + results['DEL_TRP_pump']
  sim = [results['DEL_HRO_pump'] / cfs_tafd,
       results['DEL_TRP_pump'] / cfs_tafd, 
       # (results['DEL_HRO_pump'] + results['DEL_TRP_pump']) / cfs_tafd,
       results['Combined_pump'] / cfs_tafd,
       results['SHA_storage'], 
       results['SHA_out'] / cfs_tafd,
       results['FOL_storage'],
       results['FOL_out'] / cfs_tafd,
       results['ORO_storage'],
       results['ORO_out'] / cfs_tafd,
       results['DEL_in'] / cfs_tafd,
       results['DEL_out'] / cfs_tafd,
       results['DEL_X2']]
  if plot:
    text_file = open("orca/figs/projection/datetime.txt", "w")
    text_file.write("%s\n" %now)
    text_file.write("%s" %sc)
    text_file.close()
    calibr_pts = ['HRO_pump','TRP_pump','Combined_pump','SHA_storage','SHA_out','FOL_storage','FOL_out','ORO_storage','ORO_out','DeltaIn','DeltaOut','X2']
    for f in ['D','W','M','AS-OCT']:
      for s,c in zip(sim,calibr_pts):
        plotter.plotting(s, freq=f)
        plt.savefig('orca/figs/projection/%s_%s.png' % (f,c), dpi=150)
        plt.close()  

if change_inflow_exeedance:
  write_json.modify('orca/data/json_files/SHA_properties.json', 'exceedance', SHA_exceedance_hist)
  write_json.modify('orca/data/json_files/ORO_properties.json', 'exceedance', ORO_exceedance_hist)
  write_json.modify('orca/data/json_files/FOL_properties.json', 'exceedance', FOL_exceedance_hist)

