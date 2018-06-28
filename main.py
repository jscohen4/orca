import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from orca import *
from subprocess import call
from datetime import datetime
now = datetime.now().strftime('Last run on %Y-%m-%d %H:%M:%S')

#Each of these booleans determines the actions that will be run by the model 

projection = True #True if running a single climate projection
sc = 'access1-0_rcp45_r1i1p1' #cmip5 climate scenario to use
calc_R2s = True #True if calculating R2s
plot = False #True if plotting outputs


process_hist_data = False #True if changing any historical data inputs
###Only relevant if processing historical data
cdec = True # True if downloading up-to-date cdec data
hist_indices = True #True if running calc_indices script
run_hist_forcast = True #True if running updated forecast

process_climate_data = True #only mark True if using running climate projection and processing projection input data
####### only relevant if processing projection data
climate_indices = True
climate_forecasts = True

####################### below this line data processing and model runs are executed
if process_hist_data: 
  from orca.data import *
  if cdec:
    cdec_df = cdec_scraper.scrape_cdec()
    cdec_df.to_csv('orca/data/cdec-data.csv')
  if hist_indices:
    if not cdec: 
      cdec_df = pd.read_csv('orca/data/cdec-data.csv', index_col=0, parse_dates=True)
    ind_df = process(cdec_df,'orca/data/evap_regression.json','orca/data/gains_regression.json')  
    ind_df.to_csv('orca/data/orca-data-processed.csv')
  if run_hist_forcast:
    if not hist_indices:
      ind_df = pd.read_csv('orca/data/orca-data-processed.csv', index_col=0, parse_dates=True)
    forc_df, stats_df, WYI_stats= forecast(ind_df)
    forc_df.to_csv('orca/data/orca-data-forecasted.csv')
    stats_df.to_csv('orca/data/carryover_regression_statistics.csv')
    WYI_stats.to_csv('orca/data/WYI_forcasting_regression_stats.csv')
  print('done')
if not projection:
  model = Model('orca/data/orca-data-forecasted.csv', 'orca/data/orca-data-forecasted.csv',sd='10-01-1999',projection = False, sim_gains = False) #beacuse of rolling calc in gains, we start on 10th day of
  results = model.simulate() # takes a while... save results
  results.to_csv('orca/data/results.csv')
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
         results['DEL_out'] / cfs_tafd]
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
           model.df['DeltaOut']]
    plotter.Rsquares(sim,obs)
    if plot:
      calibr_pts = ['HRO_pump','TRP_pump','Combined_pump','SHA_storage','SHA_out','FOL_storage','FOL_out','ORO_storage','ORO_out','DeltaIn','DeltaOut']
      for f in ['D','W','M','AS-OCT']:
        for s,o,c in zip(sim,obs,calibr_pts):
          plotter.compare(s, o, freq=f)
          plt.savefig('orca/figs/%s_%s_hist.pdf' % (f,c), dpi=150)
          plt.close()  
    print('done')
if process_climate_data:
  from orca.data import *
  call(['mkdir', 'orca/data/individual_projection_runs/%s'%sc])
  text_file = open("orca/data/individual_projection_runs/%s/datetime.txt"%sc, "w")
  text_file.write("%s" %now)
  text_file.close()

  call(['cp','orca/data/input_climate_files/%s_input_data.csv'%sc,'orca/data/individual_projection_runs/%s/%s_input_data.csv'%(sc,sc)]) 
  if climate_indices:
    input_df = pd.read_csv('orca/data/individual_projection_runs/%s/%s_input_data.csv'%(sc,sc), index_col = 0, parse_dates = True)
    proj_ind_df = process_projection(input_df,'orca/data/gains_regression.json')  
    proj_ind_df.to_csv('orca/data/individual_projection_runs/%s/orca-data-processed-%s.csv'%(sc,sc))
  if climate_forecasts:
    if not climate_indices:
      proj_ind_df = pd.read_csv('orca/data/individual_projection_runs/%s/orca-data-processed-%s.csv'%(sc,sc), index_col = 0, parse_dates = True)
    WYI_stats_file = pd.read_csv('orca/data/WYI_forcasting_regression_stats.csv', index_col = 0, parse_dates = True)
    carryover_stats_file = pd.read_csv('orca/data/carryover_regression_statistics.csv', index_col = 0, parse_dates = True)
    forc_df= projection_forecast(proj_ind_df,WYI_stats_file,carryover_stats_file)
    forc_df.to_csv('orca/data/individual_projection_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc))
  print('done')
if projection:
  model = Model('orca/data/individual_projection_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc), 'orca/data/results.csv',sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
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
       results['DEL_out'] / cfs_tafd]
  if plot:
    calibr_pts = ['HRO_pump','TRP_pump','Combined_pump','SHA_storage','SHA_out','FOL_storage','FOL_out','ORO_storage','ORO_out','DeltaIn','DeltaOut']
    for f in ['D','W','M','AS-OCT']:
      for s,c in zip(sim,calibr_pts):
        plotter.plotting(s, freq=f)
        plt.savefig('orca/figs/%s_%s_proj.pdf' % (f,c), dpi=150)
        plt.close()  

