import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import call
from orca import *
from orca.data import *

climate_indices = True
climate_forecasts = True
run_projection = True
consolidate_outputs = True
consolidate_inputs = False
#need climate data folders for this, which are too large for github (a few are present in repository for example)
with open('orca/data/scenario_names.txt') as f:
	scenarios = f.read().splitlines()
i = 0
for sc in scenarios:
	i+=1
	print('projection # %s' %i)
	call(['mkdir', 'orca/data/scenario_runs/%s'%sc])
	if climate_indices:
		input_df = pd.read_csv('orca/data/input_climate_files/%s_input_data.csv'%sc, index_col = 0, parse_dates = True)
		proj_ind_df = process_projection(input_df,'orca/data/json_files/gains_regression.json') 
		proj_ind_df.to_csv('orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(sc,sc))
	if climate_forecasts:
		if not climate_indices:
			proj_ind_df = pd.read_csv('orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(sc,sc),index_col = 0, parse_dates = True)
		WYI_stats_file = pd.read_csv('orca/data/forecast_regressions/WYI_forcasting_regression_stats.csv', index_col = 0, parse_dates = True)
		carryover_stats_file = pd.read_csv('orca/data/forecast_regressions/carryover_regression_statistics.csv', index_col = 0, parse_dates = True)
		forc_df= projection_forecast(proj_ind_df,WYI_stats_file,carryover_stats_file)
		forc_df.to_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc))
	if run_projection: 
		model = Model('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc), 'orca/data/historical_runs_data/results.csv',sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
		projection_results = model.simulate() # takes a while... save results
		projection_results.to_csv('orca/data/scenario_runs/%s/%s-results.csv'%(sc,sc))

if consolidate_outputs: 
	result_ids = ['SHA_storage','SHA_out','SHA_target','SHA_out_to_delta','SHA_tocs','FOL_storage','FOL_out',
								'FOL_target','FOL_out_to_delta','FOL_tocs','ORO_storage','ORO_out','ORO_target','ORO_out_to_delta',
								'ORO_tocs','DEL_in','DEL_out','DEL_TRP_pump','DEL_HRO_pump']
	for obj in result_ids:
		df = pd.DataFrame()
		print(obj)
		i = 0
		for sc in scenarios:
			i+=1
			print('projection # %s' %i)
			dfobj = pd.read_csv('orca/data/scenario_runs/%s/%s-results.csv'%(sc,sc), parse_dates = True, index_col = 0)
			df['%s'%sc] = dfobj[obj]
		df.to_csv('orca/data/climate_results/%s.csv'%obj)

if consolidate_inputs: 
	input_ids = ['TLG_fnf', 'FOL_fnf', 'MRC_fnf', 'MIL_fnf', 'NML_fnf', 'ORO_fnf',
       'MKM_fnf', 'BND_fnf', 'NHG_fnf', 'SHA_fnf', 'YRS_fnf', 'BKL_swe',
       'SHA_pr', 'ORO_pr', 'FOL_pr', 'SHA_tas', 'ORO_tas', 'FOL_tas',
       'SHA_tasmax', 'ORO_tasmax', 'FOL_tasmax', 'SHA_tasmin', 'ORO_tasmin',
       'FOL_tasmin', 'WY', 'DOWY', 'SR_WYI', 'SR_WYT', 'SR_WYT_rolling',
       'SJR_WYI', 'SJR_WYT', '8RI', 'SHA_fci', 'ORO_fci', 'FOL_fci', 'GOL_swe',
       'CSL_swe', 'HYS_swe', 'SCN_swe', 'RBB_swe', 'CAP_swe', 'RBP_swe',
       'HMB_swe', 'FOR_swe', 'RTL_swe', 'GRZ_swe', 'SDF_swe', 'SLT_swe',
       'MED_swe', 'BND_swe', 'ORO_swe', 'YRS_swe', 'FOL_swe', 'aprjul_slope',
       'aprjul_intercept', 'aprjul_mean', 'aprjul_std', 'octmar_mean',
       'octmar_std', 'octmar_intercept', 'octmar_slope', 'WYI_sim', 'WYT_sim',
       'gains_sim', 'SHA_snowpack', 'SHA_cum_flow_to_date',
       'SHA_remaining_flow', 'SHA_slope', 'SHA_intercept', 'SHA_mean',
       'SHA_std', 'ORO_snowpack', 'ORO_cum_flow_to_date', 'ORO_remaining_flow',
       'ORO_slope', 'ORO_intercept', 'ORO_mean', 'ORO_std', 'FOL_snowpack',
       'FOL_cum_flow_to_date', 'FOL_remaining_flow', 'FOL_slope',
       'FOL_intercept', 'FOL_mean', 'FOL_std']

	for obj in input_ids:
		df = pd.DataFrame()
		print(obj)
		i = 0
		for sc in scenarios:
			i+=1
			print('projection # %s' %i)
			dfobj = pd.read_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc), parse_dates = True, index_col = 0)
			df['%s'%sc] = dfobj[obj]
		df.to_csv('orca/data/climate_input_forecasts/%s.csv'%obj)		