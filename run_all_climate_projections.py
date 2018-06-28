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
		proj_ind_df = process_projection(input_df,'orca/data/gains_regression.json') 
		proj_ind_df.to_csv('orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(sc,sc))
	if climate_forecasts:
		if not climate_indices:
			proj_ind_df = pd.read_csv('orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(sc,sc),index_col = 0, parse_dates = True)
		WYI_stats_file = pd.read_csv('orca/data/WYI_forcasting_regression_stats.csv', index_col = 0, parse_dates = True)
		carryover_stats_file = pd.read_csv('orca/data/carryover_regression_statistics.csv', index_col = 0, parse_dates = True)
		forc_df= projection_forecast(proj_ind_df,WYI_stats_file,carryover_stats_file)
		forc_df.to_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc))
	if run_projection: 
		model = Model('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc), 'orca/data/results.csv',sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
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