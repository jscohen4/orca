import numpy as np
np.warnings.filterwarnings('ignore') #to not display numpy warnings... be careful
import pandas as pd
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import call
from orca import *
from orca.data import *
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# this whole script will run on all processors requested by the job script
with open('orca/data/scenario_names_all.txt') as f:
	scenarios = f.read().splitlines()
result_ids = ['SHA_storage','SHA_out','SHA_target','SHA_out_to_delta','SHA_tocs','FOL_storage','FOL_out',
				      'FOL_target','FOL_out_to_delta','FOL_tocs','ORO_storage','ORO_out','ORO_target','ORO_out_to_delta',
				      'ORO_tocs','DEL_in','DEL_out','DEL_TRP_pump','DEL_HRO_pump','SHA_sodd','SHA_spill', 
					'ORO_sodd','ORO_spill','FOL_sodd','FOL_spill', 'DEL_X2','ORO_forecast','FOL_forecast','SHA_forecast','DEL_SODD_CVP','DEL_SODD_SWP']

input_ids = ['SHA_in_tr', 'ORO_in_tr','FOL_in_tr','ORO_fnf','BND_fnf', 'SHA_fnf','FOL_fnf','SR_WYI','SR_WYT_rolling','SR_WYT', 'WYT_sim','WYI_sim',
        'SHA_fci', 'ORO_fci', 'FOL_fci', 'BND_swe', 'ORO_swe', 'YRS_swe', 'FOL_swe','SHA_remaining_flow','ORO_remaining_flow','FOL_remaining_flow',
        'octmar_flow_to_date','octmar_mean','octmar_std','aprjul_flow_to_date','aprjul_mean','aprjul_std','aprjul_slope','aprjul_intercept']

comm = MPI.COMM_WORLD # communication object
rank = comm.rank # what number processor am I?
s = scenarios[rank] 

#Define some parameters and adaptations
window_type = 'expanding'
window_length = 40
index_exceedence_sac = 0
shift = 0
SHA_shift = shift
ORO_shift = shift
FOL_shift = shift

SHA_exceedance = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}
FOL_exceedance = {"W": 10, "AN": 10, "BN": 5, "D": 2, "C": 1}
ORO_exceedance = {"W": 2, "AN": 2, "BN": 2, "D": 2, "C": 2}

#####Which files should stay on cluster
remove_processed = True
remove_forecasted = True
remove_model_results = False


call(['mkdir', 'orca/data/scenario_runs/%s'%s])
input_df = pd.read_csv('input_climate_files/%s_input_data.csv'%s, index_col = 0, parse_dates = True)
gains_loop_df = pd.read_csv('orca/data/historical_runs_data/gains_loops.csv', index_col = 0, parse_dates = True)
OMR_loop_df = pd.read_csv('orca/data/historical_runs_data/OMR_loops.csv', index_col = 0, parse_dates = True)
proj_ind_df = process_projection(input_df,gains_loop_df,OMR_loop_df,'orca/data/json_files/gains_regression.json','orca/data/json_files/inf_regression.json',window = window_type) 
proj_ind_df.to_csv('orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(s,s))
proj_ind_df = pd.read_csv('orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(s,s), index_col = 0, parse_dates = True)
WYI_stats_file = pd.read_csv('orca/data/forecast_regressions/WYI_forcasting_regression_stats.csv', index_col = 0, parse_dates = True)
carryover_stats_file = pd.read_csv('orca/data/forecast_regressions/carryover_regression_statistics.csv', index_col = 0, parse_dates = True)


forc_df= projection_forecast(proj_ind_df,WYI_stats_file,carryover_stats_file,window_type,window_length, index_exceedence_sac)
forc_df.to_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(s,s))

model = Model('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(s), 
		'orca/data/historical_runs_data/results.csv',SHA_shift, ORO_shift, FOL_shift,sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
projection_results = model.simulate() # takes a while... save results
projection_results.to_csv('orca/data/scenario_runs/%s/%s-results.csv'%(s,s))
comm.barrier()
if consolidate_outputs:
	if comm.rank <= 30: 
		obj = result_ids[comm.rank]
		dfobj = pd.DataFrame()
		for sc in scenarios: 	
			projection_results = pd.read_csv('orca/data/scenario_runs/%s/%s-results.csv'%(sc,sc), index_col = 0, parse_dates = True)
			dfobj[sc] = projection_results[obj]
		dfobj.to_csv('orca/data/climate_results/%s.csv'%(obj))

if consolidate_outputs:
	if comm.rank >= 31 and comm.rank <=59: 
		obj = input_ids[comm.rank-31]
		dfobj = pd.DataFrame()
		for sc in scenarios:
			projection_results = pd.read_csv('orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(sc,sc), index_col = 0, parse_dates = True)
			dfobj[sc] = projection_results[obj]
		dfobj.to_csv('orca/data/climate_input_forecasts/%s.csv'%(obj))
comm.barrier()	
if remove_processed:
	call(['rm', 'orca/data/scenario_runs/%s/orca-data-processed-%s.csv'%(s,s)])
if remove_forecasted:
	call(['rm', 'orca/data/scenario_runs/%s/orca-data-climate-forecasted-%s.csv'%(s,s)])
if remove_model_output:
	call(['rm', 'orca/data/scenario_runs/%s/%s-results.csv'%(s,s)])




