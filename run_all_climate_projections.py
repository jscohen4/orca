import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import call
from orca import *
#need climate data folders for this, which are too large for github (a few are present in repository for example)
with open('orca/data/scenario_names.txt') as f:
	scenarios = f.read().splitlines()
i = 0
for sc in scenarios:
	i+=1
	print(i)
	model = Model('orca/data/orca-data-climate-forecasted.csv', 'orca/data/results.csv',sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
	results = model.simulate() # takes a while... save results
	results.to_csv('orca/data/scenario_runs/%s/%s-results.csv'%(sc,sc))

result_ids = ['SHA_storage','SHA_out','SHA_target','SHA_out_to_delta','SHA_tocs','FOL_storage','FOL_out',
							'FOL_target','FOL_out_to_delta','FOL_tocs','ORO_storage','ORO_out','ORO_target','ORO_out_to_delta',
							'ORO_tocs','DEL_in','DEL_out','DEL_TRP_pump','DEL_HRO_pump']
for obj in result_ids:
	df = pd.DataFrame()
	print(obj)
	i = 0
	for sc in scenarios:
		i+=1
		print(i)
		dfobj = pd.read_csv('orca/data/scenario_runs/%s/%s-results.csv'%(sc,sc), parse_dates = True, index_col = 0)
		df['%s'%sc] = dfobj[obj]
	df.to_csv('orca/data/climate_results/%s.csv'%obj)