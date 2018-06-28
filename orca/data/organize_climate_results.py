import numpy as np
from subprocess import call
import pandas as pd
with open('scenario_names.txt') as f:
	scenarios = f.read().splitlines()
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
		dfobj = pd.read_csv('scenario_runs/%s/%s-results.csv'%(sc,sc), parse_dates = True, index_col = 0)
		df['%s'%sc] = dfobj[obj]
	df.to_csv('climate_results/%s.csv'%obj)