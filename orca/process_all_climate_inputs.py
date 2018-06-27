import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import call
# from orca import *
#need climate data folders for this, which are too large for github (a few are present in repository for example)
with open('scenario_names.txt') as f:
	scenarios = f.read().splitlines()
i = 0
for sc in scenarios:
	i+=1
	print(i)
	call(['mkdir', 'scenario_runs/%s'%sc]) 
	call(['cp','input_climate_files/%s_input_data.csv'%sc,'scenario_runs/%s/%s_input_data.csv'%(sc,sc)]) #right now folder in jc- takes up a couple GB of space in github
	# call(['cp','input_climate_files/%s_input_data.csv'%sc,'orca/data/scenario_runs/%s/%s_input_data.csv'%(sc,sc)]) 
	call(['cp','scenario_runs/%s/%s_input_data.csv'%(sc,sc), 'climate_input_data.csv']) #right now folder in jc- takes up a couple GB of space
	call(['python', 'calc_indices_climate.py'])
	call(['cp','orca-data-processed-climate.csv', 'scenario_runs/%s/orca-data-processed-%s.csv'%(sc,sc)])
	call(['python', 'forecasting_climate.py'])
	call(['cp','scenario_runs/%s/orca-data-forecasted-%s.csv'%(sc,sc),'orca-data-climate-forecasted.csv'])
	# model = Model('orca/data/orca-data-climate-forecasted.csv', 'orca/data/results.csv',sd='10-01-1999',projection = True, sim_gains = True) #climate scenario test
	# results = model.simulate() # takes a while... save results
	# results.to_csv('orca/data/scenario_runs/%s/%s-results.csv'%(sc,sc))
