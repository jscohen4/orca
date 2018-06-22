import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import call
from orca import *

with open('orca/data/scenario_names.txt') as f:
	scenarios = f.read().splitlines()
i = 0
for sc in scenarios:
	i+=1
	print(i)
	# call(['mkdir', 'orca/data/scenrio_runs/%s'%sc]) 
	# call(['cp','input_climate_files/%s_input_data.csv'%sc,'orca/data/scenrio_runs/%s/%s_input_data.csv'%(sc,sc)]) #right now folder in jc- takes up a couple GB of space in github
	# call(['cp','input_climate_files/%s_input_data.csv'%sc,'orca/data/scenrio_runs/%s/%s_input_data.csv'%(sc,sc)]) 
	# call(['cp','orca/data/scenrio_runs/%s/%s_input_data.csv'%(sc,sc), 'orca/data/climate_input_data.csv']) #right now folder in jc- takes up a couple GB of space
	# call(['python', 'orca/data/calc_indices_climate.py'])
	# call(['cp','orca/data/orca-data-processed-climate.csv', 'orca/data/scenrio_runs/%s/orca-data-processed-%s.csv'%(sc,sc)])
	# call(['python', 'orca/data/forecasting_climate.py'])
	call(['cp','orca/data/scenrio_runs/%s/orca-data-forecasted-%s.csv'%(sc,sc),'orca/data/orca-data-climate-forecasted.csv'])
	model = Model('orca/data/orca-data-climate-forecasted.csv', 'orca/data/results.csv',sd='10-01-1999',scenario = True, sim_gains = True) #climate scenario test
	results = model.simulate() # takes a while... save results
	results.to_csv('orca/data/scenrio_runs/%s/%s-results.csv'%(sc,sc))
