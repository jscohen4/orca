import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import call

with open('scenario_names.txt') as f:
	scenarios = f.read().splitlines()
i = 0
for sc in scenarios:
	i+=1
	print(i)
	# call(['mkdir', 'orca/data/scenrio_runs/%s'%sc]) 
	# call(['cp','input_climate_files/%s_input_data.csv'%sc,'orca/data/scenrio_runs/%s/%s_input_data.csv'%(sc,sc)]) #right now folder in jc- takes up a couple GB of space in github
	# call(['cp','input_climate_files/%s_input_data.csv'%sc,'orca/data/scenrio_runs/%s/%s_input_data.csv'%(sc,sc)]) 
	call(['cp','scenrio_runs/%s/%s_input_data.csv'%(sc,sc), 'climate_input_data.csv']) #right now folder in jc- takes up a couple GB of space
	call(['python', 'calc_indices_climate.py'])
	call(['cp','orca-data-processed-climate.csv', 'scenrio_runs/%s/orca-data-processed-%s.csv'%(sc,sc)])
	call(['python', 'forecasting_climate.py'])