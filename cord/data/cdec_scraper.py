import numpy as np
import pandas as pd 
from ulmo.cdec import historical as cd

# first get daily FNF
# first 8 are the 8-river index
# first 4 are the SRI and SR WYT
# http://cdec.water.ca.gov/cgi-progs/staSearch?sta=&sensor_chk=on&sensor=8
df = pd.DataFrame()
sd = '10-01-1999' # reliable start for CDEC daily data
ids = ['BND', 'ORO', 'YRS', 'FOL', 'NML', 'TLG', 'MRC', 'MIL']
       # 'GDW', 'MHB', 'MKM', 'NHG', 'SHA']
data = cd.get_data(station_ids=ids, sensor_ids=[8], 
                   resolutions=['daily'], start=sd)

for k in ids:
  df[k + '_fnf'] = data[k]['FULL NATURAL FLOW daily']['value']

# flowrates: inflow / outflow / pumping
ids = ['SHA', 'CLE', 'ORO', 'FOL', 'NML', 'DNP', 'EXC', 'MIL', 'SNL']

data = cd.get_data(station_ids=ids, sensor_ids=[15,23,76], 
                   resolutions=['daily'], start=sd)

for k in ids:
  df[k + '_in'] = data[k]['RESERVOIR INFLOW daily']['value']
  df[k + '_out'] = data[k]['RESERVOIR OUTFLOW daily']['value']
  df[k + '_storage'] = data[k]['RESERVOIR STORAGE daily']['value'] / 1000 # TAF

# observed delta outflow
data = cd.get_data(['DTO'], [23], ['daily'], start=sd)
df['DeltaOut'] = data['DTO']['RESERVOIR OUTFLOW daily']['value']

# banks & tracy pumping
ids = ['HRO', 'TRP']
data = cd.get_data(ids, [70], ['daily'], start=sd)
for k in ids:
  df[k + '_pump'] = data[k]['DISCHARGE, PUMPING daily']['value']

# estimate delta inflow from this (ignores GCD and direct precip)
df['DeltaIn'] = df['DeltaOut'] + df['HRO_pump'] + df['TRP_pump']

# cleanup
df[df < 0] = np.nan
df.interpolate(inplace=True)
df.to_csv('cord-data.csv')
