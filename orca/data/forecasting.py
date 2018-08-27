import numpy as np 
import pandas as pd
from sklearn import linear_model
from .write_json import modify
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
cfs_tafd = 2.29568411*10**-5 * 86400 / 1000
pd.options.mode.chained_assignment = None  # default='warn'
from .util import *
# import matplotlib.pyplot as plt
def WYI_to_WYT(WYI, thresholds, values):
  for t,v in zip(thresholds,values):
    if WYI > t:
      return v

def octmar_cumulative(x): #cumulative inflow from october through march

	ix = (x.index.month >= 10) | (x.index.month <= 3)
	octmar = (x[ix].sum() - x[ix].cumsum())
	ix = (x.index.month >= 4) & (x.index.month <= 7)
	x[ix] = octmar[5]
	aprjul = x[ix]
	return pd.concat([octmar,aprjul])

	# ix = (x.index.month >= 10) | (x.index.month <= 3)
	# octmar = (x[ix].sum() - x[ix].cumsum())
	# ix = (x.index.month >= 4) & (x.index.month <= 7)
	# x[ix] = octmar[5]
	# aprjul = x[ix]
	# return pd.concat([octmar,aprjul])

def aprjul_cumulative(x): #cumulative inflow from april through july
	ix = (x.index.month >= 4) & (x.index.month <= 7)
	aprjul = (x[ix].sum() - x[ix].cumsum())
	ix = (x.index.month >= 10) | (x.index.month <= 3)
	x[ix] = aprjul[0]
	octmar = x[ix]
	ix = (x.index.month == 8) | (x.index.month ==9)
	x[ix] = aprjul[3]
	augsep = x[ix]

	return pd.concat([octmar,aprjul,augsep]) 

def octmar_flow_to_date(x):
	ix = (x.index.month >= 10) | (x.index.month <= 3)
	octmar = x[ix].cumsum()
	ix = (x.index.month >= 4) & (x.index.month <= 7)
	x[ix] = octmar[5]
	aprjul =  x[ix]
	ix = (x.index.month == 8) | (x.index.month == 9)
	x[ix] = octmar[5]
	augsep = x[ix]
	return pd.concat([octmar,aprjul,augsep])

def aprjul_flow_to_date(x):
	ix = (x.index.month >= 4) & (x.index.month <= 7)
	aprjul =  x[ix].cumsum()
	return aprjul

def flow_to_date(x):
	# ix = (x.index.month >= 1) 
	cum_flow = x.cumsum()
	return cum_flow

def rem_flow(x):
	# ix = (x.index.month >= 1)
	remaining_flow = (x.sum() - x.cumsum())
	return remaining_flow

def get_forecast_WYI(df, index_exceedance_sac): #now determining forecasting regression coefficients based off perfect foresight
	flow_sites = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
	snow_sites = ['BND_swe', 'ORO_swe', 'YRS_swe', 'FOL_swe']

	Qm = (df[flow_sites].sum(axis=1).resample('MS').sum().to_frame(name='flow') * cfsd_mafd)

	Qm['WY'] = df.WY.resample('MS').first() # need this for grouping
	Qm['snow'] = df[snow_sites].sum(axis=1).resample('MS').first()
	Qm = Qm[Qm.WY != Qm.WY[-1]] 
	Qm['octmar_cumulative'] = (Qm.groupby('WY').flow
	                          .apply(octmar_cumulative)
	                          .reset_index(level=0)
	                          .drop('WY', axis=1)) 

	Qm['aprjul_cumulative'] = (Qm.groupby('WY').flow
	                           .apply(aprjul_cumulative))

	Qm['octmar_flow_to_date'] = (Qm.groupby('WY').flow
	                           .apply(octmar_flow_to_date))

	Qm['aprjul_flow_to_date'] = (Qm.groupby('WY').flow
	                           .apply(aprjul_flow_to_date)
	                           .reset_index(level=0).drop('WY', axis=1)) 

#obtaining regression coefficients for WYI forecasting
	aprjul_slopes = np.zeros(12)
	aprjul_intercepts = np.zeros(12)
	aprjul_means = np.zeros(12)
	aprjul_stds = np.zeros(12)
	octmar_means = np.zeros(12)
	octmar_stds = np.zeros(12)
	octmar_slopes = np.zeros(12)
	octmar_intercepts = np.zeros(12)
	month = Qm.index.month
	snow_ind = Qm.snow
	aprjul_cumulative_ind = Qm.aprjul_cumulative
	for m in range(1,13):
		oct_index = (month == 10)
		ix = (month == m)
		coeffs = np.polyfit(snow_ind[ix],aprjul_cumulative_ind[ix],1)
		aprjul_slopes[m-1] = coeffs[0]
		aprjul_intercepts[m-1] = coeffs[1]
		aprjul_means[m-1] = np.mean(aprjul_cumulative_ind[ix]) 
		aprjul_stds[m-1] =  np.std(aprjul_cumulative_ind[ix])
	octmar_flow_to_date_ind = Qm.octmar_flow_to_date	
	octmar_cumulative_ind = Qm.octmar_cumulative
	for m in list(range(10,13))+list(range(1,3)):
		ix = (month == m)
		flow_coeffs = np.polyfit(octmar_flow_to_date_ind[oct_index], octmar_cumulative_ind[ix], 1)
		octmar_slopes[m-1] = flow_coeffs[0]
		octmar_intercepts[m-1] = flow_coeffs[1]
		octmar_means[m-1] = np.mean(octmar_cumulative_ind[ix])
		octmar_stds[m-1] = np.std(octmar_cumulative_ind[ix])
	stats =  {'aprjul_slope': aprjul_slopes,
                    'aprjul_intercept': aprjul_intercepts, 
                    'aprjul_mean': aprjul_means,
                    'aprjul_std':  aprjul_stds, 
                    'octmar_mean': octmar_means,
                    'octmar_std':  octmar_stds, 
                    'octmar_intercept': octmar_intercepts, 
                    'octmar_slope': octmar_slopes}

	stats = pd.DataFrame(stats, columns = ['aprjul_slope',
                    'aprjul_intercept', 
                    'aprjul_mean',
                    'aprjul_std', 
                    'octmar_mean',
                    'octmar_std', 
                    'octmar_intercept', 
                    'octmar_slope'])
	# stats.to_csv('WYI_forcasting_regression_stats.csv')
	for s in stats: 
		df[s] = pd.Series(index=df.index)
		for m in range(1,13):
			Qm.loc[month==m, s] = stats[s][m-1]

	Qm = Qm.fillna(0)
	Qm['WYI'] = pd.Series(index=Qm.index)
	prev = 10
	#using water year index formula
	for index, row in Qm.iterrows():
		ix = index.month
		if (ix == 10) | (ix == 11):
			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'] + Qm.loc[index, 'octmar_std']*z_table_transform[index_exceedance_sac])\
			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + Qm.loc[index, 'aprjul_mean'] + Qm.loc[index, 'aprjul_std']*z_table_transform[index_exceedance_sac])
		elif (ix == 12) | (ix <= 3):
			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'] + Qm.loc[index, 'octmar_std']*z_table_transform[index_exceedance_sac])\
			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] +  (Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept']) + Qm.loc[index, 'aprjul_std']*z_table_transform[index_exceedance_sac])
		
		elif (ix == 4) | (ix == 5): 
			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'])\
			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + (Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept']) + Qm.loc[index, 'aprjul_std']*z_table_transform[index_exceedance_sac])
			WYI_rem = WYI
		elif (ix <= 9 & ix >= 6):
			WYI = WYI_rem
			prev = min(WYI_rem,10)
		Qm.loc[index, 'WYI'] = WYI

	# Qm.WYI = Qm.WYI.shift(periods=-1)
	return(Qm.WYI,stats)

def forecast(df,index_exceedance_sac):
	WYI_sim,WYI_stats = get_forecast_WYI(df,index_exceedance_sac) #wyt
	df['WYI_sim'] = WYI_sim
	df.WYI_sim = df.WYI_sim.fillna(method = 'bfill')
	df.loc[df['WYI_sim'].isnull(),'WYI_sim'] = df['SR_WYI']

	df['WYT_sim'] = df.WYI_sim.apply(WYI_to_WYT,
	                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
	                               values=['W', 'AN', 'BN', 'D', 'C'])

	###making forcasts for reservoir carryover:

	snow_sites = ['BND_swe', 'ORO_swe','FOL_swe']
	res_ids = ['SHA','ORO','FOL']

	SHA_inf = (df['SHA_in_tr'].to_frame(name='inf') * cfs_tafd)  
	ORO_inf = (df['ORO_in_tr'].to_frame(name='inf') * cfs_tafd)
	FOL_inf = (df['FOL_in_tr'].to_frame(name='inf') * cfs_tafd)


	res_frames = [SHA_inf,ORO_inf,FOL_inf]
	stats_file = pd.DataFrame()
	#obtaining daily inflow forecast regression coeffients
	for r, swe, res_id in zip(res_frames, snow_sites, res_ids):
		r['WY'] = df.WY
		r['DOWY'] = df.DOWY
		r['snowpack'] = df[swe]
		r = r[r.WY != r.WY[-1]]
		month = r.index.month

		r['cum_flow_to_date'] = (r.groupby('WY').inf
		                           .apply(flow_to_date))
		r['remaining_flow'] = (r.groupby('WY').inf
		                            .apply(rem_flow))
		r.cum_flow_to_date.fillna(method='ffill', inplace=True)
		slopes = np.zeros(365)
		intercepts = np.zeros(365)
		means = np.zeros(365)
		stds = np.zeros(365)
		snowpack = r.snowpack
		remaining_flow = r.remaining_flow
		DOWY = r.DOWY
		for d in range(1,365):
			ix = (DOWY == d)
			coeffs = np.polyfit(snowpack[ix],remaining_flow[ix],1)
			slopes[d] = coeffs[0]
			intercepts[d] = coeffs[1]
			means[d] = np.mean(remaining_flow[ix]) 
			stds[d] = np.std(remaining_flow[ix])
		

		res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]

		stats =  {res_stats[0]: slopes,
		                  res_stats[1]: intercepts, 
		                  res_stats[2]: means,
		                  res_stats[3]: stds}


		stats = pd.DataFrame(stats)

		if res_id == 'SHA':
			stats_file = stats
		else:
			stats_file[res_stats] = stats

		stats = stats.values.T
		for i,s in enumerate(stats):
			yearly_stats = stats[i]
			v = np.append(yearly_stats,np.tile(yearly_stats, 3)) #1997-2000 WY
			v = np.append(v,[yearly_stats[364]]) #leap year
			for y in range(4): # 14 more years
				v = np.append(v,np.tile(yearly_stats, 4))
				v = np.append(v,[yearly_stats[364]]) #leap year
			v = np.append(v,np.tile(yearly_stats, 1)) #2017 WY
			r[res_stats[i]] = pd.Series(v,index=r.index)
		r.rename(columns = {'cum_flow_to_date':'%s_cum_flow_to_date'%res_id}, inplace=True)
		r.rename(columns = {'remaining_flow':'%s_remaining_flow'%res_id}, inplace=True)
		r.rename(columns = {'snowpack':'%s_snowpack'%res_id}, inplace=True)
		r.drop(['inf','WY','DOWY'], axis=1, inplace=True)
		df = pd.concat([df, r], axis=1, join_axes=[df.index])
	return df,stats_file,WYI_stats



#############climate projection functions
def get_projection_forecast_WYI(df, stats_file,index_exceedance_sac): #now determining forecasting regression coefficients based off perfect foresight
	flow_sites = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
	snow_sites = ['BND_swe', 'ORO_swe', 'YRS_swe', 'FOL_swe']

	Qm = (df[flow_sites].sum(axis=1).resample('M').sum().to_frame(name='flow') * cfsd_mafd)

	Qm['WY'] = df.WY.resample('M').first() # need this for grouping

	Qm['snow'] = df[snow_sites].sum(axis=1).resample('M').first()

	Qm = Qm[Qm.WY != Qm.WY[-1]] 

	Qm['octmar_cumulative'] = (Qm.groupby('WY').flow
	                          .apply(octmar_cumulative)
	                          .reset_index(level=0)
	                          .drop('WY', axis=1)) 
	Qm['aprjul_cumulative'] = (Qm.groupby('WY').flow
	                           .apply(aprjul_cumulative))

	Qm['octmar_flow_to_date'] = (Qm.groupby('WY').flow
	                           .apply(octmar_flow_to_date))

	Qm['aprjul_flow_to_date'] = (Qm.groupby('WY').flow
	                           .apply(aprjul_flow_to_date)
	                           .reset_index(level=0).drop('WY', axis=1)) 

	month = Qm.index.month
	# snow_ind = Qm.snow
	# aprjul_cumulative_ind = Qm.aprjul_cumulative
	# octmar_flow_to_date_ind = Qm.octmar_flow_to_date	
	# octmar_cumulative_ind = Qm.octmar_cumulative
	# stats = pd.read_csv('WYI_forcasting_regression_stats.csv', index_col=0, parse_dates=True)
	for s in stats_file: 
		for m in range(1,13):
			Qm.loc[month==m, s] = stats_file[s][m-1]
  #Qm.octmar_cumulative = Qm.octmar_cumulative.shift(periods=1)
  #Qm.aprjul_cumulative = Qm.aprjul_cumulative.shift(periods=1)
	Qm = Qm.fillna(0)
	Qm['WYI'] = pd.Series(index=Qm.index)
	prev = 10
	for index, row in Qm.iterrows():
		ix = index.month
		if (ix == 10) | (ix == 11):
			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'] + Qm.loc[index, 'octmar_std']*z_table_transform[index_exceedance_sac])\
			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + Qm.loc[index, 'aprjul_mean'] + Qm.loc[index, 'aprjul_std']*z_table_transform[index_exceedance_sac])
		elif (ix == 12) | (ix <= 3):
			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'] + Qm.loc[index, 'octmar_mean'] + Qm.loc[index, 'octmar_std']*z_table_transform[index_exceedance_sac])\
			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] +  (Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept']) + Qm.loc[index, 'aprjul_std']*z_table_transform[index_exceedance_sac])
		
		elif (ix == 4) | (ix == 5): 
			WYI = 0.3 * prev + 0.3 * (Qm.loc[index, 'octmar_flow_to_date'])\
			+ 0.4 * (Qm.loc[index, 'aprjul_flow_to_date'] + (Qm.loc[index, 'aprjul_slope'] * Qm.loc[index, 'snow'] + Qm.loc[index,'aprjul_intercept']) + Qm.loc[index, 'aprjul_std']*z_table_transform[index_exceedance_sac])
			WYI_rem = WYI

		if (ix <= 9 & ix>=6):
			WYI = WYI_rem
			prev = min(WYI_rem,10)
		Qm.loc[index, 'WYI'] = WYI

	return Qm.WYI, Qm[['octmar_flow_to_date','octmar_mean','octmar_std','aprjul_flow_to_date','aprjul_mean','aprjul_std','aprjul_slope','aprjul_intercept']]

def projection_forecast(df,WYI_stats_file,carryover_stats_file,window_type,window_length,index_exceedance_sac):
	# WYI_sim = get_projection_forecast_WYI(df,WYI_stats_file,index_exceedance_sac) #wyt
	# decade_thresh = ['1999-10-01','2009-10-01','2019-10-01','2029-10-01','2039-10-01','2049-10-01','2059-10-01','2069-10-01',
	# '2079-10-01','2089-10-01','2099-12-31']
	if window_type == 'historical':
		WYI_stats = get_projection_forecast_WYI(df,WYI_stats_file,index_exceedance_sac)
		WYI_sim = WYI_stats[0]
		WYI_stats = WYI_stats[1]

		snow_sites = ['BND_swe', 'ORO_swe','FOL_swe']
		res_ids = ['SHA','ORO','FOL']
		SHA_inf = (df['SHA_in_tr'].to_frame(name='inf') * cfs_tafd)  
		ORO_inf = (df['ORO_in_tr'].to_frame(name='inf') * cfs_tafd)
		FOL_inf = (df['FOL_in_tr'].to_frame(name='inf') * cfs_tafd)
		res_frames = [SHA_inf,ORO_inf,FOL_inf]

		for r, swe, res_id in zip(res_frames, snow_sites, res_ids):
			r['WY'] = df.WY
			r['DOWY'] = df.DOWY
			r['snowpack'] = df[swe]
			r = r[r.WY != r.WY[-1]]
			month = r.index.month
			r['cum_flow_to_date'] = (r.groupby('WY').inf
			                           .apply(flow_to_date))
			r['remaining_flow'] = (r.groupby('WY').inf
			                            .apply(rem_flow))
			r.cum_flow_to_date.fillna(method='ffill', inplace=True)
			res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]
			# stats = pd.read_csv('carryover_regression_statistics.csv', index_col = 0)
			carryover_stats = carryover_stats_file[res_stats]
			carryover_stats = carryover_stats.values.T
			for i,s in enumerate(carryover_stats):
				yearly_stats = carryover_stats[i]
				v = np.append(yearly_stats,np.tile(yearly_stats, 1)) #2000 WY
				v = np.append(v,[yearly_stats[364]]) #leap year
				for y in range(36): # 2001-2096 WYs
					v = np.append(v,np.tile(yearly_stats, 4))
					v = np.append(v,[yearly_stats[364]]) #leap year
				v = np.append(v,np.tile(yearly_stats, 2)) #2097-2099 WY
				r[res_stats[i]] = pd.Series(v,index=r.index)
			r.rename(columns = {'cum_flow_to_date':'%s_cum_flow_to_date'%res_id}, inplace=True)
			r.rename(columns = {'remaining_flow':'%s_remaining_flow'%res_id}, inplace=True)
			r.rename(columns = {'snowpack':'%s_snowpack'%res_id}, inplace=True)
			r.drop(['inf','WY','DOWY'], axis=1, inplace=True)
			df = pd.concat([df, r], axis=1, join_axes=[df.index])
		df = pd.concat([df,WYI_stats], axis=1, join_axes=[df.index])

	if window_type == 'stationary':
		decade_thresh = pd.date_range('1951-10-01','2099-10-01',freq =' AS-OCT')
		WYI_mov_stats = get_forecast_WYI_stats(df.truncate(before = decade_thresh[2], after=decade_thresh[50]),index_exceedance_sac)
		WYI_sim = get_projection_forecast_WYI(df,WYI_mov_stats,index_exceedance_sac)

		WYI_stats = WYI_sim[1]
		WYI_sim = WYI_sim[0]
		snow_sites = ['BND_swe', 'ORO_swe','FOL_swe']
		res_ids = ['SHA','ORO','FOL']
		SHA_inf = (df['SHA_in_tr'].to_frame(name='inf') * cfs_tafd)  
		ORO_inf = (df['ORO_in_tr'].to_frame(name='inf') * cfs_tafd)
		FOL_inf = (df['FOL_in_tr'].to_frame(name='inf') * cfs_tafd)
		res_frames = [SHA_inf,ORO_inf,FOL_inf]

		for r, swe, res_id in zip(res_frames, snow_sites, res_ids):
			r['WY'] = df.WY
			r['DOWY'] = df.DOWY
			r['snowpack'] = df[swe]
			r = r[r.WY != r.WY[-1]]
			month = r.index.month
			r['cum_flow_to_date'] = (r.groupby('WY').inf
			                           .apply(flow_to_date))
			r['remaining_flow'] = (r.groupby('WY').inf
			                            .apply(rem_flow))
			r.cum_flow_to_date.fillna(method='ffill', inplace=True)
			res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]

			carryover_stats = carryover_stats_file[res_stats]
			carryover_stats = carryover_stats.values.T
			for i,s in enumerate(carryover_stats):
				yearly_stats = carryover_stats[i]
				v = np.append(yearly_stats,np.tile(yearly_stats, 1)) #2000 WY
				v = np.append(v,[yearly_stats[364]]) #leap year
				for y in range(36): # 2001-2096 WYs
					v = np.append(v,np.tile(yearly_stats, 4))
					v = np.append(v,[yearly_stats[364]]) #leap year
				v = np.append(v,np.tile(yearly_stats, 2)) #2097-2099 WY
				r[res_stats[i]] = pd.Series(v,index=r.index)
			r.rename(columns = {'cum_flow_to_date':'%s_cum_flow_to_date'%res_id}, inplace=True)
			r.rename(columns = {'remaining_flow':'%s_remaining_flow'%res_id}, inplace=True)
			r.rename(columns = {'snowpack':'%s_snowpack'%res_id}, inplace=True)
			r.drop(['inf','WY','DOWY'], axis=1, inplace=True)
			df = pd.concat([df, r], axis=1, join_axes=[df.index])
		df = pd.concat([df,WYI_stats], axis=1, join_axes=[df.index])


	elif window_type == 'rolling':

		decade_thresh = pd.date_range('1951-10-01','2099-10-01',freq =' AS-OCT')

		WYI_mov_stats = get_forecast_WYI_stats(df.truncate(before = decade_thresh[2], after=decade_thresh[50]),index_exceedance_sac)
		WYI_sim = get_projection_forecast_WYI(df.truncate(before=decade_thresh[0], after=decade_thresh[50]),WYI_mov_stats,index_exceedance_sac)
		WYI_stats = WYI_sim[1]
		WYI_sim = WYI_sim[0]

		for i in range(50,148):
			WYI_mov_stats = get_forecast_WYI_stats(df.truncate(before = decade_thresh[i-window_length], after=decade_thresh[i]),index_exceedance_sac)
			WYI_dec = get_projection_forecast_WYI(df.truncate(before=decade_thresh[i], after=decade_thresh[i+1]),WYI_mov_stats,index_exceedance_sac)
			WYI_sim = pd.concat([WYI_sim,WYI_dec[0]])
			WYI_stats = pd.concat([WYI_stats, WYI_dec[1]])
		snow_sites = ['BND_swe', 'ORO_swe','FOL_swe']
		res_ids = ['SHA','ORO','FOL']
		SHA_inf = (df['SHA_in_tr'].to_frame(name='inf') * cfs_tafd)  
		ORO_inf = (df['ORO_in_tr'].to_frame(name='inf') * cfs_tafd)
		FOL_inf = (df['FOL_in_tr'].to_frame(name='inf') * cfs_tafd)
		res_frames = [SHA_inf,ORO_inf,FOL_inf]

		for r, swe, res_id in zip(res_frames, snow_sites, res_ids):
			rstats = pd.DataFrame()
			r['WY'] = df.WY
			r['DOWY'] = df.DOWY
			r['snowpack'] = df[swe]
			r = r[r.WY != r.WY[-1]]
			month = r.index.month

			r['cum_flow_to_date'] = (r.groupby('WY').inf
			                           .apply(flow_to_date))
			r['remaining_flow'] = (r.groupby('WY').inf
			                            .apply(rem_flow))
			r.cum_flow_to_date.fillna(method='ffill', inplace=True)

			rstart = r.truncate(before = decade_thresh[2], after=decade_thresh[50])
			rstart = rstart[(rstart.index < decade_thresh[50])]
			slopes = np.zeros(365)
			intercepts = np.zeros(365)
			means = np.zeros(365)
			stds = np.zeros(365)
			snowpack = rstart.snowpack + 0.0001
			remaining_flow = rstart.remaining_flow + 0.0001
			DOWY = rstart.DOWY
			for d in range(1,365):
				ix = (DOWY == d)
				coeffs = np.polyfit(snowpack[ix],remaining_flow[ix],1)
				slopes[d] = coeffs[0]
				intercepts[d] = coeffs[1]
				means[d] = np.mean(remaining_flow[ix]) 
				stds[d] = np.std(remaining_flow[ix])
			

			res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]
			stats =  {res_stats[0]: slopes,
			                  res_stats[1]: intercepts, 
			                  res_stats[2]: means,
			                  res_stats[3]: stds}
			stats = pd.DataFrame(stats)
			if res_id == 'SHA':
				stats_file = stats
			else:
				stats_file[res_stats] = stats
			stats = stats.values.T
			rstatsind = r.truncate(before = decade_thresh[0], after=decade_thresh[50])
			rstatsind = rstatsind[(rstatsind.index < decade_thresh[50])].index
			res_slope = []
			res_int = []
			res_mean = []
			res_std = []

			for i,s in enumerate(stats):
				yearly_stats = stats[i]
				if i == 0:
					res_slope = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_slope = np.append(res_slope,np.tile(yearly_stats, 4))
						res_slope = np.append(res_slope,[yearly_stats[364]]) #leap year
					res_slope = np.append(res_slope,np.tile(yearly_stats, 1))
				
				elif i == 1:
					res_intercept = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_intercept = np.append(res_intercept,np.tile(yearly_stats, 4))
						res_intercept = np.append(res_intercept,[yearly_stats[364]]) #leap year
					res_intercept = np.append(res_intercept,np.tile(yearly_stats, 1))
				
				elif i == 2:
					res_mean = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_mean = np.append(res_mean,np.tile(yearly_stats, 4))
						res_mean = np.append(res_mean,[yearly_stats[364]]) #leap year
					res_mean = np.append(res_mean,np.tile(yearly_stats, 1))

				elif i == 3:
					res_std = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_std = np.append(res_std,np.tile(yearly_stats, 4))
						res_std = np.append(res_std,[yearly_stats[364]]) #leap year
					res_std = np.append(res_std,np.tile(yearly_stats, 1))

			for y in range(50,148):
				if window_length >= 49 & ((y == 50) | (y ==51)):
					wl = 48
				else:
					wl = window_length 
				rnext = r.truncate(before = decade_thresh[y-wl], after=decade_thresh[y])
				rnext = rnext[(rnext.index < decade_thresh[y])]
				rnow_ind = r.truncate(before = decade_thresh[y], after=decade_thresh[y+1])
				rnow_ind = rnow_ind[(rnow_ind.index < decade_thresh[y+1])].index
				slopes = np.zeros(365)
				intercepts = np.zeros(365)
				means = np.zeros(365)
				stds = np.zeros(365)
				snowpack = rnext.snowpack + 0.0001
				remaining_flow = rnext.remaining_flow + 0.0001
				DOWY = rnext.DOWY
				for d in range(1,365):
					ix = (DOWY == d)
					coeffs = np.polyfit(snowpack[ix],remaining_flow[ix],1)
					slopes[d] = coeffs[0]
					intercepts[d] = coeffs[1]
					means[d] = np.mean(remaining_flow[ix]) 
					stds[d] = np.std(remaining_flow[ix])
				
				res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]

				stats =  {res_stats[0]: slopes,
				                  res_stats[1]: intercepts, 
				                  res_stats[2]: means,
				                  res_stats[3]: stds}
				stats = pd.DataFrame(stats)
				if res_id == 'SHA':
					stats_file = stats
				else:
					stats_file[res_stats] = stats
				stats = stats.values.T
				for i,s in enumerate(stats):

					yearly_stats = stats[i]
					if i == 0:
						res_slope = np.append(res_slope,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_slope = np.append(res_slope,[yearly_stats[364]])
					elif i == 1:
						res_intercept = np.append(res_intercept,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_intercept = np.append(res_intercept,[yearly_stats[364]])
					elif i == 2:
						res_mean = np.append(res_mean,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_mean = np.append(res_mean,[yearly_stats[364]])
					elif i == 3:
						res_std = np.append(res_std,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_std = np.append(res_std,[yearly_stats[364]])

			r[res_stats[0]] = pd.Series(res_slope,index=r.index)
			r[res_stats[1]] = pd.Series(res_intercept,index=r.index)
			r[res_stats[2]] = pd.Series(res_mean,index=r.index)
			r[res_stats[3]] = pd.Series(res_std,index=r.index)

			r.rename(columns = {'cum_flow_to_date':'%s_cum_flow_to_date'%res_id}, inplace=True)
			r.rename(columns = {'remaining_flow':'%s_remaining_flow'%res_id}, inplace=True)
			r.rename(columns = {'snowpack':'%s_snowpack'%res_id}, inplace=True)
			r.drop(['inf','WY','DOWY'], axis=1, inplace=True)

			df = pd.concat([df, r], axis=1, join_axes=[df.index])
		df = pd.concat([df,WYI_stats], axis=1, join_axes=[df.index])


	elif window_type == 'expanding':
		decade_thresh = pd.date_range('1951-10-01','2099-10-01',freq ='AS-OCT')

		WYI_mov_stats = get_forecast_WYI_stats(df.truncate(before = decade_thresh[2], after=decade_thresh[50]),index_exceedance_sac)
		WYI_sim = get_projection_forecast_WYI(df.truncate(before=decade_thresh[0], after=decade_thresh[50]),WYI_mov_stats,index_exceedance_sac)
		WYI_stats = WYI_sim[1]
		WYI_sim = WYI_sim[0]

		for i in range(50,148):
			WYI_mov_stats = get_forecast_WYI_stats(df.truncate(before = decade_thresh[2], after=decade_thresh[i]),index_exceedance_sac)
			WYI_dec = get_projection_forecast_WYI(df.truncate(before=decade_thresh[i], after=decade_thresh[i+1]),WYI_mov_stats,index_exceedance_sac)
			WYI_sim = pd.concat([WYI_sim,WYI_dec[0]])
			WYI_stats = pd.concat([WYI_stats, WYI_dec[1]])
		snow_sites = ['BND_swe', 'ORO_swe','FOL_swe']
		res_ids = ['SHA','ORO','FOL']
		SHA_inf = (df['SHA_in_tr'].to_frame(name='inf') * cfs_tafd)  
		ORO_inf = (df['ORO_in_tr'].to_frame(name='inf') * cfs_tafd)
		FOL_inf = (df['FOL_in_tr'].to_frame(name='inf') * cfs_tafd)
		res_frames = [SHA_inf,ORO_inf,FOL_inf]

		for r, swe, res_id in zip(res_frames, snow_sites, res_ids):
			rstats = pd.DataFrame()
			r['WY'] = df.WY
			r['DOWY'] = df.DOWY
			r['snowpack'] = df[swe]
			r = r[r.WY != r.WY[-1]]
			month = r.index.month

			r['cum_flow_to_date'] = (r.groupby('WY').inf
			                           .apply(flow_to_date))
			r['remaining_flow'] = (r.groupby('WY').inf
			                            .apply(rem_flow))
			r.cum_flow_to_date.fillna(method='ffill', inplace=True)
			r = r.replace([np.inf, -np.inf], np.nan)
			r = r.fillna(method = 'ffill')
			rstart = r.truncate(before = decade_thresh[2], after=decade_thresh[50])
			rstart = rstart[(rstart.index < decade_thresh[50])]
			slopes = np.zeros(365)
			intercepts = np.zeros(365)
			means = np.zeros(365)
			stds = np.zeros(365)
			snowpack = rstart.snowpack + 0.0001
			remaining_flow = rstart.remaining_flow 
			DOWY = rstart.DOWY
			for d in range(1,365):
				ix = (DOWY == d)
				coeffs = np.polyfit(snowpack[ix],remaining_flow[ix],1)
				slopes[d] = coeffs[0]
				intercepts[d] = coeffs[1]
				means[d] = np.mean(remaining_flow[ix]) 
				stds[d] = np.std(remaining_flow[ix])
			

			res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]
			stats =  {res_stats[0]: slopes,
			                  res_stats[1]: intercepts, 
			                  res_stats[2]: means,
			                  res_stats[3]: stds}
			stats = pd.DataFrame(stats)
			if res_id == 'SHA':
				stats_file = stats
			else:
				stats_file[res_stats] = stats
			stats = stats.values.T
			rstatsind = r.truncate(before = decade_thresh[0], after=decade_thresh[50])
			rstatsind = rstatsind[(rstatsind.index < decade_thresh[50])].index

			for i,s in enumerate(stats):
				yearly_stats = stats[i]
				if i == 0:
					res_slope = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_slope = np.append(res_slope,np.tile(yearly_stats, 4))
						res_slope = np.append(res_slope,[yearly_stats[364]]) #leap year
					res_slope = np.append(res_slope,np.tile(yearly_stats, 1))
				
				elif i == 1:
					res_intercept = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_intercept = np.append(res_intercept,np.tile(yearly_stats, 4))
						res_intercept = np.append(res_intercept,[yearly_stats[364]]) #leap year
					res_intercept = np.append(res_intercept,np.tile(yearly_stats, 1))
				
				elif i == 2:
					res_mean = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_mean = np.append(res_mean,np.tile(yearly_stats, 4))
						res_mean = np.append(res_mean,[yearly_stats[364]]) #leap year
					res_mean = np.append(res_mean,np.tile(yearly_stats, 1))

				elif i == 3:
					res_std = np.append(yearly_stats,[yearly_stats[364]]) #1952 WY
					for y in range(12): # 1951  WYs
						res_std = np.append(res_std,np.tile(yearly_stats, 4))
						res_std = np.append(res_std,[yearly_stats[364]]) #leap year
					res_std = np.append(res_std,np.tile(yearly_stats, 1))
				
				# rstats[res_stats[i]] = pd.Series(v,index=rstatsind)
				
			for y in range(50,148):
				rnext = r.truncate(before = decade_thresh[2], after=decade_thresh[y])
				rnext = rnext[(rnext.index < decade_thresh[y])]
				rnow_ind = r.truncate(before=decade_thresh[y], after=decade_thresh[y+1])
				rnow_ind = rnow_ind[(rnow_ind.index < decade_thresh[y+1])].index
				slopes = np.zeros(365)
				intercepts = np.zeros(365)
				means = np.zeros(365)
				stds = np.zeros(365)
				snowpack = rnext.snowpack + 0.0001
				remaining_flow = rnext.remaining_flow+0.0001
				DOWY = rnext.DOWY
				for d in range(1,365):
					ix = (DOWY == d)
					coeffs = np.polyfit(snowpack[ix],remaining_flow[ix],1)
					slopes[d] = coeffs[0]
					intercepts[d] = coeffs[1]
					means[d] = np.mean(remaining_flow[ix]) 
					stds[d] = np.std(remaining_flow[ix])
				
				res_stats =['%s_slope'%res_id,'%s_intercept'%res_id,'%s_mean'%res_id,'%s_std'%res_id]

				stats =  {res_stats[0]: slopes,
				                  res_stats[1]: intercepts, 
				                  res_stats[2]: means,
				                  res_stats[3]: stds}
				stats = pd.DataFrame(stats)
				if res_id == 'SHA':
					stats_file = stats
				else:
					stats_file[res_stats] = stats
				stats = stats.values.T
				for i,s in enumerate(stats):

					yearly_stats = stats[i]
					if i == 0:
						res_slope = np.append(res_slope,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_slope = np.append(res_slope,[yearly_stats[364]])
					elif i == 1:
						res_intercept = np.append(res_intercept,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_intercept = np.append(res_intercept,[yearly_stats[364]])
					elif i == 2:
						res_mean = np.append(res_mean,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_mean = np.append(res_mean,[yearly_stats[364]])
					elif i == 3:
						res_std = np.append(res_std,np.tile(yearly_stats, 1)) #1951 
						if decade_thresh[y+1].is_leap_year:
							res_std = np.append(res_std,[yearly_stats[364]])

			r[res_stats[0]] = pd.Series(res_slope,index=r.index)
			r[res_stats[1]] = pd.Series(res_intercept,index=r.index)
			r[res_stats[2]] = pd.Series(res_mean,index=r.index)
			r[res_stats[3]] = pd.Series(res_std,index=r.index)
			
			r.rename(columns = {'cum_flow_to_date':'%s_cum_flow_to_date'%res_id}, inplace=True)
			r.rename(columns = {'remaining_flow':'%s_remaining_flow'%res_id}, inplace=True)
			r.rename(columns = {'snowpack':'%s_snowpack'%res_id}, inplace=True)
			r.drop(['inf','WY','DOWY'], axis=1, inplace=True)
			df = pd.concat([df, r], axis=1, join_axes=[df.index])
		df = pd.concat([df,WYI_stats], axis=1, join_axes=[df.index])









	df['WYI_sim'] = WYI_sim
	df.WYI_sim = df.WYI_sim.fillna(method = 'bfill')
	df.loc[df['WYI_sim'].isnull(),'WYI_sim'] = df['SR_WYI']

	df['WYT_sim'] = df.WYI_sim.apply(WYI_to_WYT,
                               thresholds=[9.2, 7.8, 6.5, 5.4, 0.0], 
                               values=['W', 'AN', 'BN', 'D', 'C'])
	df = df[(df.index < '2099-10-01')]
	return df

def get_forecast_WYI_stats(df, index_exceedance_sac): #now determining forecasting regression coefficients based off perfect foresight
	flow_sites = ['BND_fnf', 'ORO_fnf', 'YRS_fnf', 'FOL_fnf']
	snow_sites = ['BND_swe', 'ORO_swe', 'YRS_swe', 'FOL_swe']

	Qm = (df[flow_sites].sum(axis=1).resample('MS').sum().to_frame(name='flow') * cfsd_mafd)

	Qm['WY'] = df.WY.resample('MS').first() # need this for grouping
	Qm['snow'] = df[snow_sites].sum(axis=1).resample('MS').first()
	Qm = Qm[Qm.WY != Qm.WY[-1]] 
	Qm['octmar_cumulative'] = (Qm.groupby('WY').flow
	                          .apply(octmar_cumulative)
	                          .reset_index(level=0)
	                          .drop('WY', axis=1)) 

	Qm['aprjul_cumulative'] = (Qm.groupby('WY').flow
	                           .apply(aprjul_cumulative))

	Qm['octmar_flow_to_date'] = (Qm.groupby('WY').flow
	                           .apply(octmar_flow_to_date))

	Qm['aprjul_flow_to_date'] = (Qm.groupby('WY').flow
	                           .apply(aprjul_flow_to_date)
	                           .reset_index(level=0).drop('WY', axis=1)) 

#obtaining regression coefficients for WYI forecasting
	aprjul_slopes = np.zeros(12)
	aprjul_intercepts = np.zeros(12)
	aprjul_means = np.zeros(12)
	aprjul_stds = np.zeros(12)
	octmar_means = np.zeros(12)
	octmar_stds = np.zeros(12)
	octmar_slopes = np.zeros(12)
	octmar_intercepts = np.zeros(12)
	month = Qm.index.month
	snow_ind = Qm.snow
	aprjul_cumulative_ind = Qm.aprjul_cumulative
	for m in range(1,13):
		oct_index = (month == 10)
		ix = (month == m)
		coeffs = np.polyfit(snow_ind[ix],aprjul_cumulative_ind[ix],1)
		aprjul_slopes[m-1] = coeffs[0]
		aprjul_intercepts[m-1] = coeffs[1]
		aprjul_means[m-1] = np.mean(aprjul_cumulative_ind[ix]) 
		aprjul_stds[m-1] =  np.std(aprjul_cumulative_ind[ix])
	octmar_flow_to_date_ind = Qm.octmar_flow_to_date	
	octmar_cumulative_ind = Qm.octmar_cumulative
	for m in list(range(10,13))+list(range(1,3)):  
		ix = (month == m)
		flow_coeffs = np.polyfit(octmar_flow_to_date_ind[oct_index], octmar_cumulative_ind[ix], 1)
		octmar_slopes[m-1] = flow_coeffs[0]
		octmar_intercepts[m-1] = flow_coeffs[1]
		octmar_means[m-1] = np.mean(octmar_cumulative_ind[ix])

		octmar_stds[m-1] = np.std(octmar_cumulative_ind[ix])
	stats =  {'aprjul_slope': aprjul_slopes,
                    'aprjul_intercept': aprjul_intercepts, 
                    'aprjul_mean': aprjul_means,
                    'aprjul_std':  aprjul_stds, 
                    'octmar_mean': octmar_means,
                    'octmar_std':  octmar_stds, 
                    'octmar_intercept': octmar_intercepts, 
                    'octmar_slope': octmar_slopes}

	stats = pd.DataFrame(stats, columns = ['aprjul_slope',
                    'aprjul_intercept', 
                    'aprjul_mean',
                    'aprjul_std', 
                    'octmar_mean',
                    'octmar_std', 
                    'octmar_intercept', 
                    'octmar_slope'])
	return(stats)
