import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from .util import *
import seaborn as sns
def compare(res,obs,freq='D'):
  # input two pandas series and a frequency
  sns.set_style('whitegrid')
  plt.rcParams['figure.figsize'] = (12, 8)
  res = res.resample(freq).sum()
  obs = obs.resample(freq).sum()

  fig = plt.figure(figsize=(12, 3)) 
  gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  res = res
  obs = obs
  res.plot(ax=ax0, color='indianred')
  obs.plot(ax=ax0, color='k')
  ax0.set_title('%s, %s' % (res.name, obs.name), loc='left') # family='sans-serif'
  ax0.legend(['Simulated', 'Observed'], ncol=3)

  ax1 = plt.subplot(gs[1])
  NSE = 1 - sum((res.values-obs.values)**2)/sum((obs.values-np.mean(obs.values))**2)
  # r = np.corrcoef(obs.values,res.values)[0,1]
  ax1.scatter(obs.values, res.values, s=3, c='steelblue', edgecolor='none', alpha=0.7)
  ax1.set_ylabel('Simulated')
  ax1.set_xlabel('Observed')
  # ax1.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax1.annotate('$NSE = %f$'%NSE, xy=(0,0), color='0.3')

  ax1.set_xlim([0.0, ax1.get_xlim()[1]])
  ax1.set_ylim([0.0, ax1.get_ylim()[1]])

  # ax0.set_ylabel('TAF/day')
  ax0.set_xlabel('Date')
  # ax1.set_ylabel('TAF/day')
  # ax1.set_xlabel('TAF/day')

  plt.tight_layout()
  # plt.show()

def plotting(res,freq='D'):
  # input two pandas series and a frequency
  sns.set_style('whitegrid')
  plt.rcParams['figure.figsize'] = (12, 8)
  res = res.resample(freq).sum()

  fig = plt.figure(figsize=(12, 3)) 
  gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  res.plot(ax=ax0, color='indianred')
  ax0.set_title('%s' % (res.name), loc='left') # family='sans-serif'
  ax0.legend(['Simulated'], ncol=3)
  plt.tight_layout()

def Rsquares(sim,obs,R2sfile):
  text_file = open(R2sfile, 'w')
  r2s = []
  r2s.append(['Timestep','Daily','Weekly','Monthly','Water Year'])
  for s,o in zip(sim,obs):
    r2s_point = []
    r2s_point.append('%s ' %s.name)
    for i,f in enumerate(['D','W','M','AS-OCT']):
      sim = s.resample(f).sum()
      obs = o.resample(f).sum()
      r = np.corrcoef(obs.values,sim.values)[0,1]
      r2s_point.append('%.5s' %r**2)
    r2s_point.append('\n')
    r2s.append(r2s_point)
  for line in r2s: 
    text_file.write('{:<14} {:<12} {:<12} {:<12} {:<12}'.format(*line))
    text_file.write(' \n')
      #text_file.write('\n')
  text_file.close()

def NSE(sim,obs,NSEfile):
  text_file = open(NSEfile, 'w')
  NSEs = []
  NSEs.append(['Timestep','Daily','Weekly','Monthly','Water Year'])
  for s,o in zip(sim,obs):
    NSEs_point = []
    NSEs_point.append('%s ' %s.name)
    for i,f in enumerate(['D','W','M','AS-OCT']):
      sim = s.resample(f).sum()
      obs = o.resample(f).sum()
      NSE1 = 1 - sum((sim-obs)**2)/sum((obs-np.mean(obs))**2)
      NSEs_point.append('%.5s' %NSE1)
    NSEs_point.append('\n')
    NSEs.append(NSEs_point)
  for line in NSEs: 
    text_file.write('{:<14} {:<12} {:<12} {:<12} {:<12}'.format(*line))
    text_file.write(' \n')
      #text_file.write('\n')
  text_file.close()

def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    whereever the observed data contains nan
    
    this is used by all other functions, otherwise they will produce nan as 
    output
    """
    if np.sum(~np.isnan(s*o))>=1:
        data = np.array([s.flatten(),o.flatten()])
        data = np.transpose(data)
        data = data[~np.isnan(data).any(1)]
        s = data[:,0]
        o = data[:,1]
    return s, o

def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    s,o = filter_nan(s,o)
    return 1 - sum((s-o)**2)/sum((o-np.mean(o))**2)
