
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
  res.plot(ax=ax0, color='indianred')
  obs.plot(ax=ax0, color='k')
  ax0.set_title('%s, %s' % (res.name, obs.name), loc='left') # family='sans-serif'
  ax0.legend(['simulated with WRF', 'simulated with CDEC'], ncol=3)

  ax1 = plt.subplot(gs[1])
  r = np.corrcoef(obs.values,res.values)[0,1]
  ax1.scatter(obs.values, res.values, s=3, c='steelblue', edgecolor='none', alpha=0.7)
  ax1.set_ylabel('simulated with WRF')
  ax1.set_xlabel('simulated with CDEC')
  ax1.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax1.set_xlim([0.0, ax1.get_xlim()[1]])
  ax1.set_ylim([0.0, ax1.get_ylim()[1]])

  plt.tight_layout()
  # plt.show()
def compare_WRF_obs(res,obs,freq='D'):
  # input two pandas series and a frequency
  sns.set_style('whitegrid')
  plt.rcParams['figure.figsize'] = (12, 8)
  res = res.resample(freq).sum()
  obs = obs.resample(freq).sum()

  fig = plt.figure(figsize=(12, 3)) 
  gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  res.plot(ax=ax0, color='indianred')
  obs.plot(ax=ax0, color='k')
  ax0.set_title('%s, %s' % (res.name, obs.name), loc='left') # family='sans-serif'
  ax0.legend(['simulated with WRF', 'observed'], ncol=3)

  ax1 = plt.subplot(gs[1])
  r = np.corrcoef(obs.values,res.values)[0,1]
  ax1.scatter(obs.values, res.values, s=3, c='steelblue', edgecolor='none', alpha=0.7)
  ax1.set_ylabel('simulated with WRF')
  ax1.set_xlabel('observed')
  ax1.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax1.set_xlim([0.0, ax1.get_xlim()[1]])
  ax1.set_ylim([0.0, ax1.get_ylim()[1]])

  plt.tight_layout()

def compare_cdec_obs(res,obs,freq='D'):
  # input two pandas series and a frequency
  sns.set_style('whitegrid')
  plt.rcParams['figure.figsize'] = (12, 8)
  res = res.resample(freq).sum()
  obs = obs.resample(freq).sum()

  fig = plt.figure(figsize=(12, 3)) 
  gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  res.plot(ax=ax0, color='indianred')
  obs.plot(ax=ax0, color='k')
  ax0.set_title('%s, %s' % (res.name, obs.name), loc='left') # family='sans-serif'
  ax0.legend(['simulated with CDEC', 'observed'], ncol=3)

  ax1 = plt.subplot(gs[1])
  r = np.corrcoef(obs.values,res.values)[0,1]
  ax1.scatter(obs.values, res.values, s=3, c='steelblue', edgecolor='none', alpha=0.7)
  ax1.set_ylabel('simulated with CDEC')
  ax1.set_xlabel('observed')
  ax1.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax1.set_xlim([0.0, ax1.get_xlim()[1]])
  ax1.set_ylim([0.0, ax1.get_ylim()[1]])

  plt.tight_layout()

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


def compare_both(WRF,CDEC,OBS,freq='D'):
  # input two pandas series and a frequency
  sns.set_style('whitegrid')
  plt.rcParams['figure.figsize'] = (12, 8)
  WRF = WRF.resample(freq).sum()
  CDEC = CDEC.resample(freq).sum()
  OBS = OBS.resample(freq).sum()
  fig = plt.figure(figsize=(15, 3)) 
  gs = gridspec.GridSpec(1, 3, width_ratios=[3, 1, 1]) 

  ax0 = plt.subplot(gs[0])
  WRF.plot(ax=ax0, color='b')
  CDEC.plot(ax=ax0, color='indianred')
  OBS.plot(ax=ax0, color='k')
  ax0.set_title('%s' % (WRF.name), loc='left') # family='sans-serif'
  if WRF.name == 'Combined_pump':
    ax0.set_ylabel('TAF/%s'%freq)
  elif WRF.name in ['SHA_storage','ORO_storage','FOL_storage']:
    ax0.set_ylabel('TAF')

  ax0.legend(['simulated with WRF', 'simulated with CDEC','observed'], ncol=3)

  ax1 = plt.subplot(gs[1])
  r = np.corrcoef(OBS.values,WRF.values)[0,1]
  ax1.scatter(OBS.values, WRF.values, s=3, c='b', edgecolor='none', alpha=0.7)
  ax1.set_ylabel('Simulated with WRF')
  ax1.set_xlabel('Observed')
  ax1.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax1.set_xlim([0.0, ax1.get_xlim()[1]])
  ax1.set_ylim([0.0, ax1.get_ylim()[1]])

  ax2 = plt.subplot(gs[2])
  r = np.corrcoef(OBS.values,CDEC.values)[0,1]
  ax2.scatter(OBS.values, CDEC.values, s=3, c='indianred', edgecolor='none', alpha=0.7)
  ax2.set_ylabel('Simulated with CDEC')
  ax2.set_xlabel('Observed')
  ax2.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax2.set_xlim([0.0, ax2.get_xlim()[1]])
  ax2.set_ylim([0.0, ax2.get_ylim()[1]])

  plt.tight_layout()
