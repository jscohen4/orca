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
  res = res.resample(freq).mean()
  obs = obs.resample(freq).mean()

  fig = plt.figure(figsize=(12, 3)) 
  gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

  ax0 = plt.subplot(gs[0])
  res.plot(ax=ax0, color='indianred')
  obs.plot(ax=ax0, color='k')
  ax0.set_title('%s, %s' % (res.name, obs.name), loc='left') # family='sans-serif'
  ax0.legend(['Simulated', 'Observed'], ncol=3)

  ax1 = plt.subplot(gs[1])
  r = np.corrcoef(obs.values,res.values)[0,1]
  ax1.scatter(obs.values, res.values, s=3, c='steelblue', edgecolor='none', alpha=0.7)
  ax1.set_ylabel('Simulated')
  ax1.set_xlabel('Observed')
  ax1.annotate('$R^2 = %f$' % r**2, xy=(0,0), color='0.3')
  ax1.set_xlim([0.0, ax1.get_xlim()[1]])
  ax1.set_ylim([0.0, ax1.get_ylim()[1]])

  plt.tight_layout()
  # plt.show()
