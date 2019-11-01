import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.dates as mdates

def init_plotting():
  sns.set_style("darkgrid", {"axes.facecolor": "0.9"}) 
  # sns.palplot(sns.color_palette("Set2", 8))
  # current_palette = sns.color_palette()
  # sns.palplot(current_palette)

  # plt.rcParams['figure.figsize'] = (15, 8)

  plt.rcParams['figure.figsize'] = (10, 5)
  # plt.rcParams['axes.prop_cycle'] = cycler(color='#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8')

  plt.rcParams['font.size'] = 13
  plt.rcParams['lines.linewidth'] = 1.5
  plt.rcParams['lines.linestyle'] = '-'

  # plt.rcParams['font.family'] = 'OfficinaSanITCBoo'
  plt.rcParams['font.weight'] = 'bold'
  plt.rcParams['axes.labelsize'] = 1.1*plt.rcParams['font.size']
  plt.rcParams['axes.titlesize'] = 1.1*plt.rcParams['font.size']
  plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
  plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
  plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
  # plt.style.use('ggplot')
init_plotting()
fig, (ax0,ax1,ax2) = plt.subplots(3,1,sharey = True)
model = 'ccsm'
output = 'SHA_storage'
output_name = 'Shasta storage'
dfhist = pd.read_csv('%s19952005/%s19952005-results.csv'%(model,model),index_col = 0,parse_dates = True)
dfmid = pd.read_csv('%s20402050/%s20402050-results.csv'%(model,model),index_col = 0,parse_dates = True)
dfend = pd.read_csv('%s20902100/%s20902100-results.csv'%(model,model),index_col = 0,parse_dates = True)

ax0.plot(dfhist[output])
ax1.plot(dfmid[output])
ax2.plot(dfend[output])
ax0.xaxis.set_major_locator(mdates.YearLocator())
ax1.xaxis.set_major_locator(mdates.YearLocator())
ax2.xaxis.set_major_locator(mdates.YearLocator())
ax1.set_ylabel('TAF')
fig.suptitle('%s %s'%(model,output_name))
plt.tight_layout()
plt.show()