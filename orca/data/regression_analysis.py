import numpy as np 
import scipy.stats as sp
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from util import *
import matplotlib.pyplot as plt
from itertools import * 
import matplotlib.mlab as mlab
sns.set_style('whitegrid')

def plot_regressions():
	df = pd.read_csv('orca-data-cdf-generated.csv', index_col=0, parse_dates=True)
	#f, (ax1,a2,ax3,ax4,ax5,ax6,ax7,ax8)
	fig, axes = plt.subplots(nrows=4, ncols=2)

	(df['SHA_octmar_slope']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[0,0])
	(df['ORO_octmar_slope']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[0,0])
	(df['FOL_octmar_slope']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[0,0])

	(df['SHA_aprjul_slope']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[0,1])
	(df['ORO_aprjul_slope']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[0,1])
	(df['FOL_aprjul_slope']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[0,1])

	(df['SHA_octmar_intercept']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[1,0])
	(df['ORO_octmar_intercept']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[1,0])
	(df['FOL_octmar_intercept']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[1,0])

	(df['SHA_aprjul_intercept']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[1,1])
	(df['ORO_aprjul_intercept']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[1,1])
	(df['FOL_aprjul_intercept']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[1,1])

	(df['SHA_octmar_mean']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[2,0])
	(df['ORO_octmar_mean']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[2,0])
	(df['FOL_octmar_mean']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[2,0])

	(df['SHA_aprjul_mean']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True,ax=axes[2,1])
	(df['ORO_aprjul_mean']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True,ax=axes[2,1])
	(df['FOL_aprjul_mean']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True,ax=axes[2,1])

	(df['SHA_octmar_std']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[3,0])
	(df['ORO_octmar_std']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[3,0])
	(df['FOL_octmar_std']).plot(xlim = ('1996-10-01','1997-03-30'), legend = True, ax=axes[3,0])

	(df['SHA_aprjul_std']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[3,1])
	(df['ORO_aprjul_std']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[3,1])
	(df['FOL_aprjul_std']).plot(xlim = ('1997-04-01','1997-07-31'), legend = True, ax=axes[3,1])
	plt.show()

	#october.plot()
	#ax.set_xlim(pd.Timestamp('1996-10-01'), pd.Timestamp('1996-03-30'))

def plot_new_regressions():
	df = pd.read_csv('orca-data-cdf-generated-new.csv', index_col=0, parse_dates=True)
	df2 = pd.read_csv('orca-data-cdf-generated-monthly3.csv', index_col=0, parse_dates=True)
	#f, (ax1,a2,ax3,ax4,ax5,ax6,ax7,ax8)
	fig, axes = plt.subplots(nrows=2, ncols=2)

	(df['SHA_slope']).plot(xlim = ('1997-10-02','1998-09-30'), legend = True, ax=axes[0,0], color = 'r')
	(df['ORO_slope']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,0], color = 'b')
	(df['FOL_slope']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,0], color = 'g')
	(df2['SHA_slope']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,0], color = 'r', ls = ':')
	(df2['ORO_slope']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,0], color = 'b', ls = ':')
	(df2['FOL_slope']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,0], color = 'g', ls = ':')

	(df['SHA_intercept']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,1], color = 'r')
	(df['ORO_intercept']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,1], color = 'b')
	(df['FOL_intercept']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,1], color = 'g')
	(df2['SHA_intercept']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,1], color = 'r', ls = ':')
	(df2['ORO_intercept']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,1], color = 'b', ls = ':')
	(df2['FOL_intercept']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[0,1], color = 'g', ls = ':')

	(df['SHA_mean']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True,ax=axes[1,0], color = 'r')
	(df['ORO_mean']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True,ax=axes[1,0], color = 'b')
	(df['FOL_mean']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True,ax=axes[1,0], color = 'g')
	(df2['SHA_mean']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,0], color = 'r', ls = ':')
	(df2['ORO_mean']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,0], color = 'b', ls = ':')
	(df2['FOL_mean']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,0], color = 'g', ls = ':')

	(df['SHA_std']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,1], color = 'r')
	(df['ORO_std']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,1], color = 'b')
	(df['FOL_std']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,1], color = 'g')
	(df2['SHA_std']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,1], color = 'r', ls = ':')
	(df2['ORO_std']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,1], color = 'b', ls = ':')
	(df2['FOL_std']).plot(xlim = ('1997-10-01','1998-09-30'), legend = True, ax=axes[1,1], color = 'g', ls = ':')
	plt.show()

plot_new_regressions()