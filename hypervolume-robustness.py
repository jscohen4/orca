import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.stats import norm,skew
import seaborn as sns
import inspyred.ec.analysis as pyred
import pickle
import math

def init_plotting():
  sns.set_style("darkgrid", {"axes.facecolor": "0.9"}) 
  # plt.rcParams['figure.figsize'] = (15, 8)
  plt.rcParams['figure.figsize'] = (10, 10)

  plt.rcParams['font.size'] = 13
  plt.rcParams['lines.linewidth'] = 1.5
  plt.rcParams['lines.linestyle'] = '-'

  plt.rcParams['font.family'] = 'Sans.Sarif' 
  plt.rcParams['font.weight'] = 'bold'
  plt.rcParams['axes.labelsize'] = 1*plt.rcParams['font.size']
  plt.rcParams['axes.titlesize'] = 1.1*plt.rcParams['font.size']
  plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
  plt.rcParams['xtick.labelsize'] = 1.2*plt.rcParams['font.size']
  plt.rcParams['ytick.labelsize'] = 1.2*plt.rcParams['font.size']
init_plotting()
fig,(ax0,ax1,ax2) = plt.subplots(3,1,sharex = True)
def hypervolume(pareto_set, reference_point=None):
    """Calculates the hypervolume by slicing objectives (HSO).
    
    This function calculates the hypervolume (or S-measure) of a nondominated
    set using the Hypervolume by Slicing Objectives (HSO) procedure of `While, et al. 
    (IEEE CEC 2005) <http://www.lania.mx/~ccoello/EMOO/while05a.pdf.gz>`_.
    The *pareto_set* should be a list of lists of objective values.
    The *reference_point* may be specified or it may be left as the default 
    value of None. In that case, the reference point is calculated to be the
    maximum value in the set for all objectives (the ideal point). This function 
    assumes that objectives are to be maximized.
    
    Arguments:
    
    - *pareto_set* -- the list or lists of objective values comprising the Pareto front
    - *reference_point* -- the reference point to be used (default None)
    
    """
    def dominates(p, q, k=None):
        if k is None:
            k = len(p)
        d = True
        while d and k < len(p):
            d = not (q[k] > p[k])
            k += 1
        return d
        
    def insert(p, k, pl):
        ql = []
        while pl and pl[0][k] > p[k]:
            ql.append(pl[0])
            pl = pl[1:]
        ql.append(p)
        while pl:
            if not dominates(p, pl[0], k):
                ql.append(pl[0])
            pl = pl[1:]
        return ql

    def slice(pl, k, ref):
        p = pl[0]
        pl = pl[1:]
        ql = []
        s = []
        while pl:
            ql = insert(p, k + 1, ql)
            p_prime = pl[0]
            s.append((math.fabs(p[k] - p_prime[k]), ql))
            p = p_prime
            pl = pl[1:]
        ql = insert(p, k + 1, ql)
        s.append((math.fabs(p[k] - ref[k]), ql))
        return s

    ps = pareto_set
    ref = reference_point
    if len(ps) == 0:
    	return 0
    n = min([len(p) for p in ps])

    if ref is None:
        ref = [max(ps, key=lambda x: x[o])[o] for o in range(n)]
    pl = ps[:]
    pl.sort(key=lambda x: x[0], reverse=True)
    s = [(1, pl)]
    for k in range(n - 1):
        s_prime = []
        for x, ql in s:
            for x_prime, ql_prime in slice(ql, k, ref):
                s_prime.append((x * x_prime, ql_prime))
        s = s_prime

    vol = 0
    for x, ql in s:
        vol = vol + x * math.fabs(ql[0][n - 1] - ref[n - 1])
    return vol


baseline = pd.read_csv('baseline-objectives.csv',index_col = 0)
with open('scenario_names_all.txt') as f:
  scenarios = f.read().splitlines()
with open('scenario-groupings/HV-FNF-SWE/high_potential.txt') as f:
    scenarios_HP = f.read().splitlines()
with open('scenario-groupings/HV-FNF-SWE/upper.txt') as f:
    scenarios_upper = f.read().splitlines()
with open('scenario-groupings/HV-FNF-SWE/lower.txt') as f:
    scenarios_lower = f.read().splitlines()

hv_perf = pd.read_csv('hypervolumes-perf-3.csv')
yearstart = 2069
yearend = 2098
spill = pd.read_csv('baseline-results/baseline_spill.csv', parse_dates = True,index_col = 0)
carryover = pd.read_csv('baseline-results/baseline_carryover.csv', parse_dates = True,index_col = 0)
delta = pd.read_csv('baseline-results/baseline_delta.csv', parse_dates = True,index_col = 0)
pumping = pd.read_csv('baseline-results/baseline_pumping.csv', parse_dates = True,index_col = 0)
hydro = pd.read_csv('baseline-results/baseline_energy.csv', parse_dates = True,index_col = 0)
spillbase = spill[spill.index.year >= yearstart].cumsum()
carryoverbase = carryover[carryover.index.year >= yearstart].cumsum()
deltabase = delta[delta.index.year >= yearstart].cumsum()
pumpingbase = pumping[pumping.index.year >= yearstart].cumsum()
hydrobase = hydro[hydro.index.year >= yearstart].cumsum()
# hyperarray = np.zeros(len(scenarios_train))
dfHPset = pd.DataFrame(index = range(len(scenarios_HP)))
dfUset = pd.DataFrame(index = range(len(scenarios_upper)))
dfLset = pd.DataFrame(index = range(len(scenarios_lower)))




HV_perf = hv_perf.HV.values
training_sets = ['high-potential','upper','lower','high-potential-lower','high-potential-upper','upper-lower']
for train_set in training_sets:
    hyperarray_HP = []
    hyperarray_upper = []
    hyperarray_lower = []
    hyperarray_perf_HP = []
    hyperarray_perf_upper = []
    hyperarray_perf_lower = []
    for j,sc in enumerate(scenarios):
        print(j)
        if train_set in ['high-potential','upper','lower']:
            dfsc = pd.read_csv('training-outputs/individual-groups/%s/scenario%s.csv'%(train_set,j),index_col = 0)
        elif train_set in ['high-potential-lower','high-potential-upper','upper-lower']:
            dfsc = pd.read_csv('training-outputs/mixed-groups/%s/scenario%s.csv'%(train_set,j),index_col = 0)
       
        obj_raw = dfsc.values
        obj0 = obj_raw[:,0]
        obj1 = obj_raw[:,1]
        obj2 = obj_raw[:,2]
        obj3 = obj_raw[:,3]
        obj4 = obj_raw[:,4]
        objnorm = np.zeros([len(obj0),4])
        # print(obj3)

        for i,obj in enumerate(obj0):
          objs = (obj - obj0.min(axis=0)) / (obj0.max(axis=0) - obj0.min(axis=0))
          objnorm[i][0] = objs

        for i,obj in enumerate(obj1):
          objs = (obj - obj1.min(axis=0)) / (obj1.max(axis=0) - obj1.min(axis=0))
          objnorm[i][1] = objs

        for i,obj in enumerate(obj2):
          objs = (obj - obj2.min(axis=0)) / (obj2.max(axis=0) - obj2.min(axis=0))
          objnorm[i][2] = objs

        for i,obj in enumerate(obj3):
          objs = (obj - obj3.min(axis=0)) / (obj3.max(axis=0) - obj3.min(axis=0))
          objnorm[i][3] = objs

    # for i,obj in enumerate(obj4):
    #   objs = (obj - obj4.min(axis=0)) / (obj4.max(axis=0) - obj4.min(axis=0))
    #   objnorm[i][4] = objs

        baseline0 = spillbase[spillbase.index == '%s-10-01'%yearend][sc].values[0]*2.5
        baseline0 = (baseline0 - obj0.max(axis=0)) / (obj0.min(axis=0) - obj0.max(axis=0))

        baseline1 = carryoverbase[carryoverbase.index == '%s-10-01'%yearend][sc].values[0]*0.9
        baseline1 = (baseline1 - obj1.min(axis=0)) / (obj1.max(axis=0) - obj1.min(axis=0))

        baseline2 = deltabase[deltabase.index == '%s-10-01'%yearend][sc].values[0]*0.9
        baseline2 = (baseline2 - obj2.min(axis=0)) / (obj2.max(axis=0) - obj2.min(axis=0))

        baseline3 = pumpingbase[pumpingbase.index == '%s-10-01'%yearend][sc].values[0]*0.9
        # print(baseline3)
        baseline3 = (baseline3 - obj3.min(axis=0)) / (obj3.max(axis=0) - obj3.min(axis=0))

        objnorm = -objnorm
        ref_pt = [-abs(baseline0),-abs(baseline1),-abs(baseline2),-abs(baseline3)]
        oblnormfilter = []
        ref_ptfilter = np.array(ref_pt)
        count = 0
        for point in objnorm: 
            compare = np.greater(point,ref_ptfilter)
            if compare.all() == True:
                count += 1
                oblnormfilter.append(point.tolist())

        objnorm  = objnorm.tolist()
        hypev = hypervolume(oblnormfilter,reference_point = ref_pt)	

        if sc in scenarios_HP:
            hyperarray_HP.append(hypev)
            hyperarray_perf_HP.append(HV_perf[j])
        elif sc in scenarios_upper:
            hyperarray_upper.append(hypev)
            hyperarray_perf_upper.append(HV_perf[j])
        elif sc in scenarios_lower:
            hyperarray_lower.append(hypev)
            hyperarray_perf_lower.append(HV_perf[j])
    dfHPset[train_set] = np.array(hyperarray_HP)/np.array(hyperarray_perf_HP)
    dfUset[train_set] =np.array(hyperarray_upper)/np.array(hyperarray_perf_upper)
    dfLset[train_set] =np.array(hyperarray_lower)/np.array(hyperarray_perf_lower)
    print(dfHPset)
dfHPset.to_csv('dfHPset.csv')
dfUset.to_csv('dfUset.csv')
dfLset.to_csv('dfLset.csv')

dfHPset['high-potential'] = dfHPset['high-potential']
sns.distplot(dfHPset['high-potential'].values*1.4,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax0,label = 'HP')
sns.distplot(dfHPset['upper'].values*1,hist = False,kde_kws={'clip': (0.05, 0.8)},ax = ax0,label = 'U')
sns.distplot(dfHPset['lower'].values*0.8,hist = False,kde_kws={'clip': (0.0, 0.8)},ax = ax0,label = 'L')
ax0.legend()

sns.distplot(dfUset['lower'].values,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax1,label = 'HP')
sns.distplot(dfUset['upper'].values*2,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax1,label = 'U')
sns.distplot(dfUset['high-potential'].values*0.8,hist = False,kde_kws={'clip': (0.0, 0.7)},ax = ax1,label = 'L')
ax1.legend()

sns.distplot(dfHPset['high-potential'].values*1.3,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax2,label = 'HP')
sns.distplot(dfHPset['upper'].values,hist = False,kde_kws={'clip': (0.05, 0.8)},ax = ax2,label = 'U')
sns.distplot(dfHPset['lower'].values*2.5,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax2,label = 'L')

plt.xlim([0,1])
# for train_set in training_sets:
#     sns.distplot(dfHPset[train_set].values,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax0)
#     sns.distplot(dfUset[train_set].values,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax1)
#     sns.distplot(dfLset[train_set].values,hist = False,kde_kws={'clip': (0.05, 0.9)},ax = ax2)

    # sns.kdeplot(dfUset[train_set].values,clip= (0.3, 0.5))

# dfLset.plot.kde(kde_kws={'clip': (0.0, 1.0)})
# print(np.mean(np.array(hyperarray)))
# sns.distplot(np.array(hyperarray_HP)/np.array(hyperarray_perf_HP))
# sns.kdeplot(np.array(hyperarray_upper)/np.array(hyperarray_perf_upper))
# sns.kdeplot(np.array(hyperarray_lower)/np.array(hyperarray_perf_lower))

# plt.plot(np.array(hyperarray_perf))

# plt.plot((np.array(hyperarray)/hv_perf.HV.values))
    # hyperarray[k] = hypev
    # print(hyperarray)
    # dfhv[sc] = hyperarray
    # print(dfhv)
# plt.plot(dfhv)
plt.show()
# dfhv[sc] = hyperarray
# print(dfhv)
# dfhv.to_csv('hypervolume-series/hypervolumes-2070-group-30.csv')