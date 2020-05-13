import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.stats import norm,skew
import seaborn as sns
import inspyred.ec.analysis as pyred
import pickle
import math

def hypervolume(pareto_set, reference_point=None):
    #Adapted from inspyred: Bio-inspired Algorithms in Python
    #https://pythonhosted.org/inspyred/
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

spill = pd.read_csv('data/baseline-results/baseline_spill.csv', parse_dates = True,index_col = 0)
carryover = pd.read_csv('data/baseline-results/baseline_carryover.csv', parse_dates = True,index_col = 0)
delta = pd.read_csv('data/baseline-results/baseline_delta.csv', parse_dates = True,index_col = 0)
pumping = pd.read_csv('data/baseline-results/baseline_pumping.csv', parse_dates = True,index_col = 0)
hydro = pd.read_csv('data/baseline-results/baseline_energy.csv', parse_dates = True,index_col = 0)

for j,sc in enumerate(scenarios):
    print(j)
    hyperarray = np.zeros(len(yearstartlist))
    print(hyperarray)
    for k, (yearstart, yearend) in enumerate(zip(yearstartlist, yearendlist)):
        print(yearstart)
        spillbase = spill[spill.index.year >= yearstart].cumsum()
        carryoverbase = carryover[carryover.index.year >= yearstart].cumsum()
        deltabase = delta[delta.index.year >= yearstart].cumsum()
        pumpingbase = pumping[pumping.index.year >= yearstart].cumsum()
        hydrobase = hydro[hydro.index.year >= yearstart].cumsum()

        with open('pickle_files_adapt_pot/10-year-chunks/objectives-%s/objectives-sc%s.pkl'%(yearstart,j), 'rb') as f:
            obj_raw = pickle.load(f)
        obj_raw = np.array(obj_raw)
        obj0 = obj_raw[:,0]
        obj1 = obj_raw[:,1]*-1
        obj2 = obj_raw[:,2]*-1
        obj3 = obj_raw[:,3]*-1
        obj4 = obj_raw[:,4]*-1
        objnorm = np.zeros([len(obj0),4])

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

        for i,obj in enumerate(obj4):
          objs = (obj - obj4.min(axis=0)) / (obj4.max(axis=0) - obj4.min(axis=0))
          objnorm[i][4] = objs

    for i,obj in enumerate(obj4):
      objs = (obj - obj4.min(axis=0)) / (obj4.max(axis=0) - obj4.min(axis=0))
      objnorm[i][4] = objs

        baseline0 = spillbase[spillbase.index == '%s-10-01'%yearend][sc].values[0]
        baseline0 = (baseline0 - obj0.max(axis=0)) / (obj0.min(axis=0) - obj0.max(axis=0))

        baseline1 = carryoverbase[carryoverbase.index == '%s-10-01'%yearend][sc].values[0]
        baseline1 = (baseline1 - obj1.min(axis=0)) / (obj1.max(axis=0) - obj1.min(axis=0))

        baseline2 = deltabase[deltabase.index == '%s-10-01'%yearend][sc].values[0]
        baseline2 = (baseline2 - obj2.min(axis=0)) / (obj2.max(axis=0) - obj2.min(axis=0))

        baseline3 = pumpingbase[pumpingbase.index == '%s-10-01'%yearend][sc].values[0]
        baseline3 = (baseline3 - obj3.min(axis=0)) / (obj3.max(axis=0) - obj3.min(axis=0))

        baseline3 = pumpingbase[pumpingbase.index == '%s-10-01'%yearend][sc].values[0]
        baseline3 = (baseline3 - obj3.min(axis=0)) / (obj3.max(axis=0) - obj3.min(axis=0))

        baseline4 = hydrobase[hydrobase.index == '%s-10-01'%yearend][sc].values[0]
        baseline4 = (baseline4 - obj3.min(axis=0)) / (obj4.max(axis=0) - obj4.min(axis=0))

        objnorm = -objnorm
        ref_pt = [-abs(baseline0),-abs(baseline1),-abs(baseline2),-abs(baseline3),-abs(baseline4)]

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
        hyperarray[k] = hypev
    print(hyperarray)
    dfhv[sc] = hyperarray
    print(dfhv)
dfhv.to_csv('adaptation-potential.csv')