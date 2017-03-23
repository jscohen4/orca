import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cord import *

datafile = 'cord/data/cord-data.csv'

df = pd.read_csv(datafile, index_col=0, parse_dates=True)
T = len(df)

shasta = Reservoir(df, 'SHA')

for t in range(1,T):
  shasta.step(t)

results = shasta.results_as_df(df.index)

plt.subplot(2,1,1)
h = results.SHA_storage.plot()
df.SHA_storage.plot(ax=h)

plt.subplot(2,1,2)
h = results.SHA_out.plot()
(df.SHA_out * cfs_tafd).plot(ax=h)
plt.show()
