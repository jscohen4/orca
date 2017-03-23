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

plt.subplot(2,1,1)
plt.plot(shasta.S)
plt.plot(df.SHA_storage.values)

plt.subplot(2,1,2)
plt.plot(shasta.R)
plt.plot(df.SHA_out.values * cfs_tafd)
plt.show()