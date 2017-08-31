## Operation of Reservoirs in California (ORCA)

A simulation model. 

**Requirements:** [NumPy](http://www.numpy.org/), [pandas](http://pandas.pydata.org/), and [Matplotlib](http://matplotlib.org/).

Still in active development, not stable.

### Current calibration r^2 values (8/11/17):
These values will change as new input data is downloaded from CDEC, even if the model itself does not change.

| Timestep:  | Daily | Weekly | Monthly | Water Year |
| ---------- | ----- | ------ | ------- | ---------- |
| HRO_pump | 0.23 | 0.27 | 0.34 | 0.84 |
| TRP_pump | 0.34 | 0.40 | 0.49 | 0.80 |
| Combined | 0.43 | 0.48 | 0.60 | 0.89 |
| SHA_stor | 0.82 | 0.81 | 0.81 | n/a |
| SHA_out  | 0.62 | 0.70 | 0.79 | n/a | 
| FOL_stor | 0.74 | 0.74 | 0.78 | n/a |
| FOL_out  | 0.82 | 0.87 | 0.91 | n/a |
| ORO_stor | 0.89 | 0.88 | 0.88 | n/a |
| ORO_out  | 0.63 | 0.71 | 0.84 | n/a |
| DEL_in   | 0.95 | 0.96 | 0.98 |  n/a |
| DEL_out  | 0.95 | 0.96 | 0.98 | n/a |

### License
Copyright (C) 2017 CALFEWS Team. Released under the [MIT license](LICENSE.md).
