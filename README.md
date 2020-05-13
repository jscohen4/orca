# Operation of Reservoirs in California (ORCA)

A simulation mode incorporating CVP and SWP Delta exports, reservoir releases, and snow-pack to streamflow forecasts. This branch simulates historical scenarios with data from cdec cdec or data from 2006-2017 generated with WRF and and routed with WRF-WRF Hydro .

## Requirements:
[NumPy](http://www.numpy.org/), [pandas](http://pandas.pydata.org/), [Matplotlib](http://matplotlib.org/), [Scipy](http://www.scipy.org/), and [scikit-learn](http://scikit-learn.org/).

## Running ORCA

### Running ORCA in cdec mode:
1. In main.py, set maurer to False.
2. Choose whether to calculate R2s and plot resuts (set variables to True or False. If running any data processing scripts, set process_hist_data to True. If downloading up-to-date historical data, set cdec to True. Set hist_indices to True if re-processing new data. Set hist_forcast to True if re-running forecasts.
3. Run main.py and results will be in historical_runs_data folder.

### Running ORCA in with WRF data:
1. In main.py, set WRF to True.
2.  Set WRF_indices to True if re-processing new data. Set WRF_forecasts to True if re-running forecasts. Choose whether or not to plot.
3. Run main.py and results will be in WRF_runs_data folder.

## License
Copyright (C) 2020 Cohen, J.S., Zeff, H.B., & Herman, J.D. Released under the [MIT license](LICENSE.md).
