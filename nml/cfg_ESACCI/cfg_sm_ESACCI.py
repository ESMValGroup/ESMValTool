# This is a config file for CCI data and CMIP5 soil moisture diagnostics

# generell flags
regionalization = True
shape = "continents"
shapeNames = 2  # column of the name values
# start_year = 1988
# stop_year = 2000

# flags for basic diagnostics
globmeants = True
mima_globmeants = [0, 0.5]
cmap_globmeants = 'Blues'
mima_ts = [0, 0.5]
mima_mts = mima_ts
portrait = True
globmeandiff = True
mima_globmeandiff = [-0.25, 0.25]
mima_globmeandiff_r = [-1, 1]
trend = True
anomalytrend = True
trend_p = True
climatologies = True

# flags for specific diagnostics
percentile = True
percentile_pars = [0.05, 0.95, 0.15]  # start,stop,increment
mima_percentile = [0, 0.5]
anomaly = True

# plotting parameters
projection = {'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
