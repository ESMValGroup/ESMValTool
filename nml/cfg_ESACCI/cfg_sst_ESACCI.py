# This is a config file for CCI data and CMIP5 sea surface temperature
# diagnostics

# generell flags
regionalization = True
shape = "Seas_v"
shapeNames = 1  # column of the name values

# flags for basic diagnostics
globmeants = True
mima_globmeants = [255, 305]
cmap_globmeants = 'jet'
mima_ts = [288, 292]
mima_mts = [270, 310]
portrait = True
globmeandiff = True
mima_globmeandiff = [-10, 10]
mima_globmeandiff_r = [-0.03, 0.03]
trend = True
anomalytrend = True
trend_p = True
climatologies = True
hovmoeller = False
mima_hov = [0, 1]
mima_hovdiff = [-1, 1]

# flags for specific diagnostics
percentile = True
percentile_pars = [0.05, 0.95, 0.15]  # start,stop,increment
mima_percentile = [255, 305]

# plotting parameters
projection = {'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
