# This is a config file for QA4ECV and CMIP5 data diagnostics

# generell flags
regionalization = False
shape = "continents"
shapeNames = 2  # column of the name values

# flags for basic diagnostics
globmeants = True
mima_globmeants = [0, 2.5]
cmap_globmeants = 'gist_stern'  # 'YlGn_r'
mima_ts = [0, 2.5]
mima_mts = mima_ts
portrait = False
globmeandiff = False
globmeandiff_p = False
mima_globmeandiff = [-2.5, 2.5]
mima_globmeandiff_r = [-1, 1]
trend = False
anomalytrend = False
trend_p = False
climatologies = True
hovmoeller = True
mima_hov = [0, 1]
mima_hovdiff = [-1, 1]

# flags for specific diagnostics

# plotting parameters
projection = {'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
