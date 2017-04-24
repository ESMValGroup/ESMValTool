# This is a config file for QA4ECV and CMIP5 data diagnostics

# generell flags
regionalization = True
shape = "continents"
shapeNames = 2  # column of the name values

# flags for basic diagnostics
globmeants = True
mima_globmeants = [0, 2.5]
cmap_globmeants = 'RdYlGn_r'  # 'YlGn_r'
mima_ts = [0, 2.5]
mima_mts = mima_ts
portrait = True
globmeandiff = True
globmeandiff_p = True
mima_globmeandiff = [-2.5, 2.5]
mima_globmeandiff_r = [-1, 1]
trend = True
anomalytrend = True
trend_p = True
climatologies = True

# flags for specific diagnostics

# plotting parameters
projection = {'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
