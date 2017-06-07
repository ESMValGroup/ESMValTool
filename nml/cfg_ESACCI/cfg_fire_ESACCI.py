# This is a config file for CCI data and CMIP5 fire data diagnostics

# generell flags
regionalization = True
shape = "continents"
shapeNames = 2  # column of the name values

# flags for basic diagnostics
globmeants = True
mima_globmeants = [0, 1]
cmap_globmeants = 'RdYlGn_r'  # 'YlGn_r'
mima_ts = [0, 1]
mima_mts = mima_ts
portrait = True
globmeandiff = True
globmeandiff_p = True
mima_globmeandiff = [-1, 1]
mima_globmeandiff_r = [-1, 1]
trend = True
anomalytrend = True
trend_p = True
climatologies = True
hovmoeller = True
mima_hov = [0, 1]
mima_hovdiff = [-1, 1]

# flags for specific diagnostics

# plotting parameters
projection = {'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
