# This is a config file for CCI data and CMIP5 land cover diagnostics

# general flags for regionalization
regionalization = True
shape = "continents"
shapeNames = 2  # column of the name values

# flags for basic diagnostics
globmeants = True
mima_globmeants = [0, 100]
cmap_globmeants = 'YlGn'
mima_ts = [0, 100]
mima_mts = [0, 100]
portrait = True
globmeandiff = True
mima_globmeandiff = [-100, 100]
mima_globmeandiff_r = [-1, 1]
trend = True
anomalytrend = True
trend_p = True
climatologies = True
hovmoeller = True
mima_hov = [0, 100]
mima_hovdiff = [-100, 100]

# flags for specific diagnostics
single_years = True
mima_single_year = [-100, 100]
std_factor = 4

# plotting parameters
projection = {'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
