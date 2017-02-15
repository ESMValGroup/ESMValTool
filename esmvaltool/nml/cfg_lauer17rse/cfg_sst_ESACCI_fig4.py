# This is a config file for CCI data and CMIP5 sea surface temperature diagnostics

# generell flags
regionalization = True
shape = "Seas_v"
shapeNames = 1 #column of the name values 

# flags for basic diagnostics
globmeants = False
mima_globmeants=[255,305]
cmap_globmeants='jet'
mima_ts=[288,292]
mima_mts=[270,310]
portrait = True
globmeandiff =True
mima_globmeandiff=[-10,10]
mima_globmeandiff_r=[-0.03,0.03]
trend = False# True
anomalytrend =False#True
trend_p=False

# flags for specific diagnostics
percentile = False#True
mima_percentile=[255,305]

# percentile parameters
percentile_pars = [0,1,0.25] #start,stop,increment

#plotting parameters
projection={'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
