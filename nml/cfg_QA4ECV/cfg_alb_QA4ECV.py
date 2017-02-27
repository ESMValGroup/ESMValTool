# This is a config file for QA4ECV and CMIP5 data diagnostics

# generell flags
regionalization = True
shape = "climes"
shapeNames = 2 #column of the name values 
#start_year = 1988
#stop_year = 2005

#define albedo cmap (http://jdherman.github.io/colormap/)
import matplotlib as mpl
import numpy as np

C=np.matrix([[15,15,20],
[60,75,50],
[75,100,60],#new
[90,125,70],
[110,145,85],
[125,130,100],
[145,150,105],
[170,165,115],
[190,180,135],
[205,200,170],#new
[225,225,210],
[235,235,250]])

albedocmap = mpl.colors.ListedColormap(C/255.0)

# flags for basic diagnostics
globmeants = True
mima_globmeants=[0,0.6]
cmap_globmeants=albedocmap #'YlGn_r'
mima_ts=[0,1]
mima_mts=mima_ts
portrait = True
globmeandiff = True
globmeandiff_p = True
mima_globmeandiff=[-0.25,0.25]
mima_globmeandiff_r=[-1,1]
trend = True
anomalytrend = True
trend_p = True
climatologies = True

# flags for specific diagnostics

# percentile parameters

#plotting parameters
projection={'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}