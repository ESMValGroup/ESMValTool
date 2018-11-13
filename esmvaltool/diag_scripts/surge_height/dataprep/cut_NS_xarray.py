import numpy as np
import xarray as xr

def cut_NS(psl, ua, va):
	ymin, ymax = [47.5, 65.6]
	xmin, xmax = [-13.5+180,13.5+180]

	## ------------------------------------
	## I. Interpolate to ERA-Interim grid 
	## ------------------------------------
        # remap manually
	#lons = np.arange(0,360,0.75) 
	#lats = np.arange(-89.75,90,0.75) 
	#psl = psl.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
	#ua = ua.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})
	#va = va.interp(lon=lons,lat=lats,method='linear', kwargs={'fill_value': None})

	# -------------------
	# III. Cut NS region
	# -------------------
	# >> ds.sel(lon=(ds.lon < xmax) | (ds.lon > xmin), lat=(ds.lat < ymax) | (ds.lat > ymin))
	pslNS = psl.sel(lon=slice(xmin,xmax),lat=slice(ymin,ymax))
		#psl.sel(lon=(psl.lon < xmax) | (psl.lon > xmin), 
		#		    lat=(psl.lat < ymax) | (psl.lat > ymin))
	uaNS = ua.sel(lon=slice(xmin,xmax),lat=slice(ymin,ymax))
		#ua.sel(lon=(ua.lon < xmax) | (ua.lon > xmin), 
		#		    lat=(ua.lat < ymax) | (ua.lat > ymin))
	vaNS = va.sel(lon=slice(xmin,xmax),lat=slice(ymin,ymax))
		#va.sel(lon=(va.lon < xmax) | (va.lon > xmin), 
		#		    lat=(va.lat < ymax) | (va.lat > ymin))

        return pslNS, uaNS, vaNS
