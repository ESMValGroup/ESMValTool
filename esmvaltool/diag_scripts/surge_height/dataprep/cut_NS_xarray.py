def cut_NS(psl, ua, va):
	import numpy as np
	import xarray as xr

	ymin, ymax = [47.5, 65.6]
	xmin, xmax = [-13.5,13.5]

	## ------------------------------------
	## I. Interpolate to ERA-Interim grid 
	## ------------------------------------
	## (probably unnecessary as grid indicated in request script, however will that be native grid,
	##  i.e. Gaussian, or is it the 0.75x0.75 standard grid?) -> else load target field using DataArray
	## "interp_like()" method is a useful shortcut. This method interpolates an xarray object onto the
	## coordinates of another xarray object.
	## >> interpolated = DATA_ARRAY1.interp_like(DATA_ARRAY2)
	# -------------------
	# III. Cut NS region
	# -------------------
	# >> ds.sel(lon=(ds.lon < xmax) | (ds.lon > xmin), lat=(ds.lat < ymax) | (ds.lat > ymin))
	pslNS = psl.sel(lon=(psl.lon < xmax) | (psl.lon > xmin), 
				    lat=(psl.lat < ymax) | (psl.lat > ymin))
	uaNS = ua.sel(lon=(ua.lon < xmax) | (ua.lon > xmin), 
				    lat=(ua.lat < ymax) | (ua.lat > ymin))
	vaNS = va.sel(lon=(va.lon < xmax) | (va.lon > xmin), 
				    lat=(va.lat < ymax) | (va.lat > ymin))

        return pslNS, uaNS, vaNS
