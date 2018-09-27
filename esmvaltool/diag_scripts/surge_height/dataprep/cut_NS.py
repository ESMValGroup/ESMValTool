def cut_NS(psl, ua, va, dates):
	import numpy as np
	import xarray as ax
	from netCDF4 import num2date
	import xarray as xr

	ymin, ymax = [47.5, 65.6]
	xmin, xmax = [-13.5,13.5]
	lons = np.arange(xmin,xmax+0.25,0.75)	
	lats = np.arange(ymin,ymax+0.25,0.75)

	# -----------------------------
	# I. Convert list to DataArray
	# -----------------------------
	# >> xarray.DataArray(data, coords=None, dims=None, name=None, attrs=None, encoding=None, fastpath=False)
	# NR: Don't know if what kind of input we get, if we get xarray.DataArrays we can use this but else how
	#     would we get the units and calendar attributes? Also if we get xarray.DataArrays we don't need this
	#     conversion at all
	#tencoding  =  {'units': u'hours since 2035-01-01 06:00:00', 'calendar': u'standard'}
	tcoor = xr.Coordinate('time',dates_SLP)#,encoding=tencoding)
	coords = {'time':tcoor,'lon':lons,'lat':lats}
	psl_array = xr.DataArray(psl, coords=coords,dims=('time','lat','lon'))
	ua_array  = xr.DataArray(ua, coords=coords,dims=('time','lat','lon'))
	va_array  = xr.DataArray(va, coords=coords,dims=('time','lat','lon'))
	## ------------------------------------
	## II. Interpolate to ERA-Interim grid 
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
	pslNS_array = psl_array.sel(lon=(psl_array.lon < xmax) | (psl_array.lon > xmin), 
				    lat=(psl_array.lat < ymax) | (psl_array.lat > ymin))
	uaNS_array = ua_array.sel(lon=(ua_array.lon < xmax) | (ua_array.lon > xmin), 
				    lat=(ua_array.lat < ymax) | (ua_array.lat > ymin))
	vaNS_array = va_array.sel(lon=(va_array.lon < xmax) | (va_array.lon > xmin), 
				    lat=(va_array.lat < ymax) | (va_array.lat > ymin))
	# -------------------------
	# IV. Convert back to list 
	# -------------------------
	# >> lst = DataArray.values
	pslNS = pslNS_array.values
	vaNS = vaNS_array.values
	uaNS = uaNS_array.values

        return pslNS, uaNS, vaNS
