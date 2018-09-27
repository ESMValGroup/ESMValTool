def Xtrms(psl,ua,va):
        import numpy as np
	import xarray as xr
	
	dmin_psl = psl.resample(time='D').min('time')
	dmax_ua  = ua.resample(time='D').max('time')
	dmax_va  = va.resample(time='D').max('time')
	
        return dmin_psl, dmax_ua, dmax_va
