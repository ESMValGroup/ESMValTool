import xarray as xr
import os

def load_monmean_srgclim(stat,data_dir):
    filename = 'monanom_ERAintWAQUA_surge_1979-2016_speed.nc'
    
    ncpath = os.path.join(data_dir, filename)
    nc_srg_anom = xr.open_dataset(ncpath, 'r')
    srg_anom = nc_srg_anom.WAQUA_surge

    monmean_srgclim = {}
    for ifname in stat:
        monmean_srgclim[str(ifname)] = srg_anom.sel(stations=s)

        #ncpath = os.path.join(data_dir, filename + str(ifname) + '.nc')
        #nc = Dataset(ncpath, 'r')
	#nc.variables['WAQUA_surge'][:]
        #nc.close()
    #
    return monmean_srgclim
