import xarray as xr
import logging
import os
import numpy as np
from xclim.indices.fire import fire_weather_ufunc, fire_season
from xclim.core.units import convert_units_to
from xclim.indices import saturation_vapor_pressure

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

def load_xr_da(metadata, var):

    var_fil = select_metadata(metadata, short_name=var)[0]['filename']
    var_xr = xr.open_dataset(var_fil)[var]

    # to avoid issues with merging later on, drop  height which is scalar
    var_coords = var_xr.coords
    if 'height' in var_coords:
        var_xr = var_xr.drop('height')

    return var_xr


# def hursmin_fun(hussmax_sat_xr, hussmin_sat_xr, hurs_xr):
    
#     # calculating the mean saturation wv pressure 
#     hussday_sat_xr = (hussmin_sat_xr + hussmax_sat_xr)/2
#     # calculate mean wv pressure
#     huss_xr = hurs_xr/100 * hussday_sat_xr 
#     # assuming hursmin occurs during tasmax 
#     hursmin_xr = huss_xr/hussmax_sat_xr * 100

#     return hursmin_xr


def huss_fun(hussmax_sat_xr, hussmin_sat_xr, hurs_xr):
    
    # calculating the mean saturation wv pressure 
    hussday_sat_xr = (hussmin_sat_xr + hussmax_sat_xr)/2
    # calculate mean wv pressure
    huss_xr = hurs_xr/100 * hussday_sat_xr 

    return huss_xr


def hursmin_fun(hussmax_sat_xr, huss_xr):
    
    # assuming hursmin occurs during tasmax 
    hursmin_xr = huss_xr/hussmax_sat_xr * 100

    return hursmin_xr


def concatenate_xr_dt(xr_list):

    # needed to merge xr_da 
    for n_da in range(len(xr_list)-1):
        if np.all(xr_list[n_da]['time'] != xr_list[-1]['time']):
            xr_list[n_da]['time'] = xr_list[-1]['time']

    data_xr_dt = xr.merge(xr_list)

    return data_xr_dt 


def derive_hursmin(hurs_xr, tasmax_xr, tasmin_xr, dataset, cfg):

    hussmax_sat_xr = xr.map_blocks(saturation_vapor_pressure, 
                     tasmax_xr, kwargs={"method": "wmo08"}).rename("hussmax") 
    hussmin_sat_xr = xr.map_blocks(saturation_vapor_pressure, 
                     tasmin_xr, kwargs={"method": "wmo08"}).rename("hussmin")

    # calculate water vapor saturation pressure
    huss_xr = xr.apply_ufunc(huss_fun, hussmax_sat_xr, hussmin_sat_xr, 
                   hurs_xr, dask="parallelized", vectorize=True).rename("huss")

    # deriving hursmin
    hursmin_xr = xr.apply_ufunc(hursmin_fun, hussmax_sat_xr, huss_xr,
                              dask="parallelized", vectorize=True).rename("hursmin")
    hursmin_xr = xr.where(hursmin_xr > 100, 100, hursmin_xr)

    # saving vapor pressure deficit have as an area burned predictor
    vpd_filename = determine_filename(dataset, 'vpd', cfg)
    vpd_xr = hussmax_sat_xr - huss_xr
    vpd_xr.rename("vpd").to_netcdf(os.path.join(cfg['work_dir'], vpd_filename))

    # saving the 'mean' temperature as an area burned predictor
    tmean_filename = determine_filename(dataset, 'tmean', cfg)
    tmean_xr = (tasmax_xr + tasmin_xr)/2
    tmean_xr.rename("tmean").to_netcdf(os.path.join(cfg['work_dir'], tmean_filename))

    return hursmin_xr 


def determine_filename(metadata, fwi_idx, cfg):

    region = cfg['region']

    mtd = metadata[0]

    flnm = mtd['filename'].split('/')[-1]
    varn = mtd['short_name']

    filename = flnm.replace(varn, fwi_idx+'_'+region)

    return filename

def create_forest_mask(cfg, data_xr):

    # loading forested area mask, filename comes from the recipe
    mask_f = os.path.join(cfg['auxiliary_data_dir'], cfg['forest_mask'])
    mask_xr = xr.open_dataarray(mask_f).sel(band=1)

    # mask comes in m, need to put on lat/lon grid
    lat_lon_mask = mask_xr.rio.reproject("EPSG:4326")
    # adding projection to the data data
    data_xr.rio.write_crs(4326, inplace=True)

    # matching the projection of mask and the data
    mask_matched = lat_lon_mask.rio.reproject_match(data_xr.rio.set_spatial_dims(x_dim='lon', y_dim='lat'))
    # renaming the dimensions to match data
    mask_matched = mask_matched.rename({'x': 'lon', 'y': 'lat'})

    return mask_matched

def main(cfg):

    input_data = cfg['input_data']
    work_dir = cfg['work_dir']

    datasets = group_metadata(input_data.values(), 'alias', sort=True)

    for dataset in datasets.keys(): 
        # loading data and converting the data to the units 
        # needed by xclim to calculate FWIs
        hurs_xr = load_xr_da(datasets[dataset], var='hurs')
        sfcwind_xr = load_xr_da(datasets[dataset], var='sfcWind')
        sfcwind_xr = convert_units_to(sfcwind_xr, 'km/h')
        pr_xr = load_xr_da(datasets[dataset], var='pr')
        pr_xr = convert_units_to(pr_xr, 'mm/d')
        tasmax_xr = load_xr_da(datasets[dataset], var='tasmax')
        tasmax_xr = convert_units_to(tasmax_xr, 'degC')
        tasmin_xr = load_xr_da(datasets[dataset], var='tasmin')
        tasmin_xr = convert_units_to(tasmin_xr, 'degC')

        # creating forested area mask from the GFS geotiff 
        forest_mask = create_forest_mask(cfg, tasmax_xr)

        # creating one dataset from all the loaded data
        xr_list = [hurs_xr, sfcwind_xr, pr_xr, tasmax_xr, tasmin_xr]
        raw_dt = concatenate_xr_dt(xr_list)

        hursmin_xr = derive_hursmin(raw_dt.hurs, raw_dt.tasmax, raw_dt.tasmin, 
                                    datasets[dataset], cfg)

        raw_dt = xr.merge([raw_dt, hursmin_xr])
    
        # deriving fire season mask
        fire_season_mask = fire_season(raw_dt.tasmax, method='WF93', freq=None,  
                                    temp_start_thresh='12 degC', temp_end_thresh='5 degC', 
                                    temp_condition_days=3).rename('fire_season_mask') 

        # saving fire season mask, first determing the filename
        fwi_filename = determine_filename(datasets[dataset], 'FSM', cfg)
        fire_season_mask.to_netcdf(os.path.join(work_dir, fwi_filename))

        fwi_dict = fire_weather_ufunc(tas = raw_dt.tasmax,
                                    pr = raw_dt.pr,
                                    sfcWind = raw_dt.sfcWind,
                                    hurs = raw_dt.hursmin,
                                    lat = raw_dt.lat,
                                    season_mask = fire_season_mask, 
                                    overwintering = True, 
                                    carry_over_fraction=0.75, 
                                    wetting_efficiency_fraction=0.75, 
                                    dry_start = "CFS", 
                                    prec_thresh=1.5, 
                                    dmc_dry_factor=1.2)

        # no need to save winter precip, so removing from the dic 
        del(fwi_dict['winter_pr'])

        for fwi_idx in fwi_dict.keys():
            fwi_filename = determine_filename(datasets[dataset], fwi_idx, cfg)
            fwi_dict[fwi_idx].to_netcdf(os.path.join(work_dir, fwi_filename))
            assert np.all(fwi_dict[fwi_idx].lat == forest_mask.lat), 'Latitudes of mask and FWI files are different'
            assert np.all(fwi_dict[fwi_idx].lon == forest_mask.lon), 'Longitudes of mask and FWI files are different'
            # masking fwi file to the forest mask
            fwi_dict[fwi_idx].where(forest_mask>=0.1).to_netcdf(os.path.join(work_dir, 'masked_'+fwi_filename))

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)