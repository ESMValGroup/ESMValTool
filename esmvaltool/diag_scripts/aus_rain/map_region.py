"""analysis plots.

###############################################################
based on...
landcover/landcover.py
Authors ESMValToolV1 Version
    lejeune_quentin
Port to ESMValTool Version 2
    crezee_bas
###############################################################

"""


import logging
import os
from pathlib import Path
from cartopy import crs  
import xarray as xr
from shapely.geometry import Polygon
from cartopy.io import shapereader
import fiona

import matplotlib.pyplot as plt
import numpy as np


from esmvaltool.diag_scripts.shared import (group_metadata,
                                            run_diagnostic,
                                            ProvenanceLogger,
                                            get_plot_filename)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))



def get_country_polygon(country):
    resolution = '10m'
    category = 'cultural'
    name = 'admin_0_countries'

    shpfilename = shapereader.natural_earth(resolution, category, name)
    f = fiona.open(shpfilename)
    reader = shapereader.Reader(shpfilename)
    records = list(reader.records())

    for record in records:
        if record.attributes["ADMIN"] == country:
            return [record.geometry]

## exts read from bounds of xr ds
def get_ds_extents(ds):
    ds1 = ds.drop('time').reset_coords().squeeze()

    lon = ds1.coords['lon'].values
    lat = ds1.coords['lat'].values
    exts = [lon.min(),lon.max(),lat.min(),lat.max()]
    return exts


def rect_from_bound(xmin, xmax, ymin, ymax):
    """Returns list of (x,y)'s for a rectangle"""
    xs = [xmax, xmin, xmin, xmax, xmax]
    ys = [ymax, ymax, ymin, ymin, ymax]
    return [(x, y) for x, y in zip(xs, ys)]

def plot_ocean_map(nf, name, mask,preproc):
    """
    use xarray.. 
    """
    ## open dataset, if mask is land, mask is the boundary ..buffer?
    ds = xr.open_dataset(nf)

    poly = get_country_polygon("Australia")
    if mask == 'land':
        msk = poly[0].simplify(0.1)
    elif mask == 'sea':
        exts = get_ds_extents(ds)
        msk = Polygon(rect_from_bound(*exts)).difference(poly[0].simplify(0.1))

    fig = plt.figure()
    cx = plt.axes(projection=crs.PlateCarree()) ##
    cx.axes.coastlines()

    if mask in ['land','sea']:
        msk_stm  = crs.PlateCarree().project_geometry (msk, crs.PlateCarree())  # project geometry to the projection used by stamen
        cx.add_geometries(msk_stm, crs.PlateCarree(), zorder=12, facecolor='white', edgecolor='none', alpha=0.9)

    ticks ={'ticks':None}
    var = list(ds.keys())[0]
    cmap='RdBu'
    vminmax={}
    if preproc=='ocean_decadal': ## preprocessor_decadal
        ds = trend_val(ds,var)
        var = 'val_trend'
    
    if preproc=='decilemap':
        ds= decile_mapdata(ds,var)
        var = 'val'
        vminmax = {'vmin':0,'vmax':10}
        ticks ={'ticks':list(range(1,10,1))} # if decile category map
        cmap = plt.get_cmap('RdBu', 10)
        
    if 'time' in list(ds.dims):
        logger.info(" time dimension exists, sum along time")
        ds = ds.sum(dim='time') ##sum time
    
    # cx = cx.axes.contourf(ds['lon'],ds['lat'],ds[var],transform=crs.PlateCarree(),cmap='coolwarm') #contours?
    cx = cx.axes.pcolormesh(ds['lon'],ds['lat'],ds[var],transform=crs.PlateCarree(),**vminmax,cmap=cmap)
    plt.colorbar(cx, shrink=0.7,**ticks)
    
    # naem from variable group? get time? model?
    name = Path(nf).stem
    cx.axes.set_title(f'{name}')
    
    return fig

def trend_val(ds,var): ##decadal

    ds = ds.swap_dims({'time':'decade'})
    varls = list(ds.data_vars)
    varls.remove(var)
    ds = ds.drop_vars(varls)

    ds_temp = ds.polyfit(dim='decade',deg=1)
    
    ds_temp = ds_temp.sel(degree=1) ##value degree 0? ..coeff
    ds_temp = ds_temp.drop('degree').reset_coords().squeeze()
    ds_out = ds.assign(val_trend=ds_temp[f'{var}_polyfit_coefficients'])

    return ds_out

## full compare to recent ...difference ranges? (means) ...count above/below?
def decile_mapdata(ds,var):

    ## preprocessed -yearly sum, 
    ds = ds.swap_dims({'time':'year'})
    ## range -min/max and percentile calc ...sorted array for each grid pt? ##dataset for each quantile
    # ref.. https://stackoverflow.com/questions/62698837/calculating-percentile-for-each-gridpoint-in-xarray
    qt_dims = 'year'
    qt_values = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) #test, add each
    ds_qt = ds.quantile(qt_values,qt_dims) 
    
    ## decile map calc recent vs quantile of full ..mean of recent ~20years##
    ds = ds.sel(year=slice(1990, 2010)) #mean
    ds = ds.mean(dim='year')

    ## apply value logic over indexes? 
    for i in range(len(qt_values)): 
        da_temp = np.greater(ds[var],ds_qt[var].sel(quantile=qt_values[i]))
        #  #assign to temp, add arrays
        if i == 0:
            ds_count = da_temp.astype(int) ##set 0 to nan..
        else:
            ds_count = ds_count + da_temp.astype(int)
    ## increment, if greater than than +=1 ... 
    ## assign data array to dataset to map in same method
    ds_count = ds.assign(val=ds_count.where(~ds[var].isnull()))

    return ds_count #map

def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    msk = cfg['mask']
    logger.info(f'test... params for script.. mask {msk}')
    
    # Example of how to loop over variables/datasets in alphabetical order
    for dataset in input_data:
        # Load the data
        logging.info(f"{type(dataset),dataset['filename'] }")
        input_file = dataset['filename']
        name = dataset['variable_group']


        preproc = dataset['preprocessor']
        fig = plot_ocean_map(input_file,name,msk,preproc)

        # Save output
        output_file = Path(input_file).stem.replace(dataset['short_name'], name)
        output_path = get_plot_filename(output_file, cfg)
        fig.savefig(output_path)
    

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
