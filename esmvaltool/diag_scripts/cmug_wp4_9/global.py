"""
CMUG WP4.9
Use this to make plots from GLOBAL data
"""

import logging

import iris
import iris.coord_categorisation as icc
import iris.plot as iplt
import matplotlib as mpl

import matplotlib.pyplot as plt
import numpy as np
import cf_units as Unit


from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)


scale_factor = 0.001 # cant find this in the files, * this by all CCI values to get actual value

def _get_input_cubes(metadata):
    """Load the data files into cubes.
    Based on the hydrology diagnostic.
    Inputs:
    metadata = List of dictionaries made from the preprocessor config
    Outputs:
    inputs = Dictionary of cubes
    ancestors = Dictionary of filename information
    """
    # inputs = {}
    # ancestors = {}
    # for attributes in metadata:
    #     short_name = attributes['short_name']
    #     filename = attributes['filename']
    #     logger.info("Loading variable %s", short_name)
    #     cube = iris.load_cube(filename)
    #     cube.attributes.clear()
    #     inputs[short_name] = cube
    #     ancestors[short_name] = [filename]

    inputs = {}
    ancestors = {}
    print(metadata)
    for attributes in metadata:
        print(attributes)
        short_name = attributes['short_name']
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        
        try:
            key_name = f"{short_name}_{attributes['ensemble']}"
        except:
            key_name = short_name

        inputs[key_name] = cube
        ancestors[key_name] = [filename]
        
        print(inputs)
        print(ancestors)

    return inputs, ancestors

def _get_provenance_record(attributes, ancestor_files):
    """Create the provenance record dictionary.
    Inputs:
    attributes = dictionary of ensembles/models used, the region bounds
                 and years of data used.
    ancestor_files = list of data files used by the diagnostic.
    Outputs:
    record = dictionary of provenance records.
    """
    caption = "Timeseries of ESA CCI LST difference to mean of "\
        + "model ensembles calculated over region bounded by latitude "\
        + "{lat_south} to {lat_north}, longitude {lon_west} to {lon_east} "\
        + "and for model/ensembles {ensembles}. "\
        + "Shown for years {start_year} to {end_year}.".format(**attributes)

    record = {
        'caption': caption,
        'statistics': ['mean', 'stddev'],
        'domains': ['reg'],
        'plot_types': ['times'],
        'authors': ['king_robert'],
        # 'references': [],
        'ancestors': ancestor_files
    }

    return record


def _diagnostic(config):
    """Perform the control for the ESA CCI LST diagnostic
       with uncertainities included
    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.
    Returns
    -------
    figures made by make_plots.
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()

    loaded_data = {}
    ancestor_list = []
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes

    # loaded data is a nested dictionary
    # KEY1 model ESACCI-LST or something else
    # KEY2 is ts, the surface temperature
    
    # need to work out how multiple ensembles are passed in

    # ie loaded_data['ESACCI-LST']['ts'] is the CCI cube
    #    loaded_data['MultiModelMean']['ts'] is CMIP6 data, emsemble means
    #    similarly dor Std, see preprocessor

    # The Diagnostic uses CCI - MODEL

    # CMIP data had 360 day calendar, CCI data has 365 day calendar
    # Assume the loaded data is all the same shape
    print('LOADED DATA:')
    print(f'{loaded_data=}')

    
    ### There will be some cube manipulation todo
    # add year and month coords
    for KEY in loaded_data.keys():
        for ITEM in loaded_data[KEY].keys():
            icc.add_year(loaded_data[KEY][ITEM], 'time')
            icc.add_month_number(loaded_data[KEY][ITEM], 'time')

    # make model ensemble means
    models = []
    for KEY in loaded_data.keys():
        if KEY == 'CMUG_WP4_9' or KEY == 'ESACCI_LST_UNCERTS':
            continue
        models.append(KEY)

    print(f'{models=}')

    model_mean_tair = {}
    model_mean_lst  = {}
    for MODEL in models:
        tair_cubelist = iris.cube.CubeList()
        lst_cubelist  = iris.cube.CubeList()
        for KEY in loaded_data[MODEL].keys():
            if 'ts' in KEY:
                enumber = len(lst_cubelist)
                e_coord = iris.coords.AuxCoord(enumber, var_name='ensemble_number')
                this_cube = loaded_data[MODEL][KEY]
                this_cube.add_aux_coord(e_coord)
                lst_cubelist.append(this_cube)

            if 'tas' in KEY:
                enumber = len(tair_cubelist)
                e_coord = iris.coords.AuxCoord(enumber, var_name='ensemble_number')
                this_cube = loaded_data[MODEL][KEY]
                this_cube.add_aux_coord(e_coord)
                tair_cubelist.append(this_cube)

        model_mean_tair[MODEL] = tair_cubelist.merge_cube()
        model_mean_tair[MODEL] = model_mean_tair[MODEL].collapsed('ensemble_number', iris.analysis.MEAN) 
 
        model_mean_lst[MODEL]  = lst_cubelist.merge_cube()
        model_mean_lst[MODEL]  = model_mean_lst[MODEL].collapsed('ensemble_number', iris.analysis.MEAN) 


    # make mean lst from day and night
    cci_lst = loaded_data['ESACCI_LST_UNCERTS']['tsDay'] + loaded_data['ESACCI_LST_UNCERTS']['tsNight']
    cci_lst = cci_lst/2

    ### regrid obs to each model

    
    era5l_tair = loaded_data['CMUG_WP4_9']['tas']

    #### regrid tair to model grid

    # make climatologies
    #model_lst_clim  = model_lst.aggregated_by('month_number', iris.analysis.MEAN)
    #model_tair_clim = model_lst.aggregated_by('month_number', iris.analysis.MEAN)
    cci_lst_clim    = cci_lst.aggregated_by('month_number', iris.analysis.MEAN)
    era5l_tair_clim = era5l_tair.aggregated_by('month_number', iris.analysis.MEAN)

    # plots
    #make_plot_global_clim_maps(model_lst, model_ta, cci_lst, era5l_tair, config)
    # make this so either lst or tair passed in, will need vmin vmax as parameters #####
    make_plot_global_clim_maps(None, None, cci_lst, era5l_tair, config)

    # NEXT doFrance LC plots

    # more plotting functions go here

#     # Provenance
#     # Get this information form the data cubes
#     # data_attributes = {}
#     # data_attributes['start_year'] = lst_diff_cube.coord('time').units.num2date(
#     #     lst_diff_cube.coord('time').points)[0].year
#     # data_attributes['end_year'] = lst_diff_cube.coord('time').units.num2date(
#     #     lst_diff_cube.coord('time').points)[-1].year
#     # data_attributes['lat_south'] = lst_diff_cube.coord('latitude').bounds[0][0]
#     # data_attributes['lat_north'] = lst_diff_cube.coord('latitude').bounds[0][1]
#     # data_attributes['lon_west'] = lst_diff_cube.coord('longitude').bounds[0][0]
#     # data_attributes['lon_east'] = lst_diff_cube.coord('longitude').bounds[0][1]
#     # data_attributes['ensembles'] = ''

#     # for item in input_metadata:
#     #     if 'ESACCI' in item['alias'] or 'MultiModel' in item[
#     #             'alias'] or 'OBS' in item['alias']:
#     #         continue
#     #     data_attributes['ensembles'] += "%s " % item['alias']

#     # record = _get_provenance_record(data_attributes, ancestor_list)
#     # for file in ['%s/timeseries.png' % config['plot_dir']]:
#     #     with ProvenanceLogger(config) as provenance_logger:
#     #         provenance_logger.log(file, record)



def make_plot_global_clim_maps(model_lst_clim, model_tair_clim, 
                               cci_lst_clim, era_tair_clim,
                               config):
    """Create and save the output figure.
    Figure of airT from ERA5 and differences from airT in CMIP6 models (historic climatology – 2003-13 obs and model years)
    Figure of LST from CCI LST and differences from LST in CMIP6 models (historic climatology – 2003-2013 obs and model years)
    These could be for individual months or whole year TBC

    Inputs:
   
    config = The config dictionary from the preprocessor
    Outputs:
    Saved figure
    """
    outpath = config['plot_dir']


    # First figure
    # Tair from ERA5-land and ERA(OBS) minus MODEL
    # Lay out 2*num models + 1
    # columns = Jan and July
    # top row = absolute values from ERA5-L
    # next rows = OBS - MODEL for each model

    nmodels = 2
    nrows = nmodels + 1 # make this dynamic!!!

    fig = plt.figure(figsize=(12, 15))
    fig.suptitle('ERA5-Land Tair Climatology and Model Difference', 
                 y=0.95, # so not tight to top, adjust this ######
                 fontsize=24)
    
    # Add the first row of ERA5-L Tair
    VMIN = 260
    VMAX = 320
    cmap = plt.cm.viridis  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    ## force the first color entry to be grey
    #cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    # define the bins and normalize
    bounds = np.linspace(VMIN, VMAX,13)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    plt.subplot(nrows,2,1)
    plt.title('Jan',  fontsize=24)
    iplt.pcolormesh(era_tair_clim[0], #vmin=VMIN, vmax=VMAX,
                  cmap=cmap, norm=norm)
    #extend='both')
    plt.gca().coastlines()
    # get the current axes' subplot for use later on
    plt1_ax = plt.gca()

    plt.subplot(nrows,2,2)
    plt.title('Jul',  fontsize=24)
    cbar_source = iplt.pcolormesh(era_tair_clim[6], #vmin=VMIN, vmax=VMAX,
                  cmap=cmap, norm=norm)

    plt.gca().coastlines()
    # get the current axes' subplot for use later on
    plt2_ax = plt.gca()

    left, bottom, width, height = plt2_ax.get_position().bounds
    first_plot_left = plt1_ax.get_position().bounds[0]

    # the width of the colorbar should now be simple
    width = left - first_plot_left + width

    # Add axes to the figure, to place the colour bar
    colorbar_axes = fig.add_axes([first_plot_left, bottom  - 0.07,
                                  width, 0.03])

    # Add the colour bar
    cbar = plt.colorbar(cbar_source, colorbar_axes,
                        orientation='horizontal',
                        extend='both')
    
    # Label the colour bar and add ticks
    cbar.set_label('Air Temperature (K)', fontsize=18)
    cbar.ax.tick_params(labelsize=16,length=0)


    #### plots of obs - model

    VMIN = -15
    VMAX = 15
    cmap = plt.cm.RdYlBu_r  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    ## force the first color entry to be grey
    #cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    # define the bins and normalize
    bounds = np.linspace(VMIN, VMAX,11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    model_name = ['Model 1', 'Model 2','Model 3']
    ##### this will end up as loop
    for i in range(nmodels):
        plt.subplot(nrows,2,3+(2*i))
        plt.title(model_name[i], fontsize=24)
        # arbitary difference while jasmin wont find cmip data
        iplt.pcolormesh(era_tair_clim[0]-era_tair_clim[4],
                        cmap=cmap, norm=norm)

        plt1_ax = plt.gca() # this will default to the last row

        plt.subplot(nrows,2,3+(2*i)+1)

        cbar_source = iplt.pcolormesh(era_tair_clim[6]-era_tair_clim[10],
                                    cmap=cmap, norm=norm)
        plt2_ax = plt.gca()

    left, bottom, width, height = plt2_ax.get_position().bounds
    first_plot_left = plt1_ax.get_position().bounds[0]

    # the width of the colorbar should now be simple
    width = left - first_plot_left + width

    # Add axes to the figure, to place the colour bar
    colorbar_axes = fig.add_axes([first_plot_left, bottom  - 0.07,
                                  width, 0.03])

    # Add the colour bar
    cbar = plt.colorbar(cbar_source, colorbar_axes,
                        orientation='horizontal',
                        extend='both')
    
    # Label the colour bar and add ticks
    cbar.set_label('ERA5-Land - Model Air Temperature (K)', fontsize=18)
    cbar.ax.tick_params(labelsize=16,length=0)
        

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, 
                        wspace=0.02, hspace=0.01)

    # datetime_points = Unit.num2date(diffs.coord('time').points,
    #                                str(diffs.coord('time').units),
    #                                diffs.coord('time').units.calendar
    #                            )
    

    

    # import cartopy.crs as ccrs
    # air = xr.tutorial.open_dataset('air_temperature').air
    # ax = plt.axes(projection=ccrs.Orthographic(-80, 35))
    # seasonal.plot.contourf(ax=ax, transform=ccrs.PlateCarree())
    # ax.add_feature(cartopy.feature.BORDERS)
    # ax.coastlines()
    # for i in range(num_of_plots):
    #     print(i)
    #     year = datetime_points[i].year
    #     month = datetime_points[i].month
    #     print(year,month)

    #     fig, ax = plt.subplots(figsize=(20, 15))

    #     iplt.pcolormesh(diffs[i],
    #                     vmin=-10, vmax=10,
    #                     cmap='seismic')
        
    #     plt.gca().coastlines()

    #     plt.title(f'LST-Ta Model {year} {month}')
    #     plt.colorbar()
    
    plt.savefig(f'{outpath}/global_map_tair_clim.png')
    plt.close('all')  # Is this needed?


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
