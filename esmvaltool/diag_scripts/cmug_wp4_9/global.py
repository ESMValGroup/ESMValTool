"""
CMUG WP4.9
Use this to make plots from GLOBAL data
"""

import logging

import iris
import iris.coord_categorisation as icc
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib as mpl

import matplotlib.pyplot as plt
import numpy as np
import cf_units as Unit

import copy

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

# DOesnt look like this is needed with Iris?
#scale_factor = 0.001 # cant find this in the files, * this by all CCI values to get actual value


# Colours for lines on plots are from
# http://tableaufriction.blogspot.com/2012/11/finally-you-can-use-tableau-data-colors.html
# color blind 10
# with a change in the order.
LINECOLOURS = [
    (200 / 255, 82 / 255, 0),
    (255 / 255, 128 / 255, 14 / 255),
    (0, 107 / 255, 164 / 255),
    (171 / 255, 171 / 255, 171 / 255),
    (95 / 255, 158 / 255, 209 / 255),
    (89 / 255, 89 / 255, 89 / 255),
    (137 / 255, 137 / 255, 137 / 255),
    (162 / 255, 200 / 255, 236 / 255),
    (255 / 255, 188 / 255, 121 / 255),
    (207 / 255, 207 / 255, 207 / 255),
]

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

    # work over france like wp4.8
    FRANCE = [41.,51.5,-5,9] # lat S, lat N, lon W, lon E
    france_lon = iris.Constraint(longitude = lambda cell: FRANCE[2] <= cell <= FRANCE[3])
    france_lat = iris.Constraint(latitude = lambda cell: FRANCE[0] <= cell <= FRANCE[1])

    # for now to do this, long term this data should be loaded via the recipe
    lc_data = iris.load_cube('/work/scratch-nopw/morobking/lc_on_sm/*2003*.nc',
                             'land cover class'
                           )

    # no idea why this wouldn't work as one constraint
    lc_data = lc_data.extract(france_lat)
    lc_data = lc_data.extract(france_lon)

    lc_to_pft = {'NT': [70,71,72,80,81,82],
                 'BT': [50,60,61,62],
                 'G' : [20,130],
                 'SH': [120,121,122,180],
             }
    # lc_data is the dominant lc class from cci lc aggregated to the cci sm grid
    # only need to look at nt bt and g for france

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
        model_mean_lst[MODEL]  = lst_cubelist.merge_cube()

        if enumber > 0:
            model_mean_tair[MODEL] = model_mean_tair[MODEL].collapsed('ensemble_number', iris.analysis.MEAN) 
            model_mean_lst[MODEL]  = model_mean_lst[MODEL].collapsed('ensemble_number', iris.analysis.MEAN) 

    print(f'{model_mean_tair=}')
    print(f'{model_mean_lst=}')

    # make mean lst from day and night
  

    # gah, no idea why the lst data has a messed up longitude coordinate, but this is a fix
    loaded_data['ESACCI_LST_UNCERTS']['tsDay'].coord('longitude').points = np.arange(0,360,1/20.)
    loaded_data['ESACCI_LST_UNCERTS']['tsNight'].coord('longitude').points = np.arange(0,360,1/20.)
    loaded_data['ESACCI_LST_UNCERTS']['tsDay'].coord('longitude').circular=True
    loaded_data['ESACCI_LST_UNCERTS']['tsNight'].coord('longitude').circular=True

    iris.util.promote_aux_coord_to_dim_coord(loaded_data['ESACCI_LST_UNCERTS']['tsDay'], 'longitude')
    iris.util.promote_aux_coord_to_dim_coord(loaded_data['ESACCI_LST_UNCERTS']['tsNight'], 'longitude')

    cci_lst = loaded_data['ESACCI_LST_UNCERTS']['tsDay'] + loaded_data['ESACCI_LST_UNCERTS']['tsNight']
    cci_lst = cci_lst/2
    cci_lst.standard_name = 'surface_temperature'

    era5l_tair = loaded_data['CMUG_WP4_9']['tas']

    # make climatologies
    model_lst_clim = {}
    model_tair_clim = {}
    for MODEL in models:
        model_lst_clim[MODEL]  = model_mean_lst[MODEL].aggregated_by('month_number', iris.analysis.MEAN)
        model_tair_clim[MODEL]  = model_mean_tair[MODEL].aggregated_by('month_number', iris.analysis.MEAN)

    cci_lst_clim    = cci_lst.aggregated_by('month_number', iris.analysis.MEAN)
    era5l_tair_clim = era5l_tair.aggregated_by('month_number', iris.analysis.MEAN)

    #### regrid tair to model grid
    
    # Need to do this for each model :(
    cci_lst_clim_regrided = {}
    era5l_tair_clim_regrided = {}
    cci_lst_clim.coord(axis='x').attributes = None
    era5l_tair_clim.coord(axis='x').attributes = None
    cci_lst_clim.coord(axis='y').attributes = None
    era5l_tair_clim.coord(axis='y').attributes = None
    cci_lst_clim.coord(axis='x').bounds = None
    era5l_tair_clim.coord(axis='x').bounds = None
    cci_lst_clim.coord(axis='y').bounds = None
    era5l_tair_clim.coord(axis='y').bounds = None

    for i,MODEL in enumerate(models):
    # There seems to be a coord issue when regridding that is model dependant. 
    # This gets rid of all the gotcha and seems to make everything work
        cci_lst_clim.coord(axis='x').coord_system = model_mean_lst[MODEL].coord(axis='x').coord_system
        cci_lst_clim.coord(axis='y').coord_system = model_mean_lst[MODEL].coord(axis='y').coord_system
        era5l_tair_clim.coord(axis='x').coord_system = model_mean_lst[MODEL].coord(axis='x').coord_system
        era5l_tair_clim.coord(axis='y').coord_system = model_mean_lst[MODEL].coord(axis='y').coord_system

        model_lst_clim[MODEL].coord(axis='x').attributes = None
        model_tair_clim[MODEL].coord(axis='x').attributes = None

        model_lst_clim[MODEL].coord(axis='y').attributes = None
        model_tair_clim[MODEL].coord(axis='y').attributes = None

        model_lst_clim[MODEL].coord(axis='x').bounds = None
        model_tair_clim[MODEL].coord(axis='x').bounds = None
 
        model_lst_clim[MODEL].coord(axis='y').bounds = None
        model_tair_clim[MODEL].coord(axis='y').bounds = None

        cci_lst_clim_regrided[MODEL] =  cci_lst_clim.regrid(model_lst_clim[MODEL], iris.analysis.Linear())
        era5l_tair_clim_regrided[MODEL] = era5l_tair_clim.regrid(model_tair_clim[MODEL], iris.analysis.Linear())

    # LST-Tair Cube work
    # this will need a model grid dependant thing doing!
    #diffs_clim_obs = cci_lst_clim_regrided - era5l_tair_clim_regrided
    #print(f'{diffs_clim_obs=}')
    
    # plots
    #make_plot_global_clim_maps(model_lst, model_ta, cci_lst, era5l_tair, config)
    # make this so either lst or tair passed in, will need vmin vmax as parameters #####
    print('Start making plots')
    # Maps of climatology
    # make_plot_global_clim_maps(obs_data, model_data, 'LST or Tair', vmin_obs, vmax_obs, vnum_obs, vmin_diff, vmax_diff, vnum_obs,config)
    make_plot_global_clim_maps(cci_lst_clim_regrided, model_lst_clim, 'LST',270,330,15 ,-15,15, 11, config)
    make_plot_global_clim_maps(era5l_tair_clim_regrided, model_tair_clim, 'Tair',270,320,15 ,-15,15, 11, config)
    print('plot made')

    # make LST-Tair for OBS and Model
    #
    # cube work here
    obs_diff_clim = {}
    model_diff_clim = {}
    for MODEL in models:
        obs_diff_clim[MODEL] =  cci_lst_clim_regrided[MODEL] - era5l_tair_clim_regrided[MODEL]
        model_diff_clim[MODEL] = model_lst_clim[MODEL] - model_tair_clim[MODEL] ################################################################
    
    make_plot_global_clim_maps(obs_diff_clim, model_diff_clim, 'Diff',-15,15,21 ,-15,15, 11, config)

    # NEXT do France LC plots
    # look in WP4.8 plots_for_reports.py to get code to mask lst,tair,diff as needed on lc types

    # Do this just for LST-Tair diffrence
    # make climatologies per biome
    # regrid to LC grid
    # mask on pfts
    # plot result for all models and obs
    obs_on_lc_clim = {}
    model_on_lc_clim = {}
    outpath = config['plot_dir']

    for MODEL in models:
        obs_on_lc_clim[MODEL] = obs_diff_clim[MODEL].regrid(lc_data, iris.analysis.Linear())
        model_on_lc_clim[MODEL] = model_diff_clim[MODEL].regrid(lc_data, iris.analysis.Linear())
        qplt.pcolormesh(obs_on_lc_clim[MODEL][0])
        plt.gca().coastlines()
        plt.savefig(f'{outpath}/test_{MODEL}_clim.png')
        plt.close()

    print(f'{obs_on_lc_clim=}')
    print(f'{model_on_lc_clim=}')
    print(f'{lc_data=}')

    obs_clim_lc = {}
    model_clim_lc = {}
    
    for BIOME in ['NT', 'BT']:#lc_to_pft.keys():
        print(BIOME)
        obs_clim_lc[BIOME] = {}
        model_clim_lc[BIOME] = {}
        biome_mask = np.ma.masked_outside(lc_data[0].data,
                                          np.min(lc_to_pft[BIOME]),
                                          np.max(lc_to_pft[BIOME])
                                          )
        print(biome_mask.shape)
        new_mask = np.ones((12,biome_mask.shape[0], biome_mask.shape[1]))
        for i in range(12):
            new_mask[i] = biome_mask.mask
        
        plt.pcolormesh(new_mask[0])
        # a check on mask
        plt.savefig(f'{outpath}/mask_{BIOME}.png')
        plt.close()
        for MODEL in models:
            # appears deep copy is needed, not sure why, guess cubes considered nested objects?
            X = copy.deepcopy(obs_on_lc_clim[MODEL])
            X.data = np.ma.masked_array(obs_on_lc_clim[MODEL].data, mask=new_mask)
            obs_clim_lc[BIOME][MODEL] = X
            
            Y = copy.deepcopy(model_on_lc_clim[MODEL])
            Y.data = np.ma.masked_array(model_on_lc_clim[MODEL].data, mask=new_mask)
            model_clim_lc[BIOME][MODEL] = Y

    plt.close()

    # for BIOME in ['NT','BT']:
    #     for MODEL in models:
    #         qplt.pcolormesh(obs_clim_lc[BIOME][MODEL][0])
    #         plt.gca().coastlines()
    #         plt.savefig(f'{outpath}/mask_{MODEL}_obs_clim_{BIOME}.png')
    #         plt.close()

    #         qplt.pcolormesh(model_clim_lc[BIOME][MODEL][0])
    #         plt.gca().coastlines()
    #         plt.savefig(f'{outpath}/mask_{MODEL}_model_clim_{BIOME}.png')
    #         plt.close()

    plot_biome_timeseries(obs_clim_lc, model_clim_lc, True, config)


def plot_biome_timeseries(obs, model, clim, config):
    # clim = True for 12 month, False for a whole length plot

    for BIOME in obs.keys():
        
        plt.figure(figsize=(10,15))

        obs_count = 0
        for MODEL in obs[BIOME].keys():
            print(MODEL)
            if obs_count == 0:
                obs_ave = obs[BIOME][MODEL].collapsed(['latitude','longitude'],
                                                      iris.analysis.MEAN)

                print(obs[BIOME][MODEL].data)
                print(obs_ave.data)

                plt.plot(obs_ave.data, color='black', label='OBS')
                obs_count = 1

            this_model_ave = model[BIOME][MODEL].collapsed(['latitude','longitude'],
                                                           iris.analysis.MEAN)
            print(this_model_ave)
            print(this_model_ave.data)
            
            plt.plot(this_model_ave.data, label=MODEL)

        plt.legend()
        outpath = config['plot_dir']
        plt.savefig(f'{outpath}/timeseries_{BIOME}_clim.png')
        plt.close('all')  # Is this needed?


            
            

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


# sort this to work with diffs too!
def make_plot_global_clim_maps(obs_data, model_data, 
                               TYPE, # LST or Tair
                               vmin_obs, vmax_obs, vnum_obs, 
                               vmin_diff, vmax_diff, vnum_diff,
                               config):
    """Create and save the output figure.
    Figure of airT from ERA5 and differences from airT in CMIP6 models (historic climatology – 2003-13 obs and model years)
    Figure of LST from CCI LST and differences from LST in CMIP6 models (historic climatology – 2003-2013 obs and model years)
    These could be for individual months or whole year TBC
    ###### UPDATE THIS ######
    Inputs:
   
    config = The config dictionary from the preprocessor
    Outputs:
    Saved figure
    """
    outpath = config['plot_dir']

    cbar_obs_text = {'LST': 'LST (K)',
                     'Tair': 'Air Temperature (K)',
                     'Diff': 'LST-Tair (K)'
                     }
    cmaps = {'LST':  plt.cm.viridis,
             'Tair': plt.cm.viridis,
             'Diff': plt.cm.bwr
             }

    # First figure
    # Tair from ERA5-land and ERA(OBS) minus MODEL
    # Lay out 2*num models + 1
    # columns = Jan and July
    # top row = absolute values from ERA5-L
    # next rows = OBS - MODEL for each model
    print(f'In the plotting bit {TYPE}')
    nmodels = len(list(model_data.keys()))
    nrows = nmodels + 1 

    fig = plt.figure(figsize=(12, 30))
    fig.suptitle(f'ERA5-Land {TYPE} Climatology and Model Difference', 
                 y=0.95, # so not tight to top, adjust this ######
                 fontsize=24)
    
    # Add the first row of ERA5-L Tair

    cmap = cmaps[TYPE]
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    ## force the first color entry to be grey
    #cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    # define the bins and normalize
    bounds = np.linspace(vmin_obs, vmax_obs, vnum_obs)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    plt.subplot(nrows,2,1)
    plt.title('Jan',  fontsize=24)
    iplt.pcolormesh(obs_data['CESM2'][0],
                  cmap=cmap, norm=norm)
    plt.gca().coastlines()
    # get the current axes' subplot for use later on
    plt1_ax = plt.gca()

    plt.subplot(nrows,2,2)
    plt.title('Jul',  fontsize=24)
    cbar_source = iplt.pcolormesh(obs_data['CESM2'][6],
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
    cbar.set_label(cbar_obs_text[TYPE], fontsize=18)
    cbar.ax.tick_params(labelsize=16,length=0)


    #### plots of obs - model

    cmap = plt.cm.RdYlBu_r  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    ## force the first color entry to be grey
    #cmaplist[0] = (.5, .5, .5, 1.0)
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    # define the bins and normalize
    bounds = np.linspace(vmin_diff, vmax_diff, vnum_diff)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    model_names= list(model_data.keys())
    ##### this will end up as loop
    for i,MODEL in enumerate(model_names):
        print(MODEL)
        plt.subplot(nrows,2,3+(2*i))
        plt.title(model_names[i], fontsize=24)
        # arbitary difference while jasmin wont find cmip data
        iplt.pcolormesh(obs_data[MODEL][0]-model_data[MODEL][0],
                        cmap=cmap, norm=norm)
        plt.gca().coastlines()
        plt1_ax = plt.gca() # this will default to the last row

        plt.subplot(nrows,2,3+(2*i)+1)

        cbar_source = iplt.pcolormesh(obs_data[MODEL][6]-model_data[MODEL][6],
                                    cmap=cmap, norm=norm)
        plt.gca().coastlines()

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
    cbar.set_label('OBS - Model Difference (K)', fontsize=18)
    cbar.ax.tick_params(labelsize=16,length=0)
        

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, 
                        wspace=0.01, hspace=0.002)

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
    
    plt.savefig(f'{outpath}/global_map_{TYPE}_clim.png')
    plt.close('all')  # Is this needed?


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
