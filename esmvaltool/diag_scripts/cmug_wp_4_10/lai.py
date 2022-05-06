"""
ESMValTool diagnostic for ESA CCI LST data.
The code uses the all time average monthly data.
The ouptput is a timeseries plot of the mean differnce of
CCI LST to model ensemble average, with the ensemble spread
represented by a standard deviation either side of the mean.
"""

import logging

import iris
import iris.coord_categorisation as icc
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np


from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

tab_cols = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
            '#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']

month_list = np.array(['Jan','Feb','Mar','Apr','May','Jun',
                       'Jul','Aug','Sep','Oct','Nov','Dec'])


def _get_input_cubes(metadata):
    """Load the data files into cubes.
    Based on the hydrology diagnostic.
    Inputs:
    metadata = List of dictionaries made from the preprocessor config
    Outputs:
    inputs = Dictionary of cubes
    ancestors = Dictionary of filename information
    """
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
        
        if 'CMIP5' in attributes['alias']:
            data_type = 'CMIP5'
        elif 'OBS' in attributes['alias']:
            data_type = 'OBS'
        else:
            data_type = 'CMIP6' # this way meand CMIP5 doesnt get counted twice

    return inputs, ancestors, data_type

# def _get_input_cubes(metadata):
#     """Load the data files into cubes.
#     Based on the hydrology diagnostic.
#     Inputs:
#     metadata = List of dictionaries made from the preprocessor config
#     Outputs:
#     inputs = Dictionary of cubes
#     ancestors = Dictionary of filename information
#     """
#     print('################################################')
#     inputs = {}
#     ancestors = {}
#     print(metadata)
#     for attributes in metadata:
#         print(attributes)
#         short_name = attributes['short_name']
#         filename = attributes['filename']
#         logger.info("Loading variable %s", short_name)
#         cube = iris.load_cube(filename)
#         cube.attributes.clear()
        
#         try:
#             key_name = f"{short_name}_{attributes['ensemble']}"
#         except:
#             key_name = short_name

#         inputs[key_name] = cube
#         ancestors[key_name] = [filename]
        
#         print(inputs)
#         print(ancestors)

#     return inputs, ancestors

def _get_provenance_record(attributes, ancestor_files):
    """Create the provenance record dictionary.
    Inputs:
    attributes = dictionary of ensembles/models used, the region bounds
                 and years of data used.
    ancestor_files = list of data files used by the diagnostic.
    Outputs:
    record = dictionary of provenance records.
    """
    ### THIS ALL NEEDS CHANGING!!!
    caption = "CMUG WP4.10"

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
        cubes, ancestors, data_type = _get_input_cubes(metadata)
        for KEY in cubes.keys():
            cubes[KEY].coord('longitude').circular = True

        loaded_data[f'{data_type}_{dataset}'] = cubes

    # The Diagnostics
    
    #### REMEMBER TO APPLY FACTORS TO OBS DATA, DO THIS IN CMORIZER?????
    #### Iris seems to do this automatically for LAI

    print(loaded_data)
    loaded_data['OBS_CMUG_WP4_10']['lai'].coord('latitude').coord_system = None
    loaded_data['OBS_CMUG_WP4_10']['lai'].coord('longitude').circular = True

    loaded_data['OBS_CMUG_WP4_10']['lai'].coord('longitude').coord_system = None

    print( loaded_data['OBS_CMUG_WP4_10']['lai'].coord('longitude'))


    icc.add_year(loaded_data['OBS_CMUG_WP4_10']['lai'],'time')
    icc.add_month_number(loaded_data['OBS_CMUG_WP4_10']['lai'],'time')


    # make ensemble average of models
    model_means = {}
    for KEY in loaded_data.keys():
        if KEY == 'OBS_CMUG_WP4_10': # add in LAI and Veg KEYS here as they become available
            continue # dont need to do this for CCI

        # loop over ensembles
        ensemble_ts = iris.cube.CubeList()
        for e_number,ENSEMBLE in enumerate(loaded_data[KEY].keys()):

            if ENSEMBLE[0:3] != 'lai':
                continue
               
            this_cube = loaded_data[KEY][ENSEMBLE]#.collapsed(['latitude','longitude'], iris.analysis.MEAN)

            ensemble_coord = iris.coords.AuxCoord(e_number, standard_name=None, 
                                                  long_name='ensemble_number', 
                                                  var_name=None,
                                                  units='1', bounds=None, 
                                                  attributes=None, coord_system=None)
            this_cube.add_aux_coord(ensemble_coord)
            ensemble_ts.append(this_cube)

        model_means[KEY] = ensemble_ts.merge_cube()
        model_means[KEY] = model_means[KEY].collapsed('ensemble_number', iris.analysis.MEAN)
        icc.add_year(model_means[KEY], 'time')
        icc.add_month_number(model_means[KEY], 'time')
            


    # regrid obs to each model grid - they seem to be differnt
    lai_obs_regridded = {}
    for KEY in loaded_data.keys():
        if KEY == 'OBS_CMUG_WP4_10': continue

        for item in loaded_data[KEY].keys():
            print(item)
            ensemble = item
            break
        
        lai_regrid = loaded_data['OBS_CMUG_WP4_10']['lai'].regrid(
            loaded_data[KEY][ensemble], iris.analysis.Linear())

        lai_obs_regridded[KEY] = lai_regrid

   

    print("££££££££££££££££££££££££££££££")
    print(lai_obs_regridded)

    # calculate gridbox peaks for maps
    # OBS peak
    ARGMAX = iris.analysis.Aggregator("argmax", np.argmax, units_func=lambda units: 1)
    obs_peak_regridded = {}
    for KEY in lai_obs_regridded.keys():
        this_cubelist = iris.cube.CubeList()
        for YEAR in np.unique(lai_obs_regridded[KEY].coord('year').points):
            cube = lai_obs_regridded[KEY].extract(iris.Constraint(year=YEAR))
            cube = cube.collapsed('time',ARGMAX)
            
            cube.coord('year').bounds=None
            
            cube.remove_coord('time')
            cube.remove_coord('month_number')
            cube = iris.util.new_axis(cube, scalar_coord='year')
            this_cubelist.append(cube)

        obs_peak_regridded[KEY] = this_cubelist.concatenate_cube()
        obs_peak_regridded[KEY] = obs_peak_regridded[KEY].collapsed('year', iris.analysis.MEAN)

    # MODEL peak
    model_peak = {}
    for KEY in model_means.keys():
        print(KEY)
        this_cubelist = iris.cube.CubeList()
        for YEAR in np.unique(model_means[KEY].coord('year').points):
            print(YEAR)
            cube = model_means[KEY].extract(iris.Constraint(year=YEAR))
            cube = cube.collapsed('time',ARGMAX)
            
            cube.coord('year').bounds=None
            
            cube.remove_coord('time')
            cube.remove_coord('month_number')
            cube = iris.util.new_axis(cube, scalar_coord='year')
            print(cube.coord('year'))
            this_cubelist.append(cube)

        model_peak[KEY] = this_cubelist.concatenate_cube()
        model_peak[KEY] = model_peak[KEY].collapsed('year', iris.analysis.MEAN)


    plot_peak_diff_map(loaded_data, model_peak, obs_peak_regridded, config)
    #print(0/0)

    
    #plot_all_members(loaded_data, config)
    #plot_season_peaks(loaded_data, model_means, config)


def plot_peak_diff_map(loaded_data, model_peak, obs_peak_regridded, config):

    num_plot = len(model_peak.keys())+1

#    fig, ax = plt.subplots(nrows=num_plot, ncols=1, figsize=(20, 15))
    #fig = plt.figure(figsize=(20, 15))
    #plt.axes(projection=ccrs.NorthPolarStereo())
    # first plot is the OBS LAI
    #plt.subplot(num_plot,1,1)
   #
    
    im0 = iplt.pcolormesh(obs_peak_regridded['CMIP6_UKESM1-0-LL'],
                          cmap='Paired',
                          vmin=1,
                          vmax=12)
    plt.gca().coastlines()
    plt.xlim((-13,38))
    #plt.ylim((34,73))

    plt.colorbar(orientation='vertical')


    # for i,KEY in enumerate(model_peak.keys()):
    #     this_diff =  obs_peak_regridded[KEY] - model_peak[KEY]

    #     plt.subplot(num_plot,1,i+2)
    #     iplt.pcolormesh(this_diff,
    #                     cmap='PiYG',
    #                     vmin=-6,
    #                     vmax=6)
        
    #     plt.xlim((-13,38))
    #     plt.ylim((34,73))
    #     plt.gca().coastlines()
    #     plt.colorbar(orientation='vertical')
    #     plt.title(f'{KEY}',fontsize=18)


    # fig.suptitle('LAI from OBS and OBS-MODEL', fontsize=24)
    #qplt.pcolormesh(obs_peak_regridded['CMIP6_UKESM1-0-LL'])
    outpath =  config['plot_dir']
    plt.savefig(f'{outpath}/lai_peak_map.png')
    print(f'{outpath}/lai_peak_map.png')
    plt.close('all')  # Is this needed?


def plot_season_peaks(loaded_data, model_means, config):

    peaks_all = {}
    print('PEAKS ALL')
    for KEY in loaded_data.keys():
        if KEY == 'CMUG_WP4_10': continue
        peaks_all[KEY] = []
        print(KEY)
        for e_number,ENSEMBLE in enumerate(loaded_data[KEY].keys()):
            if ENSEMBLE[0:3] != 'lai': continue
            print(ENSEMBLE)
            this_result = []
            try:
                icc.add_year(loaded_data[KEY][ENSEMBLE],'time')
            except:
                pass

            for YEAR in np.unique(loaded_data[KEY][ENSEMBLE].coord('year').points):
                print(YEAR)
                cube = loaded_data[KEY][ENSEMBLE].extract(iris.Constraint(year=YEAR))
                cube_values= (cube.collapsed(['latitude','longitude'],iris.analysis.MEAN)).data                
                print(cube_values)
                this_result.append(np.argmax(cube_values))
            print(this_result)
            peaks_all[KEY].append(this_result)

            print('$$$$$$')
    print(peaks_all)
            
        
    peaks = {}
    for i, MODEL in enumerate(model_means.keys()):
        peaks[MODEL]= []
        for YEAR in np.unique(model_means[MODEL].coord('year').points):
            cube = model_means[MODEL].extract(iris.Constraint(year=YEAR))
            cube_values= cube.data
            peaks[MODEL].append(np.argmax(cube_values))

    obs_peak = []
    obs_mean = loaded_data['CMUG_WP4_10']['lai'].collapsed(['latitude','longitude'], iris.analysis.MEAN)
    for YEAR in np.unique(obs_mean.coord('year').points):
        cube = obs_mean.extract(iris.Constraint(year=YEAR))
        cube_values= cube.data
        obs_peak.append(np.argmax(cube_values))

    # basic plot with just model mean peak
    fig, ax = plt.subplots(figsize=(20, 15))

    for i,KEY in enumerate(peaks.keys()):
        ax.plot(peaks[KEY],'o-',
                label=KEY,
                color=tab_cols[i])

    ax.plot(obs_peak,'o-', label = 'OBS',
             color='black')

    year_list = np.unique(obs_mean.coord('year').points)
    ax.set_xticks(np.arange(0,len(year_list)))
    ax.set_xticklabels(year_list, fontsize=18)
    ax.set_xlim((-1,len(year_list)+1))

    ax.set_yticks(np.arange(0,12))
    ax.set_yticklabels(month_list, fontsize=18)
    ax.set_ylim((-1,12))

    ax.grid()

    ax.legend(loc='upper right')

    outpath =  config['plot_dir']
    plt.savefig(f'{outpath}/lai_peak.png')
    plt.close('all')  # Is this needed?

    
    # plot with model range of peaks
    # years on y axis, x axis is models
    # plot shaded range for each model like a violin plot, dot for obs
    # aspect will need changing
    fig, ax = plt.subplots(figsize=(20, 15))

    x_labels = []
    x_ticks  = np.array([])
    year_list = np.unique(obs_mean.coord('year').points)
    for i,KEY in enumerate(peaks_all.keys()):
        X = i*20
        x_labels.append(KEY)
        
        this_array = np.array(peaks_all[KEY])

        ax.fill_betweenx(year_list,
                         X+np.min(peaks_all[KEY],axis=0),
                         X+np.max(peaks_all[KEY],axis=0),
                         color=tab_cols[i], alpha=0.5
                         )
        ax.plot(X+np.mean(peaks_all[KEY],axis=0), year_list,
                color=tab_cols[i],
                linewidth=3
                )

        ax.plot(X+np.array(obs_peak), year_list,
                'o',
                color='black')

        ax.text(X,2016,KEY,
                fontsize=20
                )

    for j in range(i+1):
        x_ticks = np.append(x_ticks, np.arange(j*20,j*20 + 12))
    
    x_labels = np.tile(month_list,i+1)


    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=18, rotation=45)
    
    ax.set_yticks(year_list)
    ax.set_yticklabels(year_list, fontsize=18)

    ax.grid()


    plt.savefig(f'{outpath}/lai_yearly_peak.png')
    plt.close('all')  # Is this needed?


   
def plot_all_members(loaded_data, config):
    
    fig, ax = plt.subplots(figsize=(20, 15))

    for i,KEY in enumerate(loaded_data.keys()):

        if KEY == 'CMUG_WP4_10':
            COL = 'black'
            LW=1
        else:
            COL = tab_cols[i]
            LW=3

        for e_number,ENSEMBLE in enumerate(loaded_data[KEY].keys()):
            this_cube_mean = loaded_data[KEY][ENSEMBLE].collapsed(['latitude','longitude'], iris.analysis.MEAN)
        
            ax.plot(this_cube_mean.data, label = f"{KEY} {ENSEMBLE}",
                      color=COL)

    x_tick_list = []
    time_list = this_cube_mean.coord('time').units.num2date(this_cube_mean.coord('time').points)
    for item in time_list:
        if item.month == 1:
            x_tick_list.append(item.strftime('%Y %b'))
        elif item.month == 7:
            x_tick_list.append(item.strftime('%b'))
        else:
            x_tick_list.append('')

    ax.set_xticks(range(len(this_cube_mean.data)))
    ax.set_xticklabels(x_tick_list, fontsize=18, rotation=45)

    ax.set_yticks(np.arange(0,6))
    ax.set_yticklabels(np.arange(0,6), fontsize=18)
    ax.set_ylim((0,6))

    ax.set_xlabel('Date', fontsize=20)
    ax.set_ylabel('LAI', fontsize=20)


    lons = loaded_data['CMUG_WP4_10']['lai'].coord('longitude').bounds
    lats = loaded_data['CMUG_WP4_10']['lai'].coord('latitude').bounds

    ax.set_title('Area: lon %s lat %s' % (lons[0], lats[0]), fontsize=22)

    fig.suptitle('ESA LAI and CMIP6 LAI', fontsize=24)
    

    ax.grid()

    outpath =  config['plot_dir']
    plt.legend(loc='upper right')
    plt.savefig(f'{outpath}/lai.png')
    plt.close('all')  # Is this needed?


    # _make_plots(cci_lst, total_uncert, model_lst, model_std, ensemble_ts, config)

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


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
