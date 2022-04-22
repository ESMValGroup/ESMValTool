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

import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np


from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

tab_cols = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
            '#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']

month_list = ['Jan','Feb','Mar','Apr','May','Jun',
              'Jul','Aug','Sep','Oct','Nov','Dec']

def _get_input_cubes(metadata):
    """Load the data files into cubes.
    Based on the hydrology diagnostic.
    Inputs:
    metadata = List of dictionaries made from the preprocessor config
    Outputs:
    inputs = Dictionary of cubes
    ancestors = Dictionary of filename information
    """
    print('################################################')
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
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes

    # loaded data is a nested dictionary
    # KEY1 model ESACCI-LST or something else
    # KEY2 is variable or variable_ensemble
    
    # The Diagnostics

    # CMIP data had 360 day calendar, CCI data has 365 day calendar
    # Assume the loaded data is all the same shape
    print('LOADED DATA:')
    print(loaded_data)

    #### REMEMBER TO APPLY FACTORS TO OBS DATA, DO THIS IN CMORIZER?????
    #### Iris seems to do this automatically for LAI

    # make ensemble mean of area average
    model_means = {}

    icc.add_year(loaded_data['CMUG_WP4_10']['lai'],'time')
    icc.add_month_number(loaded_data['CMUG_WP4_10']['lai'],'time')

    for KEY in loaded_data.keys():
        if KEY == 'CMUG_WP4_10': # add in LAI and Veg KEYS here as they become available
            continue # dont need to do this for CCI

        # loop over ensembles
        ensemble_ts = iris.cube.CubeList()
        for e_number,ENSEMBLE in enumerate(loaded_data[KEY].keys()):

            if ENSEMBLE[0:3] != 'lai':
                continue
               
            this_cube_mean = loaded_data[KEY][ENSEMBLE].collapsed(['latitude','longitude'], iris.analysis.MEAN)

            ensemble_coord = iris.coords.AuxCoord(e_number, standard_name=None, 
                                                  long_name='ensemble_number', 
                                                  var_name=None,
                                                  units='1', bounds=None, 
                                                  attributes=None, coord_system=None)
            this_cube_mean.add_aux_coord(ensemble_coord)
            ensemble_ts.append(this_cube_mean)

        model_means[KEY] = ensemble_ts.merge_cube()
        model_means[KEY] = model_means[KEY].collapsed('ensemble_number', iris.analysis.MEAN)
        icc.add_year(model_means[KEY], 'time')
        icc.add_month_number(model_means[KEY], 'time')
                
    plot_all_members(loaded_data, config)
    plot_season_peaks(loaded_data, model_means, config)


def plot_season_peaks(loaded_data, model_means, config):

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
