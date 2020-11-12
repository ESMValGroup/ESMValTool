# ESA CCI LST Diagnostic

import logging

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__) # from OC example, dont know what this does!

def get_input_cubes(metadata):
    # This is from 
    # https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/hydrology/lisflood.py

    inputs = {}
    ancestors = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        inputs[short_name] = cube
        ancestors[short_name] = [filename]

    return inputs, ancestors

def make_plots(lst_diff_data,lst_diff_data_low,lst_diff_data_high, config, input_metadata):
    # Make a timeseries plot of the difference OBS-MODEL

    fig,ax = plt.subplots(figsize=(20,15))

    ax.plot(lst_diff_data.data, color='black', linewidth=4)
    ax.plot(lst_diff_data_low.data,'--', color='blue', linewidth=3)
    ax.plot(lst_diff_data_high.data,'--', color='blue', linewidth=3)
    ax.fill_between(range(len(lst_diff_data.data)), lst_diff_data_low.data,lst_diff_data_high.data,
                    color='blue',alpha=0.25)

    # make X ticks
    x_tick_list = []
    time_list = lst_diff_data.coord('time').units.num2date(lst_diff_data.coord('time').points)
    for item in time_list:
        if item.month == 1:
            x_tick_list.append(item.strftime('%Y %b'))
        elif item.month == 7:
            x_tick_list.append(item.strftime('%b'))
        else:
            x_tick_list.append('')

    ax.set_xticks(range(len(lst_diff_data.data)))
    ax.set_xticklabels(x_tick_list, fontsize=18, rotation = 45)

    # make Y ticks
    Y_lower = np.floor(lst_diff_data_low.data.min())
    Y_upper = np.ceil(lst_diff_data_high.data.max())
    ax.set_yticks(np.arange(Y_lower,Y_upper+0.1,2))
    ax.set_yticklabels(np.arange(Y_lower,Y_upper+0.1,2), fontsize=18)
    ax.set_ylim((Y_lower-0.1,Y_upper+0.1))

    ax.set_xlabel('Date', fontsize = 20)
    ax.set_ylabel('Difference / K', fontsize = 20)

    ax.grid()

    lons = lst_diff_data.coord('longitude').bounds
    lats = lst_diff_data.coord('latitude').bounds

    ax.set_title('Area: lon %s lat %s' % (lons[0], lats[0]), fontsize=22)

    fig.suptitle('ESACCI LST - CMIP6 Historical Ensemble Mean', fontsize=24)

    plt.savefig('%s/timeseries.png' % config['plot_dir'])
    plt.close('all') # Is this needed?

    return None

def get_provenance_record(attributes, ancestor_files):
    
    caption = ("Timeseries of ESA CCI LST difference to mean of model ensembles, calculated over region bounded by latitude {lat_south} to {lat_north}, longitude {lon_west} to {lon_east} and for model/ensembles {ensembles}. Shown for years {start_year} to {end_year}."
               .format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['mean','stddev'],
        'domains': ['reg'],
        'plot_types': ['times'],
        'authors': ['king_robert'],
        #'references': [],
        'ancestors': ancestor_files
    }

    return record

def diagnostic(config):
    
    # this loading function is based on the hydrology py above
    input_metadata = config['input_data'].values()
    logger.info(config)
    loaded_data = {}
    ancestor_list = []
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = get_input_cubes(metadata)
        loaded_data[dataset] = cubes
        ancestor_list.append(ancestors['ts'][0])

    # loaded data is a nested dictionary
    # KEY1 model ESACCI-LST or something else
    # KEY2 is ts, the surface temperature
    # ie loaded_data['ESACCI-LST']['ts'] is the CCI cube
    #    loaded_data['MultiModelMean']['ts'] is CMIP6 data, emsemble means, see preprocessor

    ##### The Diagnostic uses CCI - MODEL

    # CMIP data had 360 day calendar, CCI data has 365 day calendar
    # Assume the loaded data is all the same shape
    loaded_data['MultiModelMean']['ts'].remove_coord('time')
    loaded_data['MultiModelMean']['ts'].add_dim_coord(loaded_data['ESACCI-LST']['ts'].coord('time'),0)
    loaded_data['MultiModelStd']['ts'].remove_coord('time')
    loaded_data['MultiModelStd']['ts'].add_dim_coord(loaded_data['ESACCI-LST']['ts'].coord('time'),0)

    # Make a cube of the LST difference, and with +/- std of model variation
    lst_diff_cube = loaded_data['ESACCI-LST']['ts'] - loaded_data['MultiModelMean']['ts']
    lst_diff_cube_low = loaded_data['ESACCI-LST']['ts'] - (loaded_data['MultiModelMean']['ts']+loaded_data['MultiModelStd']['ts'])
    lst_diff_cube_high = loaded_data['ESACCI-LST']['ts'] - (loaded_data['MultiModelMean']['ts']-loaded_data['MultiModelStd']['ts'])

    # Plotting
    make_plots(lst_diff_cube,lst_diff_cube_low,lst_diff_cube_high, config, input_metadata)

    # Provenance
    # Get this information form the data cubes
    data_attributes = {}
    data_attributes['start_year'] = lst_diff_cube.coord('time').units.num2date(lst_diff_cube.coord('time').points)[0].year
    data_attributes['end_year']   = lst_diff_cube.coord('time').units.num2date(lst_diff_cube.coord('time').points)[-1].year
    data_attributes['lat_south']  = lst_diff_cube.coord('latitude').bounds[0][0]
    data_attributes['lat_north']  = lst_diff_cube.coord('latitude').bounds[0][1]
    data_attributes['lon_west']   = lst_diff_cube.coord('longitude').bounds[0][0]
    data_attributes['lon_east']   = lst_diff_cube.coord('longitude').bounds[0][1]
    data_attributes['ensembles'] = ''

    for ITEM in input_metadata:
        if 'ESACCI' in ITEM['alias'] or 'MultiModel' in ITEM['alias'] or 'OBS' in ITEM['alias']:
            continue
        data_attributes['ensembles'] += "%s " % ITEM['alias']

    record = get_provenance_record(data_attributes,ancestor_list)
    for file in ['%s/timeseries.png' % config['plot_dir']]:
        with ProvenanceLogger(config) as provenance_logger:
            provenance_logger.log(file, record)

    return None

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        diagnostic(config)
