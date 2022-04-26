"""
CMUG WP4.9
Use this to make plots from GLOBAL data
"""

import logging

import iris
import iris.plot as iplt
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
    print(loaded_data)

    
    ### There will be some cube manipulation todo

    

        
    ### PLOT is CCI LST with bars of uncertainty
    ####     with shaded MODEL MEAN +/- std
    # Plotting
    print(list(loaded_data.keys()))
    model = list(loaded_data.keys())[0]
    print(model)

    model_lst = loaded_data[model]['ts']
    model_ta  = loaded_data[model]['tas']

    _make_plot_global_map(model_lst, model_ta, config)

    # more plotting functions go here

   #  print(0/0)
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



def _make_plot_global_map(model_lst, model_ta, config):
    """Create and save the output figure.
    Plot maps of LST-Ta from model
    Inputs:
   
    config = The config dictionary from the preprocessor
    Outputs:
    Saved figure
    """


    num_of_plots = len(model_lst.coord('time').points)

    diffs = model_lst - model_ta
    datetime_points = Unit.num2date(diffs.coord('time').points,
                                   str(diffs.coord('time').units),
                                   diffs.coord('time').units.calendar
                               )
    
    print(diffs)
    
    outpath = config['plot_dir']

    for i in range(num_of_plots):
        print(i)
        year = datetime_points[i].year
        month = datetime_points[i].month
        print(year,month)

        fig, ax = plt.subplots(figsize=(20, 15))

        iplt.pcolormesh(diffs[i],
                        vmin=-10, vmax=10,
                        cmap='seismic')
        
        plt.gca().coastlines()

        plt.title(f'LST-Ta Model {year} {month}')
        plt.colorbar()
    
        plt.savefig(f'{outpath}/global_map_{year}_{month}.png')
        plt.close('all')  # Is this needed?


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
