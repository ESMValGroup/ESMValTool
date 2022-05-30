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
import cf_units as unit

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np
import datetime

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

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

LINE_STYLES = ['-','--',':']

MONTH_LIST = ['Jan','Feb','Mar','Apr','May','Jun',
              'Jul','Aug','Sep','Oct','Nov','Dec']


ARGMAX = iris.analysis.Aggregator("argmax", np.argmax, units_func=lambda units: 1)


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

    print(f'{loaded_data=}')

    lai = {}
    ts = {}
    tas = {}
    sm = {}

    # Assume single ensemble for each MODEL
    for KEY in loaded_data.keys(): # this is the model
        print(KEY)

        for ITEM in loaded_data[KEY].keys():
            print(ITEM) # this is the variable

            icc.add_year(loaded_data[KEY][ITEM], 'time')
            icc.add_month_number(loaded_data[KEY][ITEM], 'time')
            
            loaded_data[KEY][ITEM].coord('latitude').bounds = None
            loaded_data[KEY][ITEM].coord('longitude').bounds = None

            if 'ts' in ITEM:
                ts[KEY] = loaded_data[KEY][ITEM]

            if 'tas' in ITEM:
                tas[KEY] = loaded_data[KEY][ITEM]

            if 'lai' in ITEM:
                lai[KEY] = loaded_data[KEY][ITEM]
                
            if 'mrso' in ITEM:
                sm[KEY] = loaded_data[KEY][ITEM]

    print(ts)
    print(tas)
    print(lai)
    print(sm)

    diff = {}
    for KEY in ts.keys():
        diff[KEY] = ts[KEY] - tas[KEY]
    
    print(diff)

    # now plot
    ts_monthly = {}
    tas_monthly = {}
    lai_monthly = {}
    sm_monthly = {}
    diff_monthly = {}
    for KEY in ts.keys():
        ts_monthly[KEY] = ts[KEY].aggregated_by('month_number', iris.analysis.MEAN)
        tas_monthly[KEY] = tas[KEY].aggregated_by('month_number', iris.analysis.MEAN)
        lai_monthly[KEY] = lai[KEY].aggregated_by('month_number', iris.analysis.MEAN)
        sm_monthly[KEY] = sm[KEY].aggregated_by('month_number', iris.analysis.MEAN)
        diff_monthly[KEY] = diff[KEY].aggregated_by('month_number', iris.analysis.MEAN)


    outpath = config['plot_dir']

    cmap = plt.cm.YlGn
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(0, 6, 13)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # figure x=sm y=diff
    for KEY in ts.keys():
        plt.figure(figsize=(15,15))
        im = qplt.scatter(sm_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN),
                     diff_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN),
                     #edgecolors=sm_monthly[KEY].coord('month_number').points,
                     c=lai_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN).data,
                     cmap=cmap, norm=norm,
                     s=200,marker='o',
                     zorder=3
                     )
        plt.colorbar(im,label='LAI')

        im=qplt.scatter(sm_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN),
                     diff_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN),
                     c=sm_monthly[KEY].coord('month_number').points,
                     #c=lai_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN).data,
                     cmap='hsv',
                     s=300,marker='o',
                     zorder=1
                     )
        plt.colorbar(im,label='month number')

        iplt.plot(sm_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN),
                  diff_monthly[KEY].collapsed(['latitude','longitude'], iris.analysis.MEAN),
                  c='k', linewidth=2,
                  zorder=1
                  )

        plt.grid()

        plt.savefig(f'{outpath}/scatter_sm_diff_{KEY}.png')
        plt.close()
    

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
