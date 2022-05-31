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

band_width = 2


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

    # load preporcessed bands of lai obs
    obs_data = iris.load(f'/work/scratch-nopw/morobking/lai_bands_{band_width}degree.nc')
    for cube in obs_data:
        cube.attributes = None

    lai_obs_data = obs_data.concatenate_cube()

    print(f'{lai_obs_data=}')

    # assume single ensemble member for each MODEL!!!
    lai_model = {}
    sm_model = {}
    ts_model = {}
    tas_model = {}

    for KEY in loaded_data.keys(): # this is the model
        print(KEY)

        for ITEM in loaded_data[KEY].keys():
            print(ITEM) # this is the variable

            icc.add_year(loaded_data[KEY][ITEM], 'time')
            icc.add_month_number(loaded_data[KEY][ITEM], 'time')
            
            loaded_data[KEY][ITEM].coord('latitude').bounds = None
            loaded_data[KEY][ITEM].coord('longitude').bounds = None

            if 'ts' in ITEM:
                ts_model[KEY] = loaded_data[KEY][ITEM]

            if 'tas' in ITEM:
                tas_model[KEY] = loaded_data[KEY][ITEM]

            if 'lai' in ITEM:
                lai_model[KEY] = loaded_data[KEY][ITEM]
                
            if 'mrso' in ITEM:
                sm_model[KEY] = loaded_data[KEY][ITEM]          


    # for each model, average monthly LAI in 2 degree bands
    # note already done this and loading a seperate netcdf file for OBS
    lai_band_ave_model = {}
    sm_band_ave_model  = {}
    for MODEL in lai_model.keys():
        print(MODEL)
        this_banding = lat_band_ave(lai_model[MODEL])
        lai_band_ave_model[MODEL] = this_banding

        this_banding = lat_band_ave(sm_model[MODEL])
        sm_band_ave_model[MODEL] = this_banding

    print(f'{lai_band_ave_model=}')


    #### load CCI SM data and make zonal band
    sm_raw = iris.load_cube('/home/users/robking/cci_sm/cci_sm_*.nc')

    # sdepth1", taken out of CMIP6 CMOR table in Lmon in ESMValCore
    sm_band_obs = lat_band_ave(sm_raw)

    lai_obs_peak = find_peak_month(lai_obs_data)
    sm_obs_peak  = find_peak_month(sm_band_obs)

    print(f'{sm_obs_peak=}')

    lai_peak_month_model = {}
    sm_peak_month_model  = {}
    for MODEL in lai_band_ave_model.keys():
        this_peak = find_peak_month(lai_band_ave_model[MODEL])
        lai_peak_month_model[MODEL] = this_peak

        this_peak = find_peak_month(sm_band_ave_model[MODEL])

        sm_peak_month_model[MODEL] = this_peak

    print(f'{sm_peak_month_model=}')

    print(f'{lai_obs_peak=}')
    plot_mean_peak(lai_peak_month_model, lai_obs_peak, 'LAI',config)
    plot_mean_peak(sm_peak_month_model, sm_obs_peak, 'SM',config)

#    print(0/0)
def lat_band_ave(data):

    lats = np.arange(-60,80,band_width)
    output = iris.cube.CubeList()
    for item in lats:
        lat_con = iris.Constraint(latitude=lambda cell: item<=cell<item+band_width)
        
        this_cube = data.extract(lat_con)
        try:
            this_ave = this_cube.collapsed(['latitude','longitude'], iris.analysis.MEAN)
        except:
            continue
        output.append(this_ave)
        
    output = output.merge_cube()

    return output

def plot_mean_peak(model_peaks, obs_peaks, variable, config):

    plt.figure(figsize=(40,40))

    for i,MODEL in enumerate(model_peaks.keys()):
        
        plot_data = model_peaks[MODEL].collapsed('time', iris.analysis.MEAN)
        #plot_data.data = np.round(plot_data.data)
        iplt.plot(plot_data,
                  color = LINECOLOURS[i%10],
                  linestyle=LINE_STYLES[i//10],
                  linewidth=2,
                  label = MODEL[6:])

    plot_data = obs_peaks.collapsed('time', iris.analysis.MEAN)
    #plot_data.data = np.round(plot_data.data)
    iplt.plot(plot_data,
              color = 'black',
              linewidth=3,
              label = 'OBS')

    plt.xticks(range(1,13), MONTH_LIST, fontsize=20)
    plt.yticks(np.arange(-60,81,10), fontsize=20)
    plt.grid()

    plt.xlabel('Date', fontsize=24)
    plt.ylabel('Latitude', fontsize=24)

    outpath = config['plot_dir']
    plt.savefig(f'{outpath}/{variable}_peak_month_mean_clean.png',bbox_inches='tight')
    plt.legend(bbox_to_anchor=(1.01,0.0), loc='lower left', fontsize=20)            
    plt.savefig(f'{outpath}/{variable}_peak_month_mean_legend.png')
    plt.close()

def find_peak_month(data):

    try:
        icc.add_month_number(data, 'time')
    except:
        pass

    try:
        icc.add_year(data, 'time')
    except:
        pass
        
    this_cubelist = iris.cube.CubeList()
    for YEAR in np.unique(data.coord('year').points):
        cube = data.extract(iris.Constraint(year=YEAR))
        cube = cube.collapsed('time',ARGMAX)
            
        cube.coord('year').bounds=None
            
        cube.remove_coord('time')
        cube.remove_coord('month_number')
        cube = iris.util.new_axis(cube, scalar_coord='year')
        this_cubelist.append(cube)

    

    output = this_cubelist.concatenate_cube()
    output.data = np.ma.masked_outside(output.data, 1,12)


    years = output.coord('year').points
    dates = np.array([datetime.datetime(YEAR,1,1) for YEAR in years])
 
    data_dates = unit.date2num(dates,
                               'hours since 1970-01-01',
                               'gregorian'
                               )
    time_coord = iris.coords.AuxCoord(data_dates,
                                      var_name='time', standard_name='time', long_name='time',
                                      units=unit.Unit('hours since 1970-01-01', calendar='gregorian')
                                      )
    print(time_coord)
    output.add_aux_coord(time_coord,0)
    iris.util.promote_aux_coord_to_dim_coord(output,'time')

    return output

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
