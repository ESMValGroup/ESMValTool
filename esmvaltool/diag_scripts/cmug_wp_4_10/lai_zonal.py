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
    #print(metadata)
    for attributes in metadata:
        #print(attributes)
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
    np.savetxt('/home/users/robking/csv_data/LAI_value_months_OBS.csv',
                   lai_obs_data.data, delimiter=',', newline='\n', header='',
                   fmt='%5.3f')
    #
    # make mean lst from day and night
  
    print(loaded_data)
    # # gah, no idea why the lst data has a messed up longitude coordinate, but this is a fix
    # loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'].coord('longitude').points = np.arange(0,360,1/20.)
    # loaded_data['OBS_ESACCI_LST_UNCERTS']['tsNight'].coord('longitude').points = np.arange(0,360,1/20.)
    # loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'].coord('longitude').circular=True
    # loaded_data['OBS_ESACCI_LST_UNCERTS']['tsNight'].coord('longitude').circular=True

    # iris.util.promote_aux_coord_to_dim_coord(loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'], 'longitude')
    # iris.util.promote_aux_coord_to_dim_coord(loaded_data['OBS_ESACCI_LST_UNCERTS']['tsNight'], 'longitude')

    # # this makes the all time LST from Day and Night time overpasses
    # cci_lst = loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'] + loaded_data['OBS_ESACCI_LST_UNCERTS']['tsNight']
    # cci_lst = cci_lst/2
    # cci_lst.standard_name = 'surface_temperature'

    # assume single ensemble member for each MODEL!!!
    lai_model = {}
    sm_model = {}
    ts_model = {}
    tas_model = {}

    for KEY in loaded_data.keys(): # this is the model
        if 'OBS' in KEY: continue
        for ITEM in loaded_data[KEY].keys():

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
        #print(MODEL)
        this_banding = lat_band_ave(lai_model[MODEL])
        lai_band_ave_model[MODEL] = this_banding

        this_banding = lat_band_ave(sm_model[MODEL])
        sm_band_ave_model[MODEL] = this_banding


    #### load CCI SM data and make zonal band
    sm_raw = iris.load_cube('/work/scratch-nopw/morobking/cci_sm/cci_sm_*.nc')

    # sdepth1", taken out of CMIP6 CMOR table in Lmon in ESMValCore
    sm_band_obs = lat_band_ave(sm_raw)
    np.savetxt('/home/users/robking/csv_data/sm_value_months_OBS.csv',
                   sm_band_obs.data, delimiter=',', newline='\n', header='',
                   fmt='%5.3f')


    lai_obs_peak, lai_obs_min, lai_obs_max = find_peak_month(lai_obs_data)
    sm_obs_peak, sm_obs_min, sm_obs_max  = find_peak_month(sm_band_obs)

    #print(f'{sm_obs_peak=}')

    lai_peak_month_model = {}
    sm_peak_month_model  = {}

    lai_min_model = {}
    lai_max_model = {}
    sm_min_model  = {}
    sm_max_model  = {}

    for MODEL in lai_band_ave_model.keys():
        this_peak, this_min, this_max = find_peak_month(lai_band_ave_model[MODEL])
        lai_peak_month_model[MODEL] = this_peak
        lai_min_model[MODEL] = this_min
        lai_max_model[MODEL] = this_max

        this_peak, this_min, this_max = find_peak_month(sm_band_ave_model[MODEL])

        sm_peak_month_model[MODEL] = this_peak
        sm_min_model[MODEL] = this_min
        sm_max_model[MODEL] = this_max


        np.savetxt(f'/home/users/robking/csv_data/lai_peak_value_months_{MODEL}.csv',
                   lai_max_model[MODEL].data, delimiter=',', newline='\n', header='',
                   fmt='%5.3f')

        np.savetxt(f'/home/users/robking/csv_data/sm_MIN_value_months_{MODEL}.csv',
                   sm_min_model[MODEL].data, delimiter=',', newline='\n', header='',
                   fmt='%5.3f')

    #iris.save(lai_obs_peak,'/home/users/robking/lai_obs_peak.nc')
    np.savetxt('/home/users/robking/csv_data/lai_obs_peak_months.csv',
                  lai_obs_peak.data, delimiter=',', newline='\n', header='',
               fmt='%2i')
    for KEY in lai_peak_month_model.keys():
        np.savetxt('/home/users/robking/csv_data/lai_model_peak_months_{KEY}.csv',
                   lai_peak_month_model[KEY].data, delimiter=',', newline='\n', header='',
                   fmt='%2i')

    # find average months where needed
    
    lai_peak_months_ave_obs = month_ave_cube_2d_to_1d(lai_obs_peak)
    lai_peak_months_ave_model = {}
    for KEY in lai_peak_month_model.keys():
        lai_peak_months_ave_model[KEY] = month_ave_cube_2d_to_1d(lai_peak_month_model[KEY])

        
    np.savetxt('/home/users/robking/csv_data/lai_obs_peak_month_average.csv',
                  lai_peak_months_ave_obs.data, delimiter=',', newline='\n', header='',
               fmt='%2i')
    for KEY in lai_peak_month_model.keys():
        np.savetxt(f'/home/users/robking/csv_data/lai_model_peak_month_average_{KEY}.csv',
                   lai_peak_months_ave_model[KEY].data, delimiter=',', newline='\n', header='',
                   fmt='%2i')


 
    plot_mean_peak(lai_peak_months_ave_model, lai_peak_months_ave_obs, 'LAI_PEAK',config, False)
#    plot_mean_peak(lai_peak_month_model, lai_peak_months_obs, 'LAI_PEAK',config)
#    plot_mean_peak(sm_peak_month_model, sm_obs_peak, 'SM_PEAK',config)


#    plot_mean_peak(lai_max_model, lai_obs_max, 'LAI_MAX',config)
#    plot_mean_peak(sm_min_model, sm_obs_min, 'SM_MIN',config)
#    plot_mean_peak(sm_max_model, sm_obs_max, 'SM_MAX',config)



    
#    print(0/0)



def lat_band_ave(data):

    lats = np.arange(-60,70,band_width)
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

def plot_mean_peak(model_peaks, obs_peaks, variable, config, mean=True):

    plt.figure(figsize=(40,40))

    for i,MODEL in enumerate(model_peaks.keys()):
        
        if mean:
            plot_data = model_peaks[MODEL].collapsed('time', iris.analysis.MEAN)
        else:
            plot_data = model_peaks[MODEL]
        #plot_data.data = np.round(plot_data.data)
        iplt.plot(plot_data,
                  color = LINECOLOURS[i%10],
                  linestyle=LINE_STYLES[i//10],
                  linewidth=2,
                  label = MODEL[6:])

    if mean:
        plot_data = obs_peaks.collapsed('time', iris.analysis.MEAN)
    else:
        plot_data = obs_peaks
    #plot_data.data = np.round(plot_data.data)
    iplt.plot(plot_data,
              color = 'black',
              linewidth=3,
              label = 'OBS')


    if 'PEAK' in variable:
        plt.xticks(range(1,13), MONTH_LIST, fontsize=20)
        plt.xlabel('Date', fontsize=24)

    if variable == 'LAI_MAX':
        plt.xticks(range(0,7), fontsize=20)
        plt.xlabel('LAI', fontsize=20)


    if variable == 'SM_MAX':
        plt.xlabel('Soil Moisture', fontsize=20)

    plt.yticks(np.arange(-60,81,10), fontsize=20)
    plt.grid()

    
    
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
    this_min      = iris.cube.CubeList()
    this_max      = iris.cube.CubeList()
    for YEAR in np.unique(data.coord('year').points):
        this_data = data.extract(iris.Constraint(year=YEAR))
        
        cube = this_data.collapsed('time',ARGMAX)
        
        cube_min = this_data.collapsed('time', iris.analysis.MIN)
        cube_max = this_data.collapsed('time', iris.analysis.MAX)

        cube.coord('year').bounds=None
        cube_min.coord('year').bounds=None
        cube_max.coord('year').bounds=None
            
        cube.remove_coord('time')
        cube.remove_coord('month_number')
        cube = iris.util.new_axis(cube, scalar_coord='year')
        this_cubelist.append(cube)

        cube_min.remove_coord('time')
        cube_min.remove_coord('month_number')
        cube_min = iris.util.new_axis(cube_min, scalar_coord='year')
        this_min.append(cube_min)

        cube_max.remove_coord('time')
        cube_max.remove_coord('month_number')
        cube_max = iris.util.new_axis(cube_max, scalar_coord='year')
        this_max.append(cube_max)    

    output = this_cubelist.concatenate_cube()
    

    output_min = this_min.concatenate_cube()
    output_max = this_max.concatenate_cube()


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
    output.add_aux_coord(time_coord,0)
    iris.util.promote_aux_coord_to_dim_coord(output,'time')
    output_min.add_aux_coord(time_coord,0)
    iris.util.promote_aux_coord_to_dim_coord(output_min,'time')
    output_max.add_aux_coord(time_coord,0)
    iris.util.promote_aux_coord_to_dim_coord(output_max,'time')

    return output, output_min, output_max


def month_ave_cube_2d_to_1d(cube):
    # collapse a cube on time (years) with month_ave_calc value
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxx')


    print(cube)
    lat_coord = cube.coord('latitude')
    lat_coord.bounds = None
    print(lat_coord.points)
    
    output_cube = iris.cube.Cube(np.zeros(len(lat_coord.points)),
                                 #dim_coords_and_dims = (lat_coord,0),
                                 var_name = 'average_month'
                                 )
    print(output_cube)

    for i in range(len(lat_coord.points)):
        months = cube[:,i].data
        ave_months = month_ave_calc(months)
        
        output_cube.data[i] = ave_months

    print(output_cube)
    output_cube.add_aux_coord(lat_coord,0)
    iris.util.promote_aux_coord_to_dim_coord(output_cube,'latitude')
    print(output_cube.data)

    return output_cube
    
    

def month_ave_calc(month_list):
    # input is a list of integer months 0-11 corresponding to Jan-Dec
    # dont have to be integers for this to work but integers come from 
    # the underlying data

    # idea from Laura at AVD surgey 1st June 2022, Thank You!
    # consider each month as unit vector, then aveage individual
    # x & y components
    # this 'average' vector correcpond to the angle of the average month
    # this will handle inputs that cross Dec-Jan ie 11-0

    # 1) convert each month to a unit vecor
    X = [] 
    Y = []
    for item in month_list:
        x_comp = np.cos(item*2*np.pi/12)
        y_comp = np.sin(item*2*np.pi/12)

        X.append(x_comp)
        Y.append(y_comp)

    # 2) calculate the mean vector
    X_mean = np.mean(np.array(X))
    Y_mean = np.mean(np.array(Y))

    # angle is the arctan of Y_mean and X_mean
    # use acrtan2 to get sector correct,
    # convert radians to 0-11
    # modulo 12 and %12 mean output is in 0 <= month < 12
    # ie doesnt inc 12 itself    

    return np.mod(np.arctan2(Y_mean,X_mean)*12/np.pi/2,12)%12

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
