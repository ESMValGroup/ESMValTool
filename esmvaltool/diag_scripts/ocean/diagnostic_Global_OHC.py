"""
Global Ocean Heat Content Time series diagnostics
=======================

Diagnostic to produce figures of the time development of globally-integrated
ocean heat content (OHC; units = Joules). These plots show time on the x-axis and
global OHC on the y-axis. In terms of preprocessors: the global OHC is integrated
over the chosen depth interval in _volume.py and it is then integrated over the
global ocean area in _area.py [using area_average with a SUM rather than MEAN].

Model output is grouped into specific sets (historical/RCP84/RCP25 etc) and plotted
in a manner suitable for IPCC reports. Gloal OHC is presented as an anomaly relative
to the mean across a selected historical period [anom_start anom_end] for each model.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, no depth component, and no
latitude or longitude coordinates.

An approproate preprocessor for the 3D temperature (thetao) field would be::

  preprocessors:
      prep_timeseries_1:
        custom_order: true
        extract_volume:
          z_min: 0.
          z_max: 700.
        depth_integration: # No coordz required in new versions
        average_region:
          coord1: longitude
          coord2: latitude
          operator: 'sum'

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Brodie Pearson (Brown - brodie_pearson@brown.edu)
based upon code by Lee de Mora (PML) ledm@pml.ac.uk
"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def timeplot(cube, **kwargs):
    """
        Make a time series plot
        
        Needed because iris version 1.13 fails due to the time axis.
        """
    if iris.__version__ > '2.0':
        qplt.plot(cube, kwargs)
    else:
        times = diagtools.cube_time_to_float(cube)
        plt.plot(times, cube.data, **kwargs)


def timeplotmultimodel(cube, **kwargs):
    """
    Create a time series plot for the multi model state.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.pyplot.plot function.
    These arguments are typically, color, linewidth, linestyle, etc...

    If there's only one datapoint in the cube, it is plotted as a
    horizontal line.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube

    """
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return
    
    hist_end_time=2004
    times = diagtools.cube_time_to_float(cube)
    #print(times)
    #print(cubedata)
    timesarray = np.asarray(times)
    #print(timesarray)
    #print(type(timesarray))
    #print(timesarray[(timesarray>=2004)])
    #print(cubedata[(timesarray<=2004)])
    
          
    #plt.plot(timesarray, cubedata, **kwargs)
    plt.plot(timesarray[(timesarray<=hist_end_time)], cubedata[(timesarray<=hist_end_time)], 'k')
    plt.plot(timesarray[(timesarray>=hist_end_time)], cubedata[(timesarray>=hist_end_time)], 'r')


def timeplotspread(cube, cube2, **kwargs):
    """
        Make a SHADED time series plot (BCP)
        This is used for plotting standard deviations about a mean
        
        Needed because iris version 1.13 fails due to the time axis.
        """
    if iris.__version__ > '2.0':
        qplt.plot(cube, kwargs)
    else:
        times = diagtools.cube_time_to_float(cube)
        timesarray = np.asarray(times)
        hist_end_time=2004
        cubedata = np.ma.array(cube.data)
        cube2data = np.ma.array(cube2.data)
        
        plt.fill_between(timesarray[(timesarray<=hist_end_time)], \
                         cube.data[(timesarray<=hist_end_time)],  \
                         cube2.data[(timesarray<=hist_end_time)], \
                         alpha=0.25, color='k')
        plt.fill_between(timesarray[(timesarray>=hist_end_time)], \
                         cubedata[(timesarray>=hist_end_time)], \
                         cube2data[(timesarray>=hist_end_time)], \
                         alpha=0.25, color='r')


def moving_average(cube, window):
    """
    Calculate a moving average.

    The window is a string which is a number and a measuremet of time.
    For instance, the following are acceptable window strings:

    * ``5 days``
    * ``12 years``
    * ``1 month``
    * ``5 yr``

    Also note the the value used is the total width of the window.
    For instance, if the window provided was '10 years', the the moving
    average returned would be the average of all values within 5 years
    of the central value.

    In the case of edge conditions, at the start an end of the data, they
    only include the average of the data available. Ie the first value
    in the moving average of a ``10 year`` window will only include the average
    of the five subsequent years.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube
    window: str
        A description of the window to use for the

    Returns
    ----------
    iris.cube.Cube:
        A cube with the movinage average set as the data points.

    """
    window = window.split()
    window_len = int(window[0]) / 2.
    win_units = str(window[1])

    if win_units not in [
            'days', 'day', 'dy', 'months', 'month', 'mn', 'years', 'yrs',
            'year', 'yr'
    ]:
        raise ValueError("Moving average window units not recognised: " +
                         "{}".format(win_units))

    times = cube.coord('time').units.num2date(cube.coord('time').points)

    datetime = diagtools.guess_calendar_datetime(cube)

    output = []

    times = np.array([
        datetime(time_itr.year, time_itr.month, time_itr.day, time_itr.hour,
                 time_itr.minute) for time_itr in times
    ])

    for time_itr in times:
        if win_units in ['years', 'yrs', 'year', 'yr']:
            tmin = datetime(time_itr.year - window_len, time_itr.month,
                            time_itr.day, time_itr.hour, time_itr.minute)
            tmax = datetime(time_itr.year + window_len, time_itr.month,
                            time_itr.day, time_itr.hour, time_itr.minute)

        if win_units in ['months', 'month', 'mn']:
            tmin = datetime(time_itr.year, time_itr.month - window_len,
                            time_itr.day, time_itr.hour, time_itr.minute)
            tmax = datetime(time_itr.year, time_itr.month + window_len,
                            time_itr.day, time_itr.hour, time_itr.minute)

        if win_units in ['days', 'day', 'dy']:
            tmin = datetime(time_itr.year, time_itr.month,
                            time_itr.day - window_len, time_itr.hour,
                            time_itr.minute)
            tmax = datetime(time_itr.year, time_itr.month,
                            time_itr.day + window_len, time_itr.hour,
                            time_itr.minute)

        arr = np.ma.masked_where((times < tmin) + (times > tmax), cube.data)
        output.append(arr.mean())
    cube.data = np.array(output)
    return cube


def make_time_series_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple time series plot for an indivudual model 1D cube.

    This tool loads the cube from the file, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)
        if 'moving_average' in cfg:
            cube_layer = moving_average(cube_layer, cfg['moving_average'])

        if multi_model:
            timeplot(cube_layer, label=metadata['dataset'], ls=':')
        else:
            timeplot(cube_layer, label=metadata['dataset'])

        # Add title, legend to plots
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer != '':
            if cube_layer.coords('depth'):
                z_units = cube_layer.coord('depth').units
            else:
                z_units = ''
            title = ' '.join([title, '(', layer, str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')
        plt.ylabel(str(cube_layer.units))

        # Determine image filename:
        if multi_model:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                prefix='MultiModel',
                suffix='_'.join(['timeseries',
                                 str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='timeseries_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multi_model_time_series(
        cfg,
        metadata,
):
    """
    Make a time series plot showing several preprocesssed datasets.

    This tool loads several cubes from the files, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.

    """

    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    layers = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])

        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[filename] = cubes
        for layer in cubes:
            layers[layer] = True

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)
    
    # Initialize switch for standard deviation shading
    cubep1=0
    cuben1=0
    
    serieshist = None
    seriesrcp26 = None
    seriesrcp45 = None
    seriesrcp60 = None
    seriesrcp85 = None
    anomaly_hist = 1
    anom_start = 1980
    anom_end = 2000
    seriesobs=None
    
    
    if anomaly_hist is not None:
        anomaly_dataset = None
        anomaly_value = None
        for layer in layers:
            for index, filename in enumerate(sorted(metadata)):
                if 'CanESM2' in metadata[filename]['dataset']: # CanESM2 is dummy observation for now
                    cube_temp=model_cubes[filename][layer]
                    times_anom_obs = np.asarray(diagtools.cube_time_to_float(cube_temp))
                    series_anom_obs = np.ma.array(cube_temp.data)
                    if anomaly_dataset is None:
                        anomaly_dataset = np.ma.array(metadata[filename]['dataset'])
                        anomaly_value = np.ma.array(np.mean(series_anom_obs[(times_anom_obs<anom_end) + (times_anom_obs>anom_start)], axis=0))
                    else:
                        anomaly_dataset = np.column_stack((anomaly_dataset, np.ma.array(metadata[filename]['dataset'])))
                        anomaly_value = np.column_stack((anomaly_value, np.ma.array(np.mean(series_anom_obs[(times_anom_obs<anom_end) + (times_anom_obs>anom_start)], axis=0))))
                else:
                    if 'historical' in metadata[filename]['exp']:
                        cube_temp=model_cubes[filename][layer]
                        times_anom = np.asarray(diagtools.cube_time_to_float(cube_temp))
                        series_anom = np.ma.array(cube_temp.data)
                        if anomaly_dataset is None:
                            anomaly_dataset = np.ma.array(metadata[filename]['dataset'])
                            anomaly_value = np.ma.array(np.mean(series_anom[(times_anom<anom_end) + (times_anom>anom_start)], axis=0))
                        else:
                            anomaly_dataset = np.column_stack((anomaly_dataset, np.ma.array(metadata[filename]['dataset'])))
                            anomaly_value = np.column_stack((anomaly_value, np.ma.array(np.mean(series_anom[(times_anom<anom_end) + (times_anom>anom_start)], axis=0))))

    
    # Make a plot for each layer (each layer is a group of simulations for a given RCP or for historical)
    for layer in layers:

        title = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('viridis')

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'

            # Take a moving average, if needed.
            if 'moving_average' in cfg:
                cube = moving_average(model_cubes[filename][layer],
                                      cfg['moving_average'])
            else:
                cube = model_cubes[filename][layer]
            
            # Select hist, rcp26, etc.
            #if 'rcp26' in metadata[filename]['exp']:
            #    cube_temp=model_cubes[filename][layer]
            #seriesrcp26 = np.append(seriesrcp26, np.ma.array(cube.data))
                    #print("seriesrcp26 = " seriesrcp26)

            if 'CanESM2' in metadata[filename]['dataset']: # CanESM2 is dummy OBS (project=OBS)
                cube_temp=model_cubes[filename][layer]
                times_obs = np.asarray(diagtools.cube_time_to_float(cube_temp))
                if seriesobs is None:
                    if anomaly_hist is None:
                        seriesobs = np.ma.array(cube_temp.data)
                    else:
                        seriesobs = np.ma.array(cube_temp.data) -anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]
                else:
                    if anomaly_hist is None:
                        seriesobs = np.column_stack((seriesobs, np.ma.array(cube_temp.data)))
                    else:
                        seriesobs = np.column_stack((seriesobs, np.ma.array(cube_temp.data)-anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]))
            else:
                if 'rcp45' in metadata[filename]['exp']:
                    cube_temp=model_cubes[filename][layer]
                    if seriesrcp45 is None:
                        if anomaly_hist is None:
                            seriesrcp45 = np.ma.array(cube_temp.data)
                        else:
                        #print(anomaly_value)
                        #print(anomaly_dataset)
                        #seriesrcp45 = np.ma.array(cube_temp.data)
                            seriesrcp45 = np.ma.array(cube_temp.data) -anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]
                    else:
                        if anomaly_hist is None:
                            seriesrcp45 = np.column_stack((seriesrcp45, np.ma.array(cube_temp.data)))
                        else:
                            seriesrcp45 = np.column_stack((seriesrcp45, np.ma.array(cube_temp.data)-anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]))
                elif 'rcp26' in metadata[filename]['exp']:
                    cube_temp=model_cubes[filename][layer]
                    times_rcp = np.asarray(diagtools.cube_time_to_float(cube_temp))
                    if seriesrcp26 is None:
                        if anomaly_hist is None:
                            seriesrcp26 = np.ma.array(cube_temp.data)
                        else:
                            seriesrcp26 = np.ma.array(cube_temp.data) -anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]
                    else:
                        if anomaly_hist is None:
                            seriesrcp26 = np.column_stack((seriesrcp26, np.ma.array(cube_temp.data)))
                        else:
                            seriesrcp26 = np.column_stack((seriesrcp26, np.ma.array(cube_temp.data)-anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]))
                elif 'rcp60' in metadata[filename]['exp']:
                    cube_temp=model_cubes[filename][layer]
                    if seriesrcp60 is None:
                        if anomaly_hist is None:
                            seriesrcp60 = np.ma.array(cube_temp.data)
                        else:
                            seriesrcp60 = np.ma.array(cube_temp.data) -anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]
                    else:
                        if anomaly_hist is None:
                            seriesrcp60 = np.column_stack((seriesrcp60, np.ma.array(cube_temp.data)))
                        else:
                            seriesrcp60 = np.column_stack((seriesrcp60, np.ma.array(cube_temp.data)-anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]))
                elif 'rcp85' in metadata[filename]['exp']:
                    cube_temp=model_cubes[filename][layer]
                    if seriesrcp85 is None:
                        if anomaly_hist is None:
                            seriesrcp85 = np.ma.array(cube_temp.data)
                        else:
                            seriesrcp85 = np.ma.array(cube_temp.data) -anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]
                    else:
                        if anomaly_hist is None:
                            seriesrcp85 = np.column_stack((seriesrcp85, np.ma.array(cube_temp.data)))
                        else:
                            seriesrcp85 = np.column_stack((seriesrcp85, np.ma.array(cube_temp.data)-anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]))
                elif 'historical' in metadata[filename]['exp']:
                    cube_temp=model_cubes[filename][layer]
                    times_hist = np.asarray(diagtools.cube_time_to_float(cube_temp))
                    if serieshist is None:
                        if anomaly_hist is None:
                            serieshist = np.ma.array(cube_temp.data)
                        else:
                            serieshist = np.ma.array(cube_temp.data) -anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]
                    else:
                        if anomaly_hist is None:
                            serieshist = np.column_stack((serieshist, np.ma.array(cube_temp.data)))
                        else:
                            serieshist = np.column_stack((serieshist, np.ma.array(cube_temp.data)-anomaly_value[(anomaly_dataset==metadata[filename]['dataset'])]))
                   


            title = metadata[filename]['long_name'] + ' Anomaly (1980-2000)'
            if layer != '':
                if model_cubes[filename][layer].coords('depth'):
                    z_units = model_cubes[filename][layer].coord('depth').units
                else:
                    z_units = ''

        # Constant to convert heat to Joules
        # approximate by assuming constant density for these purposes
        density_of_water = 1029
        specific_heat_capacity_water=3990
        OHC_factor=density_of_water*specific_heat_capacity_water

        serieshist=serieshist*OHC_factor
        seriesrcp26=seriesrcp26*OHC_factor
        seriesrcp45=seriesrcp45*OHC_factor
        seriesrcp60=seriesrcp60*OHC_factor
        seriesrcp85=seriesrcp85*OHC_factor
        seriesobs=seriesobs*OHC_factor

                        
        # Calculate the means and STDs of each scenario and plot
        if serieshist is not None:
            #plt.plot(times_hist, serieshist, c='k', ls='-', lw=0.25)
            logger.info('times_hist is', times_hist)
            logger.info('serieshist is', serieshist)
            logger.info('serieshistmean is', np.mean(serieshist, axis=1))
            plt.plot(times_hist, np.mean(serieshist, axis=1), c='k', ls='-', lw=1)
            plt.fill_between(times_hist, \
                         np.mean(serieshist, axis=1)+np.std(serieshist, axis=1), \
                         np.mean(serieshist, axis=1)-np.std(serieshist, axis=1), \
                         alpha=0.25, color='k')
            plot_details[1] = {
                                'c': 'k',
                                'ls': '-',
                                'lw': 1,
                                'label': 'Historical'
                                    }
        if seriesrcp26 is not None:
            #plt.plot(times_rcp, seriesrcp26, c='C0', ls='-', lw=0.25)
            plt.fill_between(times_rcp, \
                             np.mean(seriesrcp26, axis=1)+np.std(seriesrcp26, axis=1), \
                             np.mean(seriesrcp26, axis=1)-np.std(seriesrcp26, axis=1), \
                             alpha=0.25, color=(67./255, 147./255, 195./255))
            plt.plot(times_rcp, np.mean(seriesrcp26, axis=1), c=(0./255, 52./255, 102./255), ls='-', lw=1)
            plot_details[2] = {
                                'c': (0./255, 52./255, 102./255),
                                'ls': '-',
                                'lw': 1,
                                'label': 'RCP2.6'
                                }
        if seriesrcp45 is not None:
            #plt.plot(times_rcp, seriesrcp45, c='C9', ls='-', lw=0.25)
            plt.fill_between(times_rcp, \
                             np.mean(seriesrcp45, axis=1)+np.std(seriesrcp45, axis=1), \
                             np.mean(seriesrcp45, axis=1)-np.std(seriesrcp45, axis=1), \
                             alpha=0.25, color=(146./255, 197./255, 146./255))
            plt.plot(times_rcp, np.mean(seriesrcp45, axis=1), c=(112./255, 160./255, 205./255), ls='-', lw=1)
            plot_details[3] = {
                            'c': (112./255, 160./255, 205./255),
                            'ls': '-',
                            'lw': 1,
                            'label': 'RCP4.5'
                            }
        if seriesrcp60 is not None:
            #plt.plot(times_rcp, seriesrcp60, c='C9', ls='-', lw=0.25)
            plt.fill_between(times_rcp, \
                             np.mean(seriesrcp60, axis=1)+np.std(seriesrcp60, axis=1), \
                             np.mean(seriesrcp60, axis=1)-np.std(seriesrcp60, axis=1), \
                             alpha=0.25, color=(204./255, 174./255, 113./255))
            plt.plot(times_rcp, np.mean(seriesrcp60, axis=1), c=(196./255, 121./255, 0./255), ls='-', lw=1)
            plot_details[4] = {
                            'c': (196./255, 121./255, 0./255),
                            'ls': '-',
                            'lw': 1,
                            'label': 'RCP6.0'
                            }
        if seriesrcp85 is not None:
            #plt.plot(times_rcp, seriesrcp85, c='C1', ls='-', lw=0.25)
            plt.fill_between(times_rcp, \
                             np.mean(seriesrcp85, axis=1)+np.std(seriesrcp85, axis=1), \
                             np.mean(seriesrcp85, axis=1)-np.std(seriesrcp85, axis=1), \
                             alpha=0.25, color=(252./255, 209./255, 197./255))
                             #alpha=0.25, color='C1')
            plt.plot(times_rcp, np.mean(seriesrcp85, axis=1), c=(153./255, 0, 2./255), ls='-', lw=1)
            plot_details[5] = {
                            'c': (153./255, 0, 2./255),
                            'ls': '-',
                            'lw': 1,
                            'label': 'RCP8.5'
                            }
        if seriesobs is not None:
            plt.plot(times_obs, seriesobs, c='k', ls='-', lw=3)
            #plt.plot(times_rcp, np.mean(seriesrcp85, axis=1), c=(153./255, 0, 2./255), ls='-', lw=1)
            plot_details[6] = {
                            'c': 'k',
                            'ls': '-',
                            'lw': 3,
                            'label': 'Observations'
                            }


        
        # Add title, legend to plots
        if layer:
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])
        #plt.title(title)
        plt.legend(loc=1)
        plt.ylabel('Ocean Heat Content Anomaly (10$^{24}$ Joules)')
        plt.xlabel('Year')
        plt.ylim(-0.1e24, 1.1e24)
        #plt.ylim(-0.1*700, 0.1*700)
        plt.yticks((0,0.2e24,0.4e24,0.6e24,0.8e24,1e24),('0','0.2','0.4','0.6','0.8','1'))

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels_',
                suffix='_'.join(['timeseries',
                                 str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        # Resize and add legend outside thew axes.
        plt.gcf().set_size_inches(9., 6.)
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.15)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()


def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)
        
        metadatas = diagtools.get_input_files(cfg, index=index)
        
        #######
        # Multi model time series
        multi_model_time_series(
                                cfg,
                                metadatas,
                                )

        for filename in sorted(metadatas):
        
            logger.info('-----------------')
            logger.info(
                    'model filenames:\t%s',
                    filename,
                    )

            ######
            # Time series of individual model
            #make_time_series_plots(cfg, metadatas[filename], filename)
    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
