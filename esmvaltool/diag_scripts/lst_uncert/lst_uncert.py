"""
ESMValTool diagnostic for ESA CCI LST data.

The code uses the seperate Day and Night overpass monthly data.
The diagnostic calculates the average CCI LST togive an 'all time' value,
as well as propagating the gridbox uncertainity values to regional values.
The ouptput is a set of timeseries plots showing the CMIP5 and CMIP6 LST
and the CCI LST with its error plotted, as well as plot of if the model(s)
are warmer or cooler than CCI, or within the uncertainity bound.
"""

import logging

import iris
import matplotlib.pyplot as plt

import numpy as np

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

scalefactor = 0.001
# cant find this in the files, * this by all CCI values to get actual value


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
            data_type = 'CMIP6'
            # this way round means CMIP5 doesnt get counted twice

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
    """ESA CCI LST diagnostic with uncertainities included.

    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.

    Returns
    -------
    Figures made by make_plots.
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()

    loaded_data = {}

    data_ready = {'OBS': {},
                  'CMIP5': [],
                  'CMIP6': []
                  }
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors, data_type = _get_input_cubes(metadata)
        loaded_data[f'{data_type}_{dataset}'] = cubes

    # sort data into ensemble aves if needed
    for key in loaded_data.keys():
        if 'OBS' in key:
            continue
        if 'CMIP5' in key:
            model = 'CMIP5'
        if 'CMIP6' in key:
            model = 'CMIP6'

        cubes = iris.cube.CubeList()
        for i, item in enumerate(loaded_data[key].keys()):
            ensemble_coord = iris.coords.AuxCoord(i, var_name='realisation')
            loaded_data[key][item].add_aux_coord(ensemble_coord)
            cubes.append(loaded_data[key][item])
        data_ready[model] = cubes.merge_cube()

    # make model means if necessary
    data_means = {'CMIP5': [],
                  'CMIP6': [],
                  }

    for key in data_means.keys():
        if data_ready[key].ndim == 3:
            # no need to average
            data_means[key] = data_ready[key].collapsed(['latitude',
                                                         'longitude'],
                                                         iris.analysis.MEAN
                                                        )
        else:
            data_means[key] = data_ready[key].collapsed('realisation',
                                                        iris.analysis.MEAN
                                                        )
            data_means[key] = data_means[key].collapsed(['latitude',
                                                         'longitude'],
                                                         iris.analysis.MEAN
                                                        )

    # The Diagnostic uses CCI - MODEL
    # CMIP data had 360 day calendar, CCI data has 365 day calendar
    # Assume the loaded data is all the same shape

    # Calc CCI LST uncertainty
    uncerts = {'Day': iris.cube.CubeList(),
               'Night': iris.cube.CubeList()
               }

    # make the 'all time' LST average
    cci_lst = (loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'] +
               loaded_data['OBS_ESACCI_LST_UNCERTS']['tsNight']) / 2
    cci_lst_area_ave = cci_lst.collapsed(['latitude', 'longitude'],
                                         iris.analysis.MEAN
                                         )

    # make the gridbox total uncertainity
    # following conversation with Lizzie

    # make random uncert for gridbox from:
    # tsUnCorErr, tsLocalAtmErr and tsLocalSfcErr
    # then use tsLSSysErr with the random uncert to give
    # a gridbox total uncertanty.
    # sum in quadrature all region's total uncer to give
    # regional day/night uncertanty.
    # sum in quadrature the day and night regional total uncerts
    # to give the value to use.

    # this is gridbox total RANDOM uncert
    uncerts = {}
    for time in ['Day', 'Night']:
        num_times = len(loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'].coord('time').points)
        shape_2d = np.product(np.shape(loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'][0].data))
        n_array = np.array([shape_2d -
                            np.sum(loaded_data['OBS_ESACCI_LST_UNCERTS']['tsDay'][i].data.mask)
                            for i in range(num_times)]
                       )

        this_cube = (scalefactor * \
                     loaded_data['OBS_ESACCI_LST_UNCERTS'][f'tsUnCorErr{time}'])**2 + \
                    (scalefactor * \
                     loaded_data['OBS_ESACCI_LST_UNCERTS'][f'tsLocalAtmErr{time}'])**2 + \
                    (scalefactor * \
                     loaded_data['OBS_ESACCI_LST_UNCERTS'][f'tsLocalSfcErr{time}'])**2

        iris.analysis.maths.exponentiate(this_cube, 0.5, in_place=True)
        # this_cube is now the random uncert for each grid box

        # now get area random uncert
        this_cube = (this_cube**2).collapsed(['latitude', 'longitude'],
                                             iris.analysis.SUM
                                             )
        for i in range(num_times):
            this_cube[i].data = (1 / n_array[i]**2) * this_cube[i].data
        iris.analysis.maths.exponentiate(this_cube, 0.5, in_place=True)

        # this_cube is now the area random uncert

        # sum in quadrature to get area total uncert using 
        # the 'constant' tsLSSysErr
        this_cube = this_cube**2 + \
                  loaded_data['OBS_ESACCI_LST_UNCERTS'][f'tsLSSysErr{time}']**2
        this_cube = iris.analysis.maths.exponentiate(this_cube,
                                                     0.5,
                                                     in_place=True
                                                     )

        uncerts[time] = this_cube

    # now sum in quadrature day and night values to give total all time uncert
    total_uncert = uncerts['Day']**2 + uncerts['Night']**2
    iris.analysis.maths.exponentiate(total_uncert, 0.5, in_place=True)

    make_plots(cci_lst_area_ave, data_means, total_uncert, config)


def make_plots(cci_lst, data_means, total_uncert, config):
    """Create and save the output figure.

    PLOT 1
    The plot both CMIP5 and CMIP6 model LST and the CCI LST with
    +/- one total  error.

    PLOT 2
    This is schematic of whether the CMIP5 and CMIP6 model LST is warmer,
    cooler, or within the CCI error.

    PLOT 3
    Combines plot 1 and plot2 in a single figure.

    Inputs:
    config = The config dictionary from the preprocessor
    Outputs:
    Saved figure
    """
    num_times = len(cci_lst.coord('time').points)

    tab_cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    colours = {'OBS': tab_cols[0],
               'CMIP5': tab_cols[1],
               'CMIP6': tab_cols[2],
               'high': 'red',
               'inside': 'grey',
               'low': 'blue',
               }

    cci_high = (cci_lst + total_uncert).data
    cci_low = (cci_lst - total_uncert).data

    # calc overlaps
    overlaps5 = []
    overlaps6 = []
    for i in range(num_times):
        if data_means['CMIP5'][i].data < cci_low[i]:
            overlaps5.append(-1)
        elif data_means['CMIP5'][i].data > cci_high[i]:
            overlaps5.append(1)
        else:
            overlaps5.append(0)

    for i in range(num_times):
        if data_means['CMIP6'][i].data < cci_low[i]:
            overlaps6.append(-1)
        elif data_means['CMIP6'][i].data > cci_high[i]:
            overlaps6.append(1)
        else:
            overlaps6.append(0)

    # fig 1 = just timeseries
    # fig 2 = just overlaps
    # fig 3 = both on a shared x axis
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True,
                             figsize=(20, 15))
    fig2, ax2 = plt.subplots(nrows=1, ncols=1, sharex=True,
                             figsize=(20, 15))
    fig3, ax3 = plt.subplots(nrows=2, ncols=1, sharex=True,
                             figsize=(20, 15))

    ax1.plot(range(num_times), cci_lst.data, label='OBS',
             color=colours['OBS'],
             linewidth=2
             )
    ax1.fill_between(range(num_times), cci_low, cci_high,
                     color=colours['OBS'],
                     alpha=0.5
                     )

    ax1.plot(range(num_times), data_means['CMIP5'].data,
             label='CMIP5',
             color=colours['CMIP5']
             )
    ax1.plot(range(num_times), data_means['CMIP6'].data,
             label='CMIP6',
             color=colours['CMIP6']
             )
    ax1.legend(fontsize=18)

    # for figure 3
    ax3[0].plot(range(num_times), cci_lst.data, label='OBS',
                color=colours['OBS'],
                linewidth=2
                )
    ax3[0].fill_between(range(num_times), cci_low, cci_high,
                        color=colours['OBS'], alpha=0.5
                        )

    ax3[0].plot(range(num_times), data_means['CMIP5'].data,
                label='CMIP5',
                color=colours['CMIP5']
                )
    ax3[0].plot(range(num_times), data_means['CMIP6'].data,
                label='CMIP6',
                color=colours['CMIP6']
                )
    ax3[0].legend(fontsize=18)

    min_point = np.min([cci_low, cci_high,
                        data_means['CMIP5'].data, data_means['CMIP6'].data])
    max_point = np.max([cci_low, cci_high,
                        data_means['CMIP5'].data, data_means['CMIP6'].data])

    ax1.set_yticks(np.arange(5 * ((min_point - 5) // 5),
                             5 * ((max_point + 5) // 5) + 7,
                             5))
    ax1.set_yticklabels(np.arange(5 * ((min_point - 5) // 5),
                                  5 * ((max_point + 5 + 1) // 5) + 7,
                                  5),
                        fontsize=18)
    ax1.set_ylim((min_point - 6, max_point + 7))

    ax3[0].set_yticks(np.arange(5 * ((min_point - 5) // 5),
                                5 * ((max_point + 5) // 5) + 7,
                                5))
    ax3[0].set_yticklabels(np.arange(5 * ((min_point - 5) // 5),
                                     5 * ((max_point + 5 + 1) // 5) + 7,
                                     5),
                           fontsize=18)
    ax3[0].set_ylim((min_point - 6, max_point + 7))

    facecolor5 = []
    for i in range(num_times):
        if overlaps5[i] == -1:
            facecolor5.append(colours['low'])
        elif overlaps5[i] == 1:
            facecolor5.append(colours['high'])
        else:
            facecolor5.append(colours['inside'])

    facecolor6 = []
    for i in range(num_times):
        if overlaps6[i] == -1:
            facecolor6.append(colours['low'])
        elif overlaps6[i] == 1:
            facecolor6.append(colours['high'])
        else:
            facecolor6.append(colours['inside'])

    marker_size = 150
    marker_shape = 'o'

    ax2.scatter(range(num_times), [0 for i in range(num_times)],
                c=facecolor5,
                s=marker_size,
                edgecolor='black',
                marker=marker_shape
                )
    ax2.scatter(range(num_times), [1 for i in range(num_times)],
                c=facecolor6,
                s=marker_size,
                edgecolor='black',
                marker=marker_shape
                )

    ax3[1].scatter(range(num_times), [0 for i in range(num_times)],
                   c=facecolor5,
                   s=marker_size,
                   edgecolor='black',
                   marker=marker_shape
                   )
    ax3[1].scatter(range(num_times), [1 for i in range(num_times)],
                   c=facecolor6,
                   s=marker_size,
                   edgecolor='black',
                   marker=marker_shape
                   )

    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(['CMIP5', 'CMIP6'], fontsize=18)
    ax3[1].set_yticks([0, 1])
    ax3[1].set_yticklabels(['CMIP5', 'CMIP6'], fontsize=18)

    # make X ticks
    x_tick_list = []
    time_list = cci_lst.coord('time').units.num2date(
        cci_lst.coord('time').points)
    for item in time_list:
        if item.month == 1:
            x_tick_list.append(item.strftime('%Y %b'))
        elif item.month == 7:
            x_tick_list.append(item.strftime('%b'))
        else:
            x_tick_list.append('')

    # common stuff for all three plots
    ax1.set_xlim((-1, num_times + 1))
    ax2.set_xlim((-1, num_times + 1))
    ax2.set_ylim((-0.5, 1.5))
    ax3[1].set_xlim((-1, num_times + 1))
    ax3[1].set_ylim((-0.5, 1.5))

    ax1.set_xticks(range(num_times))
    ax1.set_xticklabels(x_tick_list, fontsize=18)
    ax2.set_xticks(range(num_times))
    ax2.set_xticklabels(x_tick_list, fontsize=18)
    ax3[0].set_xticks(range(num_times))
    ax3[1].set_xticks(range(num_times))
    ax3[0].set_xticklabels(['' for x in x_tick_list])
    ax3[1].set_xticklabels(x_tick_list, fontsize=18)

    ax1.set_xlabel('Date', fontsize=22)
    ax2.set_xlabel('Date', fontsize=22)
    ax3[1].set_xlabel('Date', fontsize=22)

    ax2.set_ylabel('Model', fontsize=22)
    ax3[1].set_ylabel('Model', fontsize=22)

    ax1.set_ylabel('LST (K)', fontsize=22)
    ax3[0].set_ylabel('LST (K)', fontsize=22)

    ax1.grid()
    ax2.grid()
    ax3[0].grid()
    ax3[1].grid()

    ax2.text(17, 0.60, 'Model cooler than Obs',
             fontsize=20,
             color=colours['low']
             )
    ax2.text(17, 0.50,
             'Model inside Obs uncertainty',
             fontsize=20,
             color=colours['inside']
             )
    ax2.text(17, 0.40,
             'Model warmer than Obs',
             fontsize=20,
             color=colours['high']
             )

    ax3[1].text(17, 0.60, 'Model cooler than Obs',
                fontsize=20,
                color=colours['low']
                )
    ax3[1].text(17, 0.50, 'Model inside Obs uncertainty',
                fontsize=20,
                color=colours['inside']
                )
    ax3[1].text(17, 0.40, 'Model warmer than Obs',
                fontsize=20,
                color=colours['high']
                )

    fig1.suptitle('ESACCI LST, CMIP5 and CMIP6 LST', fontsize=24)
    fig2.suptitle('ESACCI LST, CMIP5 and CMIP6 LST', fontsize=24)
    fig3.suptitle('ESACCI LST, CMIP5 and CMIP6 LST', fontsize=24)

    outpath = config['plot_dir']
    plt.figure(1)
    plt.savefig(f'{outpath}/timeseries_plot1.png')
    plt.close()
    plt.figure(2)
    plt.savefig(f'{outpath}/timeseries_plot2.png')
    plt.close()
    plt.figure(3)
    plt.savefig(f'{outpath}/timeseries_plot3.png')
    plt.close()


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
