#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to plot figure 10 of IPCC AR5 SPM.

Description
-----------
Calculate and plot the global mean surface air temperature increase as a
function of cumulative total global CO2 emissions (see IPCC AR5 WG1 fig.
SPM10).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
save : dict, optional
    Keyword arguments for the `fig.saveplot()` function.
axes_functions : dict, optional
    Keyword arguments for the plot appearance functions.
dataset_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
matplotlib_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python.matplotlib`).

"""

import logging
import os
from pprint import pformat

import iris
import iris.coord_categorisation as iris_cat
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger, extract_variables, get_diagnostic_filename,
    get_plot_filename, group_metadata, io, plot, run_diagnostic,
    variables_available)

logger = logging.getLogger(os.path.basename(__file__))

# # ESMValTool python packages
# from auxiliary import info, warning
# from constants import Constant
# from esmval_lib import ESMValProject

# # Iris
# import iris
# import iris.coord_categorisation as iris_cat

# # Basic python packages
# from collections import OrderedDict
# import numpy as np
# import os
# import scipy.spatial

# # Matplotlib with 'Agg' backend
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.lines as mlines
# import matplotlib.patches as mpatches


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'short_name')

    # Check if necessary variables are available
    necessary_vars = ['co2mass', 'tas', 'nbp_grid', 'fgco2_grid']
    if not variables_available(cfg, necessary_vars):
        raise ValueError(
            f"This diagnostic needs the variables {necessary_vars}")

    # Iterate over all variables
    for (var, datasets) in grouped_data.items():
        for dataset in datasets:
            cube = iris.load_cube(dataset['filename'])

            # co2mass: no preprocessing necessary
            if var == 'co2mass':
                pass

            # tas: weighted mean
            elif var == 'tas':
                if dataset.get('fx_files', {}).get('areacella'):
                    area_cube = iris.load_cube(
                        dataset['fx_files']['areacella'])
                    area_weights = np.broadcast_to(area_cube.data, cube.shape)
                else:
                    area_weights = iris.analysis.cartography.area_weights(cube)
                    logger.warning(
                        "'areacella' file not found, estimated cell areas by "
                        "bounds")
                cube = cube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.MEAN,
                                      weights=area_weights)

            # nbp_grid and fgco2_grid: simple sum
            elif var in ('nbp_grid', 'fgco2_grid'):
                cube = cube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.SUM)

            # Other variables not supported
            else:
                logger.warning("Variable '%s' not supported, skipping", var)
                continue

            # Add aux coordinate for for decade

            def time_to_decade(coord, value):
                """Convert time into decade, e.g. 2050 means 2040-2049."""
                dec = 10
                year = coord.units.num2date(value).year
                return year + dec - year % dec

            iris_cat.add_categorised_coord(
                cube, 'decade', 'time', time_to_decade, units='yr')

            # Decadal aggregation
            if var == 'co2mass':
                initial_co2mass = cube.data[0]

                def end_of_decade_value(data, axis):
                    """Return last element of decade."""
                    return numpy.take(data, -1, axis=axis)

                aggregator = iris.analysis.Aggregator('decade: point',
                                                      end_of_decade_value)
                cube = cube.aggregated_by('decade', aggregator)
                cube.data -= np.hstack(initial_co2mass, cube.data[:-1])
            else:
                cube = cube.aggregated_by('decade', iris.analysis.MEAN)
            dataset['cube'] = cube

            print(cube)
            print(cube.data)
            print(cube.coord('time'))
            print(cube.coord('decade'))
    logger.info(pformat(grouped_data))

    # Get maximum total number of decades
    # try:
    #     max_decades[exp] = max(max_decades[exp],
    #                            len(cube.coord(DECADE).points))
    # except KeyError:
    #     max_decades[exp] = len(cube.coord(DECADE).points)

    # ###########################################################################
    # # Process data
    # ###########################################################################
    # CO2EMISS = 'co2emiss'
    # CO2EMISS_CUM = 'co2emiss_cum'
    # MM_MEAN = 'Multimodel mean'

    # # Append RCP85 scenarios to end of historical period (2005-2009)
    # # tas: use weighted mean
    # for m in data[VARS['TAS']][EXPS['HIST']]:
    #     if (m not in data[VARS['TAS']][EXPS['RCP85']]):
    #         continue
    #     data[VARS['TAS']][EXPS['HIST']][m][-1] = (
    #         hist_last_dec[m] * data[VARS['TAS']][EXPS['HIST']][m][-1] +
    #         rcp85_first_dec[m] * data[VARS['TAS']][EXPS['RCP85']][m][0]) / \
    #         (hist_last_dec[m] + rcp85_first_dec[m])
    #     data[VARS['TAS']][EXPS['RCP85']][m] = np.delete(
    #         data[VARS['TAS']][EXPS['RCP85']][m], 0)
    #     decades[VARS['TAS']][EXPS['RCP85']][m] = np.delete(
    #         decades[VARS['TAS']][EXPS['RCP85']][m], 0)
    # max_decades[EXPS['RCP85']] -= 1

    # # Calculate temperature anomaly relative to base period
    # # (including multi model mean and decades)
    # for exp in data[VARS['TAS']]:
    #     decades_mean = np.zeros(max_decades[exp])
    #     tas_mean = np.zeros(max_decades[exp])
    #     number_of_models = np.zeros(max_decades[exp])
    #     for model in data[VARS['TAS']][exp]:
    #         len_data = len(data[VARS['TAS']][exp][model])
    #         tas_data = np.pad(data[VARS['TAS']][exp][model],
    #                           (0, max_decades[exp]-len_data),
    #                           'constant')
    #         decade_data = np.pad(decades[VARS['TAS']][exp][model],
    #                              (0, max_decades[exp]-len_data),
    #                              'constant')
    #         number_of_models += np.pad(np.ones(len_data),
    #                                    (0, max_decades[exp]-len_data),
    #                                    'constant')

    #         # Calculate mean
    #         decades_mean += decade_data
    #         tas_mean += tas_data

    #     decades_mean /= number_of_models
    #     tas_mean /= number_of_models

    #     # Save tas relative to base period
    #     base_tas = tas_mean[0]
    #     for model in data[VARS['TAS']][exp]:
    #         E.add_to_data(data, VARS['TAS'], exp, model,
    #                       data[VARS['TAS']][exp][model]-base_tas, warn=False)
    #     E.add_to_data(decades, VARS['TAS'], exp, MM_MEAN, decades_mean)
    #     E.add_to_data(data, VARS['TAS'], exp, MM_MEAN, tas_mean-base_tas)

    # # Convert carbon fluxes to annual changes
    # seconds = DEC * 365.0 * 24.0 * 3600.0  # Decade in seconds

    # # nbp (Note: units already in kgC, not kgCO2)
    # for exp in data[VARS['NBP']]:
    #     for model in data[VARS['NBP']][exp]:
    #         data[VARS['NBP']][exp][model] *= seconds * \
    #             areas[VARS['NBP']][exp][model] / 1E12
    # units[VARS['NBP']] = 'GtC'

    # # fgco2 (Note: units already in kgC, not kgCO2)
    # for exp in data[VARS['FGCO2']]:
    #     for model in data[VARS['FGCO2']][exp]:
    #         data[VARS['FGCO2']][exp][model] *= seconds * \
    #             areas[VARS['FGCO2']][exp][model] / 1E12
    # units[VARS['FGCO2']] = 'GtC'

    # # Convert CO2 mass changes from [kgCO2] to [GtC]
    # for exp in data[VARS['CO2MASS']]:
    #     data[VARS['CO2MASS']][exp][VARS['CO2MASS']] *= Constant.kgCO2_to_GtC
    # units[VARS['CO2MASS']] = 'GtC'

    # # Calculate annual CO2 emissions
    # # dC_E = dC_A + dC_L + dC_L = co2mass + nbp + fgco2
    # for exp in data[VARS['CO2MASS']]:
    #     for model in data[VARS['NBP']][exp]:
    #         co2emiss_data = data[VARS['CO2MASS']][exp][VARS['CO2MASS']] + \
    #                         data[VARS['NBP']][exp][model] + \
    #                         data[VARS['FGCO2']][exp][model]
    #         E.add_to_data(data, CO2EMISS, exp, model, co2emiss_data)
    # units[CO2EMISS] = 'GtC'

    # # Append RCP85 scenarios to end of historical period (2005-2009)
    # # co2emiss: simply add years
    # for m in data[CO2EMISS][EXPS['HIST']]:
    #     if (m not in data[CO2EMISS][EXPS['RCP85']]):
    #         continue
    #     data[CO2EMISS][EXPS['HIST']][m][-1] = (
    #         data[CO2EMISS][EXPS['HIST']][m][-1] +
    #         data[CO2EMISS][EXPS['RCP85']][m][0])
    #     data[CO2EMISS][EXPS['RCP85']][m] = np.delete(
    #         data[CO2EMISS][EXPS['RCP85']][m], 0)

    # # Calculate cumulative CO2 emissions (including multi model mean)
    # for exp in data[CO2EMISS]:
    #     co2emiss_cum_mean = np.zeros(max_decades[exp])
    #     for model in data[CO2EMISS][exp]:
    #         co2emiss_cum_data = np.cumsum(data[CO2EMISS][exp][model])
    #         len_data = len(co2emiss_cum_data)
    #         if (len_data < max_decades[exp]):
    #             warning("Model {0} [{1}] does not ".format(model, exp) +
    #                     "cover all decadces", verbosity, 0, exit_on_warning)
    #             co2emiss_cum_data = np.pad(co2emiss_cum_data,
    #                                        (0, max_decades[exp]-len_data),
    #                                        'constant')
    #         co2emiss_cum_mean += co2emiss_cum_data
    #         E.add_to_data(data, CO2EMISS_CUM, exp, model, co2emiss_cum_data)
    #     co2emiss_cum_mean /= len(data[CO2EMISS_CUM][exp])

    #     # Save co2emiss_cum relative to base period
    #     base_co2emiss_cum = co2emiss_cum_mean[0]
    #     for model in data[CO2EMISS_CUM][exp]:
    #         E.add_to_data(data, CO2EMISS_CUM, exp, model,
    #                       data[CO2EMISS_CUM][exp][model]-base_co2emiss_cum,
    #                       warn=False)
    #     E.add_to_data(decades, CO2EMISS_CUM, exp, MM_MEAN,
    #                   decades[VARS['TAS']][exp][MM_MEAN])
    #     E.add_to_data(data, CO2EMISS_CUM, exp, MM_MEAN,
    #                   co2emiss_cum_mean-base_co2emiss_cum)
    # units[CO2EMISS_CUM] = 'GtC'

    # # Start RCP scenarios at end of historical period
    # for exp in data[VARS['TAS']]:
    #     if (not exp.startswith('rcp')):
    #         continue
    #     for model in data[VARS['TAS']][exp]:
    #         if (model not in data[VARS['TAS']][EXPS['HIST']]):
    #             continue
    #         tas_hist = data[VARS['TAS']][EXPS['HIST']][model][-1]
    #         tas_rcp = data[VARS['TAS']][exp][model]
    #         tas_rcp += tas_hist - tas_rcp[0]
    # for exp in data[CO2EMISS_CUM]:
    #     if (not exp.startswith('rcp')):
    #         continue
    #     for model in data[CO2EMISS_CUM][exp]:
    #         if (model not in data[CO2EMISS_CUM][EXPS['HIST']]):
    #             continue
    #         co2emiss_cum_hist = data[CO2EMISS_CUM][EXPS['HIST']][model][-1]
    #         co2emiss_cum_rcp = data[CO2EMISS_CUM][exp][model]
    #         co2emiss_cum_rcp += co2emiss_cum_hist - co2emiss_cum_rcp[0]

    # # Calculate ranges (with given model range)
    # model_range = E.get_config_option(modelconfig, 'plot', 'model_range', 0.9)

    # # Calculate historical and RCP range
    # tas_vs_co2emiss_cum_rcp = []
    # for exp in EXPS.values():
    #     if (exp == EXPS['1PCTCO2']):
    #         continue
    #     valid_models_tas = []
    #     valid_models_co2emiss_cum = []
    #     for model in data[VARS['TAS']][exp]:
    #         if (model not in data[CO2EMISS_CUM][exp]):
    #             continue
    #         valid_models_tas.append(data[VARS['TAS']][exp][model])
    #         valid_models_co2emiss_cum.append(data[CO2EMISS_CUM][exp][model])
    #     valid_models_tas = np.stack(valid_models_tas, axis=-1)
    #     valid_models_co2emiss_cum = np.stack(valid_models_co2emiss_cum,
    #                                          axis=-1)

    #     # Get desired range for every time step
    #     for time_idx in xrange(len(valid_models_tas)):
    #         time_slice = valid_models_tas[time_idx]
    #         mean = np.mean(time_slice)
    #         err = model_range * (np.amax(time_slice) - np.amin(time_slice)) * \
    #             0.5
    #         for idx in xrange(len(time_slice)):
    #             tas_val = time_slice[idx]
    #             if (tas_val < mean-err or tas_val > mean+err):
    #                 continue
    #             co2emiss_cum_val = valid_models_co2emiss_cum[time_idx][idx]
    #             tas_vs_co2emiss_cum_rcp.append(np.array([co2emiss_cum_val,
    #                                                      tas_val]))

    # tas_vs_co2emiss_cum_rcp = np.array(tas_vs_co2emiss_cum_rcp)
    # hist_rcp_range = scipy.spatial.ConvexHull(tas_vs_co2emiss_cum_rcp)

    # # Calculate 1pctCO2 range
    # tas_vs_co2emiss_cum_1pct = []
    # for exp in EXPS.values():
    #     if (exp != EXPS['1PCTCO2']):
    #         continue
    #     valid_models_tas = []
    #     valid_models_co2emiss_cum = []
    #     for model in data[VARS['TAS']][exp]:
    #         if (model not in data[CO2EMISS_CUM][exp]):
    #             continue
    #         valid_models_tas.append(data[VARS['TAS']][exp][model])
    #         valid_models_co2emiss_cum.append(data[CO2EMISS_CUM][exp][model])
    #     valid_models_tas = np.stack(valid_models_tas, axis=-1)
    #     valid_models_co2emiss_cum = np.stack(valid_models_co2emiss_cum,
    #                                          axis=-1)

    #     # Get desired range for every time step
    #     for time_idx in xrange(len(valid_models_tas)):
    #         time_slice = valid_models_tas[time_idx]
    #         mean = np.mean(time_slice)
    #         err = model_range * (np.amax(time_slice) - np.amin(time_slice)) * \
    #             0.5
    #         for idx in xrange(len(time_slice)):
    #             tas_val = time_slice[idx]
    #             if (tas_val < mean-err or tas_val > mean+err):
    #                 continue
    #             co2emiss_cum_val = valid_models_co2emiss_cum[time_idx][idx]
    #             tas_vs_co2emiss_cum_1pct.append(np.array([co2emiss_cum_val,
    #                                                       tas_val]))
    # tas_vs_co2emiss_cum_1pct = np.array(tas_vs_co2emiss_cum_1pct)
    # onepctco2_range = scipy.spatial.ConvexHull(tas_vs_co2emiss_cum_1pct)

    # ###########################################################################
    # # Plot data
    # ###########################################################################
    # if (write_plots):
    #     E.ensure_directory(plot_dir)
    #     style_file = E.get_path_to_mpl_style('small_font.mplstyle')
    #     plt.style.use(style_file)
    #     fig, ax = plt.subplots()
    #     handles = []

    #     # Plot 1pctco2 range
    #     alph_1pctco2_range = 0.3
    #     c_1pctco2_range = 'black'
    #     ax.fill(tas_vs_co2emiss_cum_1pct[:, 0][onepctco2_range.vertices],
    #             tas_vs_co2emiss_cum_1pct[:, 1][onepctco2_range.vertices],
    #             color=c_1pctco2_range, alpha=alph_1pctco2_range, linewidth=0.0)

    #     # Plot historical and rcp range
    #     alph_hist_rcp_range = 0.3
    #     c_hist_rcp_range = 'red'
    #     ax.fill(tas_vs_co2emiss_cum_rcp[:, 0][hist_rcp_range.vertices],
    #             tas_vs_co2emiss_cum_rcp[:, 1][hist_rcp_range.vertices],
    #             color=c_hist_rcp_range, alpha=alph_hist_rcp_range,
    #             linewidth=0.0)

    #     # Years which should be plotted
    #     text_years = {EXPS['1PCTCO2']: [],
    #                   EXPS['HIST']: [1890, 1950, 1980, 2000, 2010],
    #                   EXPS['RCP26']: [2030, 2050, 2100],
    #                   EXPS['RCP45']: [2030, 2050, 2100],
    #                   EXPS['RCP60']: [2050, 2100],
    #                   EXPS['RCP85']: [2050, 2100]}
    #     text_shifts = {EXPS['1PCTCO2']: [],
    #                    EXPS['HIST']: [(20, -0.1), (-63, 0.14), (20, -0.12),
    #                                   (-160, 0.1), (-160, 0.1)],
    #                    EXPS['RCP26']: [(-170, 0.1), (-150, 0.1), (0, -0.18)],
    #                    EXPS['RCP45']: [(20, -0.1), (30, -0.08), (20, -0.1)],
    #                    EXPS['RCP60']: [(20, -0.1), (40, 0)],
    #                    EXPS['RCP85']: [(-170, 0.1), (40, 0)]}

    #     # Plot all the selected models for all experiments
    #     for exp in EXPS.values():
    #         c = E.get_config_option(modelconfig, exp, 'color', 'black')
    #         lw = E.get_config_option(modelconfig, exp, 'linewidth', 1)
    #         m = E.get_config_option(modelconfig, exp, 'marker', 'o')
    #         lab = E.get_config_option(modelconfig, exp, 'label', 'Unknowm')
    #         ax.plot(data[CO2EMISS_CUM][exp][MM_MEAN],
    #                 data[VARS['TAS']][exp][MM_MEAN],
    #                 color=c, linewidth=lw, marker=m,
    #                 markerfacecolor=c, markeredgecolor=c)

    #         # Append historical and RCPs to legend
    #         if (exp != EXPS['1PCTCO2']):
    #             handles.append(mlines.Line2D([], [], color=c, label=lab,
    #                                          linewidth=lw))

    #         # Add years
    #         for idx in xrange(len(text_years[exp])):
    #             txt_year = text_years[exp][idx]
    #             try:
    #                 x_idx = np.where(
    #                     decades[CO2EMISS_CUM][exp][MM_MEAN] == txt_year)[0]
    #                 y_idx = np.where(
    #                     decades[VARS['TAS']][exp][MM_MEAN] == txt_year)[0]
    #                 x_pos = data[CO2EMISS_CUM][exp][MM_MEAN][x_idx][0] + \
    #                     text_shifts[exp][idx][0]
    #                 y_pos = data[VARS['TAS']][exp][MM_MEAN][y_idx][0] + \
    #                     text_shifts[exp][idx][1]
    #                 ax.text(x_pos, y_pos, '{:d}'.format(txt_year), ha='left',
    #                         va='center', color=c, fontsize='smaller')
    #             except IndexError:
    #                 pass

    #     # Append historical and RCP range to legend
    #     handles.append(mpatches.Patch(color=c_hist_rcp_range,
    #                                   alpha=alph_hist_rcp_range,
    #                                   linewidth=0.0, label='RCP range'))

    #     # Append 1pctCO2 to legend
    #     c = E.get_config_option(modelconfig, EXPS['1PCTCO2'], 'color', 'black')
    #     lw = E.get_config_option(modelconfig, EXPS['1PCTCO2'], 'linewidth', 1)
    #     m = E.get_config_option(modelconfig, EXPS['1PCTCO2'], 'marker', 'o')
    #     lab = E.get_config_option(modelconfig, EXPS['1PCTCO2'], 'label',
    #                               EXPS['1PCTCO2'])
    #     handles.append(mlines.Line2D([], [], color=c, label=lab, linewidth=lw))
    #     handles.append(mpatches.Patch(color=c_1pctco2_range,
    #                                   alpha=alph_1pctco2_range,
    #                                   linewidth=0.0, label=lab+' range'))

    #     # Create legend
    #     legend = ax.legend(handles=handles, loc='lower right',
    #                        fontsize='smaller', ncol=2)

    #     # General plot appearance
    #     sec = 'plot'
    #     xlim_left = E.get_config_option(modelconfig, sec, 'xlim_left',
    #                                     0.0)
    #     xlim_right = E.get_config_option(modelconfig, sec, 'xlim_right',
    #                                      2500.0)
    #     ylim_top = E.get_config_option(modelconfig, sec, 'ylim_top',
    #                                    5.0)
    #     ax.set_xlim(left=xlim_left, right=xlim_right)
    #     ax.set_ylim(top=ylim_top)
    #     ax.set_xlabel(r"Cumulative total anthropogenic CO$_2$ emissions " +
    #                   "from 1870 [GtC]")
    #     ax.set_ylabel(r"Temperature anomaly relative to 1861-1880 [$^\circ$C]")

    #     # Save file
    #     filename = "TCRE" + "." + plot_file_type
    #     filepath = os.path.join(plot_dir, filename)
    #     info("Creating {0}".format(filepath), verbosity, 1)
    #     fig.savefig(filepath, additional_artists=[legend],
    #                 bbox_inches='tight', orientation='landscape')

    #     plt.close()


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
