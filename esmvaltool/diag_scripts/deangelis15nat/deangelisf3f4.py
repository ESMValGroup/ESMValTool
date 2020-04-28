""".

Calculates radiative constraint on hydrologic cycle intensification
following DeAngelis et al. (2015).

###############################################################################
testkw/deangelis3.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################


Description

-----------
    Calculates radiative constraint on hydrologic cycle intensification
    following DeAngelis et al. (2015).
    Creates figure 3b (with or without observations)
    Based on diag_scripts/climate_metrics/ecs.py by Manuel Schlund

Configuration options
---------------------

###############################################################################

"""


import logging
import os
from collections import OrderedDict
from pprint import pformat

import iris
import iris.coord_categorisation as cat
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as scisi
from scipy import stats

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, group_metadata,
                                            plot, run_diagnostic,
                                            select_metadata,
                                            variables_available)

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Run the diagnostic.

    Parameters :
    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    ###########################################################################
    # Read recipe data
    ###########################################################################
    # Get Input data
    input_data = (cfg['input_data'].values())
    required_vars = ('tas', 'rsnstcs', 'rsnstcsnorm', 'prw')

    if not variables_available(cfg, required_vars):
        raise ValueError("This diagnostic needs {required_vars} variables")

    available_exp = list(group_metadata(input_data, 'exp'))

    if 'abrupt-4xCO2' not in available_exp:
        if 'abrupt4xCO2' not in available_exp:
            raise ValueError("The diagnostic needs an experiment with " +
                             "4 times CO2.")

    if 'piControl' not in available_exp:
        raise ValueError("The diagnostic needs a pre industrial control " +
                         "experiment.")

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset and save annual means
    cubes = reform_data_iris_deangelis3b4(input_data)

    meas_tub_rsnstcsnorm = []
    meas_tub_prw = []
    for ctub in cubes:
        if ctub[2] == 'nomodel':
            if ctub[1] == 'prw':
                meas_tub_prw.append(ctub)
            if ctub[1] == 'rsnstcsnorm':
                meas_tub_rsnstcsnorm.append(ctub)

    if len(meas_tub_rsnstcsnorm) > 1:
        raise ValueError(
            "This diagnostic expects one (or no) observational "
            "dataset for rsnstcsnorm"
        )

    ###########################################################################
    # Process data
    ###########################################################################

    [grid_pw, reg_prw_obs] = set_grid_pw_reg_obs(cubes, meas_tub_rsnstcsnorm,
                                                 meas_tub_prw)

    data_model = substract_and_reg_deangelis(cfg, cubes, grid_pw,
                                             reg_prw_obs)

    plot_deangelis_fig3b4(cfg, data_model, reg_prw_obs)


def _get_sel_files_var(cfg, varnames):
    """Get filenames from cfg for all model mean and differen variables."""
    selection = []

    for var in varnames:
        for hlp in select_metadata(cfg['input_data'].values(), short_name=var):
            selection.append(hlp['filename'])

    return selection


def set_grid_pw_reg_obs(cubes, meas_tub_rsnstcsnorm, meas_tub_prw):
    """Set prw grid and calculate regression for observational data."""
    # Observations (one rsnstcsnorm, X PRW data set, DeAngelis contains 3
    data_prw_obs = OrderedDict()
    grid_pw = {}
    grid_pw["x"] = np.arange(12.0, 59.0, 2, dtype=float)
    grid_pw["yobs"] = OrderedDict()
    reg_prw_obs = OrderedDict()

    if np.min(np.array([len(meas_tub_rsnstcsnorm), len(meas_tub_prw)])) > 0:
        data_rsnstcsnorm_obs = cubes[meas_tub_rsnstcsnorm[0]]
        for kmeas_tub_prw in meas_tub_prw:
            data_prw_obs[kmeas_tub_prw[0]] = \
                cubes[kmeas_tub_prw]
            (grid_pw["yobs"])[kmeas_tub_prw[0]] = \
                make_grid_prw(grid_pw["x"],
                              data_prw_obs[kmeas_tub_prw[0]],
                              data_rsnstcsnorm_obs)
            reg_prw_obs[kmeas_tub_prw[0]] = \
                stats.linregress(grid_pw["x"],
                                 (grid_pw["yobs"])[kmeas_tub_prw[0]])
    else:
        logger.info('No observations, only model data used')

    return [grid_pw, reg_prw_obs]


def cube_to_save_matrix(var1, name):
    """Create cubes to prepare scatter plot data for saving to netCDF."""
    cubes = iris.cube.CubeList([iris.cube.Cube(var1,
                                               var_name=name['var_name'],
                                               long_name=name['long_name'],
                                               units=name['units'])])

    return cubes


def cube_to_save_vars(list_dict):
    """Create cubes to prepare bar plot data for saving to netCDF."""
    # cubes = iris.cube.CubeList()
    for iii, var in enumerate(list_dict["data"]):
        if iii == 0:
            cubes = iris.cube.CubeList([
                iris.cube.Cube(var,
                               var_name=list_dict["name"][iii]['var_name'],
                               long_name=list_dict["name"][iii]['long_name'],
                               units=list_dict["name"][iii]['units'])])
        else:
            cubes.append(
                iris.cube.Cube(var,
                               var_name=list_dict["name"][iii]['var_name'],
                               long_name=list_dict["name"][iii]['long_name'],
                               units=list_dict["name"][iii]['units']))

    return cubes


def cube_to_save_scatter(var1, var2, names):
    """Create cubes to prepare scatter plot data for saving to netCDF."""
    cubes = iris.cube.CubeList([iris.cube.Cube(var1,
                                               var_name=names['var_name1'],
                                               long_name=names['long_name1'],
                                               units=names['units1'])])
    cubes.append(iris.cube.Cube(var2, var_name=names['var_name2'],
                                long_name=names['long_name2'],
                                units=names['units2']))

    return cubes


def get_provenance_record(ancestor_files, caption, statistics,
                          domains, plot_type='other'):
    """Get Provenance record."""
    record = {
        'caption': caption,
        'statistics': statistics,
        'domains': domains,
        'plot_type': plot_type,
        'themes': ['phys'],
        'authors': [
            'weigel_katja',
        ],
        'references': [
            'deangelis15nat',
        ],
        'ancestors': ancestor_files,
    }
    return record


def set_axx_deangelis4(axx, ylen, ytickstrs, x_obs, dx_obs):
    """Axis settings for deangelis 4."""
    axx.set_xlabel(r'$\delta$rsnst / $\delta$prw (% kg$^{-1}$ m$^2$)')
    axx.set_xlim([0.02, 0.13])
    axx.set_xticks(np.linspace(0.02, 0.12, 6))
    axx.set_ylim([-0.5, ylen + 0.5])
    axx.set_yticks(np.linspace(0, ylen - 1, ylen))
    axx.legend(loc=2)
    axx.set_yticklabels(ytickstrs)

    # Observations
    if dx_obs != 0:
        axx.text(x_obs - dx_obs * 0.95, 11.8, 'Obs.', color='k')
    axx.vlines([x_obs - dx_obs, x_obs + dx_obs], -1, ylen + 1,
               colors='k', linestyle='dashed')
    if dx_obs != 0:
        axx.text(x_obs - dx_obs * 0.95, 11.8, 'Obs.', color='k')
        axx.arrow(x_obs - dx_obs, 11.5, 2 * dx_obs - 0.002, 0, head_width=0.2,
                  head_length=0.002, facecolor='k', edgecolor='k')
        axx.arrow(x_obs + dx_obs, 11.5, -2 * dx_obs + 0.002, 0,
                  head_width=0.2, head_length=0.002, facecolor='k',
                  edgecolor='k')

    return axx


def set_axx_deangelis3b(axx, x_obs, dx_obs):
    """Axis settings for deangelis 3b."""
    axx.set_xlabel(r'$\delta$rsnstcs / $\delta$prw (% kg$^{-1}$ m$^2$)')
    axx.set_ylabel(r'$\delta$rsnstcs / $\delta$tas (W m$^{-2}$ K$^{-1}$)')
    axx.set_xlim([0.0, 0.13])
    axx.set_xticks(np.linspace(0.0, 0.12, 7))
    axx.set_ylim([0.45, 1.15])
    axx.set_yticks(np.linspace(0.5, 1.1, 7))
    axx.legend(loc=2)

    # Observations
    axx.vlines([x_obs - dx_obs, x_obs + dx_obs], 0.4, 1.2, colors='k',
               linestyle='dashed')
    if dx_obs != 0:
        axx.text(x_obs - dx_obs * 0.95, 0.78, 'Obs.', color='k')
        axx.arrow(x_obs - dx_obs, 0.75, 2 * dx_obs - 0.002, 0,
                  head_width=0.02, head_length=0.002, facecolor='k',
                  edgecolor='k')
        axx.arrow(x_obs + dx_obs, 0.75, -2 * dx_obs + 0.002, 0,
                  head_width=0.02, head_length=0.002, facecolor='k',
                  edgecolor='k')

    return axx


def plot_deangelis_fig3a(cfg, dataset_name, data, reg_prw, reg_obs):
    """Plot DeAngelis Fig. 3a."""
    reg_prw_dict = {}
    reg_prw_dict["x"] = np.linspace(0.0, 65, 2)
    reg_prw_dict["y"] = reg_prw.slope * reg_prw_dict["x"] + reg_prw.intercept
    reg_prw_dict["text"] = dataset_name + \
        r', $\delta$rsnst / $\delta$prw ={:.2f}'.format(reg_prw.slope)

    fig, axx = plt.subplots(figsize=(8, 8))
    axx.plot(reg_prw_dict["x"], reg_prw_dict["y"], linestyle='solid',
             color='k',
             label=reg_prw_dict["text"])
    axx.plot(data["x"], data["ypic"], linestyle='none',
             color='k', marker='d')

    ccc = 0.0
    for kobs in reg_obs.keys():
        reg_prw_dict["y_obs"] = (reg_obs[kobs].slope * reg_prw_dict["x"] +
                                 reg_obs[kobs].intercept)
        reg_prw_dict["text_obs"] = ('CERES-EBAF/' + kobs +
                                    r', $\delta$rsnst / $\delta$' +
                                    'prw ={:.2f}'.format(reg_obs[kobs].slope))
        axx.plot(reg_prw_dict["x"], reg_prw_dict["y_obs"], linestyle='solid',
                 color=(0.25 * ccc, 1.0 - 0.25 * ccc, 0.7),
                 label=reg_prw_dict["text_obs"])
        axx.plot(data["x"], (data["yobs"])[kobs], linestyle='none',
                 color=(0.25 * ccc, 1.0 - 0.25 * ccc, 0.7), marker='<')
        ccc = ccc + 1.0
        if ccc > 4:
            ccc = 0.5

    axx.set_title(' ')
    axx.set_xlabel(r'Water Vapor Path (prw) (kg m$^{-2}$)')
    axx.set_ylabel('Normalized Netto Short Wave Radiation for clear syke (%)')
    axx.set_xlim([9, 61])
    axx.set_xticks(np.linspace(10, 60, 6))
    axx.set_ylim([16.5, 26.5])
    axx.set_yticks(np.linspace(17, 26, 10))
    axx.legend(loc=1)

    fig.tight_layout()
    fig.savefig(get_plot_filename('fig3a_' + dataset_name, cfg), dpi=300)
    plt.close()

    caption = 'Scatter plot and regression lines the between the ' + \
        'netto short wave radiation for clear skye normalized by ' + \
        'normalized by incoming solar flux (rsnstcsnorm) and the ' + \
        'Water Vapor Path (prw) in the pre-industrial climate.'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, ['prw', 'rsnstcsnorm']), caption, ['other'],
        ['global'])

    diagnostic_file = get_diagnostic_filename('fig3a_' + dataset_name, cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = {}
    list_dict["data"] = [data["x"], data["ypic"]]
    list_dict["name"] = [{'var_name': 'prw_' + dataset_name,
                          'long_name': 'Water Vapor Path ' + dataset_name,
                          'units': 'kg m-2'},
                         {'var_name': 'rsnstdtsnorm',
                          'long_name': 'Normalized Netto Short ' +
                                       'Wave Radiation fo' +
                                       'clear syke',
                          'units': 'percent'}]
    for kobs in reg_obs.keys():
        list_dict["data"].append((data["yobs"])[kobs])
        list_dict["name"].append({'var_name': 'prw_' + kobs,
                                  'long_name': 'Water Vapor Path ' +
                                               'CERES-EBAF/' + kobs,
                                  'units': 'kg m-2'})

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def plot_deangelis_fig4(cfg, data_model, mdrsnstdts, prw):
    """Plot DeAngelis Fig. 4."""
    # Set dictionary modelkey vs number from DeAngelis et al., 2015
    # all models which are in the paper, other dictonaries contain only
    # the currently selected ones
    model_dict = {
        'ACCESS1-0': ('Correlated-k-distribution, 10 < N < 20', 'limegreen'),
        'ACCESS1-3': ('Correlated-k-distribution, 10 < N < 20', 'limegreen'),
        'bcc-csm1-1': ('7-band parameterization (N = 7)', 'crimson'),
        'bcc-csm1-1-m': ('7-band parameterization (N = 7)', 'crimson'),
        'CanESM2': ('Correlated-k-distribution, 10 < N < 20', 'limegreen'),
        'CCSM4': ('7-band parameterization (N = 7)', 'crimson'),
        'CNRM-CM5': ('Pade approximants, higher res.', 'slategrey'),
        'CNRM-CM5-2': ('Pade approximants, higher res.', 'slategrey'),
        'GFDL-CM3': ('Pseudo-k-distribution, N >= 20', 'gold'),
        'GFDL-ESM2G': ('Pseudo-k-distribution, N >= 20', 'gold'),
        'GFDL-ESM2M': ('Pseudo-k-distribution, N >= 20', 'gold'),
        'GISS-E2-H': ('Pseudo-k-distribution, 10 < N < 20', 'darkorange'),
        'GISS-E2-R': ('Pseudo-k-distribution, 10 < N < 20', 'darkorange'),
        'HadGEM2-ES': ('Correlated-k-distribution, N >= 20', 'dodgerblue'),
        'inmcm4': ('7-band parameterization (N = 7)', 'crimson'),
        'IPSL-CM5A-LR': ('Pade approximants, lower res.', 'silver'),
        'IPSL-CM5A-MR': ('Pade approximants, lower res.', 'silver'),
        'IPSL-CM5B-LR': ('Pade approximants, lower res.', 'silver'),
        'MIROC-ESM': ('Correlated-k-distribution, N >= 20', 'dodgerblue'),
        'MIROC5': ('Correlated-k-distribution, N >= 20', 'dodgerblue'),
        'MPI-ESM-LR': ('Correlated-k-distribution, N >= 20', 'dodgerblue'),
        'MPI-ESM-MR': ('Correlated-k-distribution, N >= 20', 'dodgerblue'),
        'MPI-ESM-P': ('Correlated-k-distribution, N >= 20', 'dodgerblue'),
        'MRI-CGCM3': ('7-band parameterization (N = 7)', 'crimson'),
        'NorESM1-M': ('7-band parameterization (N = 7)', 'crimson'),
        'default': ('unknown', 'xkcd:pale purple'),
    }

    ytickstrs_and_schemes = {}
    ytickstrs_and_schemes["ytickstrs"] = []
    ytickstrs_and_schemes["schemes"] = []
    # used_schemes = []
    fig, axx = plt.subplots(figsize=(8, 8))

    for iii, jjj in enumerate(np.argsort(mdrsnstdts[:, 1])):
        modelkey = list(data_model.keys())[jjj]
        if modelkey not in model_dict.keys():
            modelkey = 'default'
        if (model_dict[modelkey])[0] not in ytickstrs_and_schemes["schemes"]:
            axx.fill([mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2]],
                     [iii - 0.3, iii - 0.3, iii + 0.3, iii + 0.3],
                     color=model_dict[modelkey][1],
                     label=model_dict[modelkey][0])
            ytickstrs_and_schemes["schemes"].append(model_dict[modelkey][0])
        else:
            axx.fill([mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2]],
                     [iii - 0.3, iii - 0.3, iii + 0.3, iii + 0.3],
                     color=model_dict[modelkey][1])

        ytickstrs_and_schemes["ytickstrs"].append(list(data_model.keys())[jjj])

    axx = set_axx_deangelis4(axx, len(data_model.keys()),
                             ytickstrs_and_schemes["ytickstrs"],
                             prw["mean"], prw["diff"])

    fig.tight_layout()

    fig.savefig(get_plot_filename('fig4', cfg), dpi=300)
    plt.close()

    caption = 'The relationship between the ratio of the change of ' + \
        'netto short wave radiation (rsnst) and the change of the ' + \
        'Water Vapor Path (prw) and characteristics of the ' + \
        'parameterization scheme for solar absorption by water vapour ' + \
        'in a cloud-free atmosphere, with colours for each model ' + \
        'referring to different types of parameterizations as described ' + \
        'in the key (N refers to the number of exponential terms ' + \
        'representing water vapour absorption). The width of horizontal ' + \
        'shading for models and the vertical dashed lines for ' + \
        'observations (Obs.) represent statistical uncertainties of ' + \
        'the ratio, as the 95% confidence interval (CI) of the regression ' + \
        'slope to the rsnst versus prw curve. For the observations ' + \
        'the minimum of the lower bounds of all CIs to the maximum of ' + \
        'the upper bounds of all CIs is shown.'

    provenance_record = get_provenance_record(_get_sel_files_var(cfg,
                                                                 ['prw',
                                                                  'rsnst']),
                                              caption, ['other'], ['global'])

    diagnostic_file = get_diagnostic_filename('fig4', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_matrix(mdrsnstdts, {'var_name': 'mdrsnstdts',
                                               'long_name': 'Change of Netto' +
                                                            'Short Wave ' +
                                                            'Radiation ' +
                                                            'with Change ' +
                                                            'of the ' +
                                                            'Water Vapor Path',
                                               'units': '% kg-1 m2'}),
              target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def plot_deangelis_fig3b4(cfg, data_model, reg_prw_obs):
    """Plot DeAngelis Fig. 3b and prepare 4."""
    # Fig 4
#    model_scheme_name = dict([(1, ['Correlated-k-distribution, N >= 20',
#                                   'dodgerblue']),
#                              (2, ['Correlated-k-distribution, 10 < N < 20',
#                                   'limegreen']),
#                              (3, ['Pseudo-k-distribution, N >= 20',
#                                   'gold']),
#                              (4, ['Pseudo-k-distribution, 10 < N < 20',
#                                   'darkorange']),
#                              (5, ['7-band parameterization (N = 7)',
#                                   'crimson']),
#                              (6, ['Pade approximants, higher res.',
#                                   'slategrey']),
#                              (7, ['Pade approximants, lower res.',
#                                   'silver']),
#                              (8, ['unknown',
#                                   'xkcd:pale purple'])])

    # Plot data
    fig, axx = plt.subplots(figsize=(8, 8))

    mdrsnstdts = np.zeros((len(data_model.keys()), 3))

    for iii, modelkey in enumerate(data_model.keys()):
        mdrsnstdts[iii, 0] = (list(data_model[modelkey]))[0]
        mdrsnstdts[iii, 1] = (list(data_model[modelkey]))[2]
        mdrsnstdts[iii, 2] = (list(data_model[modelkey]))[3]
        axx.fill([mdrsnstdts[iii, 1] - 2.0 * mdrsnstdts[iii, 2],
                  mdrsnstdts[iii, 1] + 2.0 * mdrsnstdts[iii, 2],
                  mdrsnstdts[iii, 1] + 2.0 * mdrsnstdts[iii, 2],
                  mdrsnstdts[iii, 1] - 2.0 * mdrsnstdts[iii, 2]],
                 [mdrsnstdts[iii, 0] - 0.01, mdrsnstdts[iii, 0] - 0.01,
                  mdrsnstdts[iii, 0] + 0.01, mdrsnstdts[iii, 0] + 0.01],
                 color=(0.8, 0.8, 0.8))
        style = (select_metadata(cfg['input_data'].values(),
                                 dataset=modelkey))[0]['project']
        style = plot.get_dataset_style(modelkey, style_file=style.lower())
        axx.plot(mdrsnstdts[iii, 1],
                 mdrsnstdts[iii, 0],
                 marker=style['mark'],
                 color=style['color'],
                 markerfacecolor=style['facecolor'], linestyle='none',
                 markersize=10, markeredgewidth=2.0, label=modelkey)

    prw = {}
    if reg_prw_obs:
        prw["min"] = np.zeros(len(reg_prw_obs.keys()))
        prw["max"] = np.zeros(len(reg_prw_obs.keys()))
        for iii, modelkey in enumerate(reg_prw_obs.keys()):
            (prw["min"])[iii] = reg_prw_obs[modelkey].slope - 2.0 * \
                reg_prw_obs[modelkey].stderr
            (prw["max"])[iii] = reg_prw_obs[modelkey].slope + 2.0 * \
                reg_prw_obs[modelkey].stderr
        prw["mean"] = (np.min(prw["min"]) + np.max(prw["max"])) / 2.0
        prw["diff"] = np.max(prw["max"]) - prw["mean"]
    else:
        prw["mean"] = 0.0
        prw["diff"] = 0.0

    axx = set_axx_deangelis3b(axx, prw["mean"], prw["diff"])

    # Regression line
    reg = stats.linregress(mdrsnstdts[:, 1], mdrsnstdts[:, 0])

    axx.plot(np.linspace(0.0, 0.15, 2), reg.slope *
             np.linspace(0.0, 0.15, 2) + reg.intercept,
             linestyle='solid', color='r',
             label='Fit (r ={:.2f})'.format(reg.rvalue))

    axx.legend(ncol=2, loc=2)
    fig.tight_layout()
    fig.savefig(get_plot_filename('fig3b', cfg), dpi=300)
    plt.close()

    caption = 'Scatter plot and regression line the between the ratio ' + \
        'of the change of ' + \
        'netto short wave radiation (rsnst) and the change of the ' + \
        'Water Vapor Path (prw) against the ratio of the change of ' + \
        'netto short wave radiation for clear skye (rsnstcs) and the ' + \
        'the change of surface temperature (tas).' + \
        'The width of horizontal ' + \
        'shading for models and the vertical dashed lines for ' + \
        'observations (Obs.) represent statistical uncertainties of ' + \
        'the ratio, as the 95% confidence interval (CI) of the regression ' + \
        'slope to the rsnst versus prw curve. For the observations ' + \
        'the minimum of the lower bounds of all CIs to the maximum of ' + \
        'the upper bounds of all CIs is shown.'

    provenance_record = get_provenance_record(_get_sel_files_var(cfg,
                                                                 ['prw',
                                                                  'tas',
                                                                  'rsnst',
                                                                  'rsnstcs']),
                                              caption, ['other'], ['global'])

    diagnostic_file = get_diagnostic_filename('fig3b', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_matrix(mdrsnstdts, {'var_name': 'mdrsnstdts',
                                               'long_name': 'Change of Netto' +
                                                            'Short Wave ' +
                                                            'Radiation ' +
                                                            'with Change ' +
                                                            'of the ' +
                                                            'Water Vapor Path',
                                               'units': '% kg-1 m2'}),
              target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)

    # Fig 4
    plot_deangelis_fig4(cfg, data_model, mdrsnstdts, prw)


def make_grid_prw(grid_pwx, data_prw_obs, data_rsnstcsnorm_obs):
    """Grid rsnstcsnorm based on prw grid."""
    gridded_rsnstcsnorm_obs = np.zeros(len(grid_pwx), dtype=float)

    for jjj, bincenter in enumerate(grid_pwx):
        index_obs = np.where((data_prw_obs >= bincenter - 1.0) &
                             (data_prw_obs < bincenter + 1.0))
        gridded_rsnstcsnorm_obs[jjj] = np.mean(data_rsnstcsnorm_obs[index_obs])

    return gridded_rsnstcsnorm_obs


def reform_data_iris_deangelis3b4(input_data):
    """Extract data from IRIS cubes and average or reformat them."""
    # Model data for 'tas', 'rsnstcs'
    cubes = {}
    for my_short_name in ['tas', 'rsnstcs']:
        # my_data: List of dictionaries
        my_data = select_metadata(input_data, short_name=my_short_name)
        # subdata: dictionary
        for subdata in my_data:
            cube = iris.load(subdata['filename'])[0]
            cat.add_year(cube, 'time', name='year')
            cube = cube.aggregated_by('year', iris.analysis.MEAN)
            experiment = subdata['exp']
            if experiment == 'abrupt-4xCO2':
                experiment = 'abrupt4xCO2'
            dataset = subdata['dataset']
            cubetuple = (dataset, my_short_name, experiment)
            if experiment == 'piControl':
                # DeAngelis use a 21 month running mean on piControl but the
                # full extend of 150 years abrupt4xCO2. I could not find out,
                # how they tread the edges, currently I just skip the mean for
                # the edges. This is not exacly the same as done in the paper,
                # small differences remain in extended data Fig 1,
                # but closer than other methods I
                # tried, e.g. skipping the edges.
                # For most data sets it would be also possible to
                # extend the piControl for 20 years, but then it would
                # not be centered means of piControl for each year of
                # abrupt4xCO2 any more.
                # cube_new = cube.rolling_window('time',iris.analysis.MEAN, 21)
                # endm10 = len(cube.coord('time').points) - 10
                # cube.data[10:endm10] = cube_new.data
                cube.data = scisi.savgol_filter(cube.data, 21, 1, axis=0)
            cubes[cubetuple] = cube.data

    # Model data and observations for 'rsnstcsnorm', and 'prw'
    for my_short_name in ['rsnstcsnorm', 'prw']:
        # my_data: List of dictionaries
        my_data = select_metadata(input_data, short_name=my_short_name)
        # subdata: dictionary
        for subdata in my_data:
            if 'exp' in subdata.keys():
                experiment = subdata['exp']
            else:
                experiment = 'nomodel'
            dataset = subdata['dataset']
            cubetuple = (dataset, my_short_name, experiment)
            if experiment in ['piControl', 'nomodel']:
                cube = iris.load(subdata['filename'])[0]
                total_len = len(cube.coord('time').points) * \
                    len(cube.coord('latitude').points) * \
                    len(cube.coord('longitude').points)
                data_new = np.reshape(cube.data, total_len)
                cubes[cubetuple] = data_new

    return cubes


def substract_and_reg_deangelis(cfg, cubes, grid_pw,
                                reg_prw_obs):
    """Substract piControl from abrupt4xCO2 for all models and variables."""
    data_model = OrderedDict()

    model_tub_tas_pi = []

    for tub in cubes.keys():
        if tub[2] == 'piControl':
            if tub[1] == 'tas':
                model_tub_tas_pi.append(tub)

    # Models
    for model_run in model_tub_tas_pi:
        # Substract piControl experiment from abrupt4xCO2 experiment
        data_prw_pic = (cubes[(model_run[0], 'prw', model_run[2])])
        data_rsnstcsnorm_pic = (cubes[(model_run[0],
                                       'rsnstcsnorm', model_run[2])])
        data_tas = (cubes[(model_run[0],
                           model_run[1], 'abrupt4xCO2')]).data - \
            (cubes[model_run])
        data_rsnstcs = (cubes[(model_run[0], 'rsnstcs', 'abrupt4xCO2')]) - \
            (cubes[(model_run[0], 'rsnstcs', model_run[2])])

        grid_pw["ypic"] = make_grid_prw(grid_pw["x"], data_prw_pic,
                                        data_rsnstcsnorm_pic)

        reg6 = stats.linregress(data_tas, data_rsnstcs)

        reg_prw = stats.linregress(grid_pw["x"], grid_pw["ypic"])

        data_model[model_run[0]] = [reg6.slope, reg6.stderr, reg_prw.slope,
                                    reg_prw.stderr]

        plot_deangelis_fig3a(cfg, model_run[0],
                             grid_pw, reg_prw,
                             reg_prw_obs)

    return data_model


###############################################################################
# Setup diagnostic
###############################################################################

# Experiments
# 'piControl' and 'abrupt4xCO2'

    ###########################################################################
    # Write data
    ###########################################################################
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
