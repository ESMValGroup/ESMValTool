#!/usr/bin/env python
# -*- coding: utf-8 -*-


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
import iris
import numpy as np
from scipy import stats
import scipy.signal as scisi
import matplotlib.pyplot as plt
from esmvaltool.diag_scripts.shared import (
    run_diagnostic, select_metadata,
    variables_available, plot)

logger = logging.getLogger(os.path.basename(__file__))


def set_axx_deangelis4(axx, ylen, ytickstrs, x_obs, dx_obs):
    """Axis settings for deangelis 4."""
    axx.set_xlabel(r'clr-dSWA/dPW (% kg$^{-1}$ m$^2$)')
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
    axx.set_xlabel(r'clr-dSWA/dPW (% kg$^{-1}$ m$^2$)')
    axx.set_ylabel(r'clr-dSWA/dT (W m$^{-2}$ K$^{-1}$)')
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
    filepath = os.path.join(cfg['plot_dir'], dataset_name + 'fig3a.' +
                            cfg['output_file_type'])

    x_reg_prw = np.linspace(0.0, 65, 2)
    y_reg_prw = reg_prw.slope * x_reg_prw + reg_prw.intercept
    text_reg_prw = dataset_name + ', dSWA/dPW ={:.2f}'.format(reg_prw.slope)

    fig, axx = plt.subplots(figsize=(8, 8))
    axx.plot(x_reg_prw, y_reg_prw, linestyle='solid', color='k',
             label=text_reg_prw)
    axx.plot(data["x"], data["ypic"], linestyle='none',
             color='k', marker='d')

    ccc = 0.0
    for kobs in reg_obs.keys():
        y_reg_obs = (reg_obs[kobs].slope * x_reg_prw +
                     reg_obs[kobs].intercept)
        text_reg_obs = ('CERES-EBAF/' + kobs +
                        ', dSWA/dPW ={:.2f}'.format(reg_obs[kobs].slope))
        axx.plot(x_reg_prw, y_reg_obs, linestyle='solid',
                 color=(0.25 * ccc, 1.0 - 0.25 * ccc, 0.7),
                 label=text_reg_obs)
        axx.plot(data["x"], (data["yobs"])[kobs], linestyle='none',
                 color=(0.25 * ccc, 1.0 - 0.25 * ccc, 0.7), marker='<')
        ccc = ccc + 1.0
        if ccc > 4:
            ccc = 0.5

    axx.set_title(' ')
    axx.set_xlabel(r'Column precipitable water, PW (kg m$^{-2}$)')
    axx.set_ylabel('Normalized clr-SWA(%)')
    axx.set_xlim([9, 61])
    axx.set_xticks(np.linspace(10, 60, 6))
    axx.set_ylim([16.5, 26.5])
    axx.set_yticks(np.linspace(17, 26, 10))
    axx.legend(loc=1)

    fig.tight_layout()
    fig.savefig(filepath, dpi=300)
    plt.close()

    # def plot_deangelis_fig3b(cfg, m_drsnst_dts, m_drsnstcs_dprws,
    # m_drsnstcs_dprws_obs,
    # s_drsnstcs_dprws, s_drsnstcs_dprws_obs):


def plot_deangelis_fig3b(cfg, data_model, reg_prw_obs):
    """Plot DeAngelis Fig. 3a and 4."""
    # Fig 4

    # Set dictionary modelkey vs number from DeAngelis et al., 2015
    # all models which are in the paper, other dictonaries contain only
    # the currently selected ones
    model_nr = dict([('ACCESS1-0', [1, 2]), ('ACCESS1-3', [2, 2]),
                     ('bcc-csm1-1', [3, 5]),
                     ('bcc-csm1-1-m', [4, 5]), ('CanESM2', [5, 2]),
                     ('CCSM4', [6, 5]), ('CNRM-CM5', [7, 6]),
                     ('CNRM-CM5-2', [8, 6]), ('GFDL-CM3', [9, 3]),
                     ('GFDL-ESM2G', [10, 3]),
                     ('GFDL-ESM2M', [11, 3]), ('GISS-E2-H', [12, 4]),
                     ('GISS-E2-R', [13, 4]),
                     ('HadGEM2-ES', [14, 1]), ('inmcm4', [15, 5]),
                     ('IPSL-CM5A-LR', [16, 7]),
                     ('IPSL-CM5A-MR', [17, 7]), ('IPSL-CM5B-LR', [18, 7]),
                     ('MIROC-ESM', [19, 1]), ('MIROC5', [20, 1]),
                     ('MPI-ESM-LR', [21, 1]),
                     ('MPI-ESM-MR', [22, 1]), ('MPI-ESM-P', [23, 1]),
                     ('MRI-CGCM3', [24, 5]), ('NorESM1-M', [25, 5]),
                     ('default', [0, 8])])

    model_scheme_name = dict([(1, ['Correlated-k-distribution, N >= 20',
                                   'dodgerblue']),
                              (2, ['Correlated-k-distribution, 10 < N < 20',
                                   'limegreen']),
                              (3, ['Pseudo-k-distribution, N >= 20',
                                   'gold']),
                              (4, ['Pseudo-k-distribution, 10 < N < 20',
                                   'darkorange']),
                              (5, ['7-band parameterization (N = 7)',
                                   'crimson']),
                              (6, ['Pade approximants, higher res.',
                                   'slategrey']),
                              (7, ['Pade approximants, lower res.',
                                   'silver']),
                              (8, ['unknown',
                                   'xkcd:pale purple'])])

    # Plot data
    fig, axx = plt.subplots(figsize=(8, 8))

    mdrsnstdts = np.zeros((len(data_model.keys()), 3))
    # mdrsnstcsdprws = np.zeros(len(data_model.keys()))
    # sdrsnstcsdprws = np.zeros(len(data_model.keys()))

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
        axx.plot(mdrsnstdts[iii, 1],
                 mdrsnstdts[iii, 0],
                 marker=(plot.get_dataset_style(modelkey))['mark'],
                 color=(plot.get_dataset_style(modelkey))['color'],
                 markerfacecolor=(plot.get_dataset_style(modelkey))
                 ['facecolor'], linestyle='none',
                 markersize=10, markeredgewidth=2.0, label=modelkey)
        # axx.text(mdrsnstdts[iii, 1] - 0.005, mdrsnstdts[iii, 0] - 0.5 * 0.01,
        #          modelkey, color='k')

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
    fig.savefig(os.path.join(cfg['plot_dir'], 'fig3b.' +
                             cfg['output_file_type']), dpi=300)
    plt.close()

    # Fig 4
    ytickstrs = []

    fig, axx = plt.subplots(figsize=(8, 8))

    used_schemes = []

    for iii, jjj in enumerate(np.argsort(mdrsnstdts[:, 1])):
        modelkey = list(data_model.keys())[jjj]
        if modelkey not in model_nr.keys():
            modelkey = 'default'
        if (model_nr[modelkey])[1] not in used_schemes:
            axx.fill([mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2]],
                     [iii - 0.3, iii - 0.3, iii + 0.3, iii + 0.3],
                     color=(model_scheme_name[(model_nr[modelkey])[1]])[1],
                     label=(model_scheme_name[(model_nr[modelkey])[1]])[0])
            used_schemes.append((model_nr[modelkey])[1])
        else:
            axx.fill([mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] + 2.0 * mdrsnstdts[jjj, 2],
                      mdrsnstdts[jjj, 1] - 2.0 * mdrsnstdts[jjj, 2]],
                     [iii - 0.3, iii - 0.3, iii + 0.3, iii + 0.3],
                     color=(model_scheme_name[(model_nr[modelkey])[1]])[1])

        # axx.text(mdrsnstdts[jjj, 1] - 0.002, iii - 0.55 * 0.3,
        #          str((model_nr[modelkey])[0]), color='k')
        ytickstrs.append(list(data_model.keys())[jjj])

    axx = set_axx_deangelis4(axx, len(data_model.keys()),
                             ytickstrs,
                             prw["mean"], prw["diff"])

    fig.tight_layout()
    fig.savefig(os.path.join(cfg['plot_dir'], 'fig4.' +
                             cfg['output_file_type']), dpi=300)
    plt.close()


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

    if not variables_available(cfg, ['tas', 'rsnstcs',
                                     'rsnstcsnorm', 'prw']):
        raise ValueError("This diagnostic needs 'tas', " +
                         "'rsnstcs', 'rsnstcsnorm', and 'prw' variables")

    ###########################################################################
    # Read data
    ###########################################################################

    grid_pw = {}
    grid_pw["x"] = np.arange(12.0, 59.0, 2, dtype=float)
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
        raise ValueError("This diagnostic expects one (or no) observation " +
                         "data set for rsnstcsnorm")

    ###########################################################################
    # Process data
    ###########################################################################

    # Observations (one rsnstcsnorm, X PRW data set, DeAngelis contains 3
    data_prw_obs = OrderedDict()
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

    data_model = substract_and_reg_deangelis(cfg, cubes, grid_pw,
                                             reg_prw_obs)

    plot_deangelis_fig3b(cfg, data_model, reg_prw_obs)

    ###########################################################################
    # Write data
    ###########################################################################


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
