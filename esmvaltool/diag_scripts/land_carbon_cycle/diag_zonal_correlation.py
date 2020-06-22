"""
Evaluate the correlation between turnover time of carbon and climate.

Use partial correlations between turnover time with precipitation and
temperature
"""

import sys
import iris
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

import esmvaltool.diag_scripts.land_carbon_cycle.plot_utils as plut
from esmvaltool.diag_scripts.land_carbon_cycle.shared import (
    _apply_common_mask,
    _get_obs_data_zonal,
    _load_variable,
    _remove_invalid,
)
from esmvaltool.diag_scripts.land_carbon_cycle.provenance import (
    _get_ancestor_files,
    _get_provenance_record,
)


def _get_fig_config(diag_config):
    """
    Get figure setting and configurations.

    default settings of the figure, and replace default with
    runtime settings from recipe

    Argument:
    --------
        diag_config - nested dictionary of metadata

    Return:
    ------
        a dictionary of settings
    """
    fig_config = {
        'fill_value': np.nan,
        'correlation_method': 'pearson',
        'min_points_frac': 0.125,
        # define the data and information for plotting ratios
        'ax_fs': 7.1,
        'valrange_x': (-1, 1),
        'valrange_y': (-70, 90),
        'bandsize': 9.5,
        'gpp_threshold': 0.01
    }
    fig_config.update(diag_config.get('fig_config'))
    return fig_config


def partial_corr(dat_columns, fig_config):
    """
    Calculate the linear partial correlation.

    The correlation between variables in the first and second column of
    dat_columns is controlled for the covariation with that in the third
    column.

    Argument:
    --------
        dat_columns - an array with different variables in different columns
        fig_config - configuration with correlation_method. Uses the scipy
        stats module to calculate correlation using either pearsons linear
        (http://tiny.cc/pearsonr) or spearmans rank (http://tiny.cc/spearmanr)
        correlation coefficients.

    Return:
    ------
        r123 - correlation between variables 1 and 2 controlled for 3
    """
    dat_x = dat_columns[:, 0]
    dat_y = dat_columns[:, 1]
    dat_z = dat_columns[:, 2]
    if fig_config['correlation_method'] == 'pearson':
        r12 = stats.pearsonr(dat_x, dat_y)[0]
        r13 = stats.pearsonr(dat_x, dat_z)[0]
        r23 = stats.pearsonr(dat_y, dat_z)[0]
    elif fig_config['correlation_method'] == 'spearman':
        r12 = stats.spearmanr(dat_x, dat_y)[0]
        r13 = stats.spearmanr(dat_x, dat_z)[0]
        r23 = stats.spearmanr(dat_y, dat_z)[0]
    else:
        sys.exit('set a valid correlation_method [pearson/spearman]')
    # calculate the partial correlation coefficient as,
    # rxy,z = (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
    # https://en.wikipedia.org/wiki/Partial_correlation
    r123 = (r12 - r13 * r23) / np.sqrt((1 - r13**2) * (1 - r23**2))
    return r123


def _calc_zonal_correlation(dat_tau, dat_pr, dat_tas, dat_lats, fig_config):
    """
    Calculate zonal partial correlations for sliding windows.

    Argument:
    --------
        dat_tau - data of global tau
        dat_pr - precipitation
        dat_tas - air temperature
        dat_lats - latitude of the given model
        fig_config - figure/diagnostic configurations

    Return:
    ------
        corr_dat zonal correlations
    """
    # get the interval of latitude and create array for partial correlation
    lat_int = abs(dat_lats[1] - dat_lats[0])
    corr_dat = np.ones((np.shape(dat_tau)[0], 2)) * np.nan

    # get the size of the sliding window based on the bandsize in degrees
    window_size = round(fig_config['bandsize'] / (lat_int * 2.))

    dat_tau, dat_pr, dat_tas = _apply_common_mask(dat_tau, dat_pr, dat_tas)
    # minimum 1/8 of the given window has valid data points
    min_points = np.shape(dat_tau)[1] * fig_config['min_points_frac']
    for lat_index in range(len(corr_dat)):
        istart = np.int(max(0, lat_index - window_size))
        iend = np.int(min(np.size(dat_lats), lat_index + window_size + 1))
        dat_tau_zone = dat_tau[istart:iend, :]
        dat_pr_zone = dat_pr[istart:iend, :]
        dat_tas_zone = dat_tas[istart:iend, :]
        dat_x = np.ma.masked_invalid(dat_tau_zone).compressed().flatten()
        dat_y = np.ma.masked_invalid(dat_pr_zone).compressed().flatten()
        dat_z = np.ma.masked_invalid(dat_tas_zone).compressed().flatten()
        num_valid_points = sum(~np.isnan(dat_x + dat_y + dat_z))
        if num_valid_points > min_points:
            corr_dat[lat_index, 1] = partial_corr(
                np.vstack((dat_x, dat_y, dat_z)).T, fig_config)
            corr_dat[lat_index, 0] = partial_corr(
                np.vstack((dat_x, dat_z, dat_y)).T, fig_config)
    return corr_dat


def _get_multimodel_stats(r_multimodel):
    """
    Compute mean, low and high correlations of all models.

    Uses the fisher's z transformation

    Argument:
    --------
        r_multimodel - zonal correlation from the models in the column
        dimensions

    Return:
    ------
        mean, mean - std, and mean + std correlations
    """
    # set the threshold of correlation to avoid infinities
    r_multimodel[r_multimodel > 0.99] = 0.99
    r_multimodel[r_multimodel < -0.99] = -0.99

    # z tranform the correlation
    z_multimodel = 0.5 * (np.log(1 + r_multimodel) - np.log(1 - r_multimodel))
    z_multimodel[np.isinf(z_multimodel)] = np.nan
    zmm_ens = np.nanmean(z_multimodel, axis=1)
    zmm_ens_std = np.nanstd(z_multimodel, axis=1)

    # get the mean correlation using inverse of fisher's z transformation
    r_mean = (np.exp(2 * zmm_ens) - 1) / (np.exp(2 * zmm_ens) + 1)

    # get the lower bound of correlation using inverse of fisher's z
    # transformation
    z_low = zmm_ens - zmm_ens_std
    r_low = (np.exp(2 * z_low) - 1) / (np.exp(2 * z_low) + 1)

    # get the upper bound of correlation using inverse of fisher's z
    # transformation
    z_high = zmm_ens + zmm_ens_std
    r_hi = (np.exp(2 * z_high) - 1) / (np.exp(2 * z_high) + 1)
    return r_mean, r_low, r_hi


def _fix_axis(x_lab, fig_config, axlw=0.4, rem_list=('top', 'right')):
    """
    Fix the axis limits, labels and lines.

    Argument:
    --------
        x_lab - axis labels
        fig_config - figure configurations
        ax_fs - fontsize for axis and tick labels
        ax_lw - linewidth of axis lines
        rem_list - list of axis lines to remove
    """
    plt.xlim(fig_config['valrange_x'][0], fig_config['valrange_x'][1])
    plt.ylim(fig_config['valrange_y'][0], fig_config['valrange_y'][1])
    plt.axhline(y=0, lw=0.48, color='grey')
    plt.axvline(x=0, lw=0.48, color='grey')
    plt.xlabel(x_lab, fontsize=fig_config['ax_fs'])
    _ax = plt.gca()
    for loc, spine in _ax.spines.items():
        if loc in rem_list:
            spine.set_position(('outward', 0))
            spine.set_linewidth(0.)
        else:
            spine.set_linewidth(axlw)


def _plot_zonal_correlation(plot_path, zonal_correlation_mod,
                            zonal_correlation_obs, diag_config):
    """
    Make the line plots of zonal correlations from all models.

    Argument:
    --------
        diag_config - nested dictionary of metadata
        zonal_correlation_mod - dictionary of correlations from all models
        zonal_correlation_obs - dictionary of correlations and ranges from
        observation
    """
    fig_config = _get_fig_config(diag_config)
    models = list(zonal_correlation_mod.keys())
    nmodels = len(models)
    models = sorted(models, key=str.casefold)
    multimodel_stats = 'MultiModelMedian MultiModelMean'.split()
    for _mm in multimodel_stats:
        if _mm in models:
            models.append(models.pop(models.index(_mm)))

    plt.figure(figsize=(5, 4))
    # tau-tas correlations
    sp1 = plt.subplot(1, 2, 1)

    # get the observations out of the dictionary
    lats_obs = zonal_correlation_obs['latitude']
    obs_var = diag_config.get('obs_variable')[0]
    r_tau_ctotal_tas = zonal_correlation_obs[obs_var]
    r_tau_ctotal_tas_5 = zonal_correlation_obs[obs_var + '_5']
    r_tau_ctotal_tas_95 = zonal_correlation_obs[obs_var + '_95']
    # plot the correlations from observation

    _fix_axis(obs_var, fig_config)
    plt.ylabel('{name}\n({unit})'.format(name=lats_obs.long_name,
                                         unit=lats_obs.units),
               fontsize=fig_config['ax_fs'],
               ma='center')

    sp1.plot(r_tau_ctotal_tas.data,
             lats_obs.points,
             color='k',
             lw=1.1,
             label=diag_config['obs_info']['source_label'])
    sp1.fill_betweenx(lats_obs.points,
                      r_tau_ctotal_tas_5.data,
                      r_tau_ctotal_tas_95.data,
                      facecolor='grey',
                      alpha=0.40)

    # tau-pr correlations
    sp2 = plt.subplot(1, 2, 2)

    # get the observations out of the dictionary
    obs_var = diag_config.get('obs_variable')[1]
    r_tau_ctotal_pr = zonal_correlation_obs[obs_var]
    r_tau_ctotal_pr_5 = zonal_correlation_obs[obs_var + '_5']
    r_tau_ctotal_pr_95 = zonal_correlation_obs[obs_var + '_95']
    _fix_axis(obs_var, fig_config)

    # plot the correlations from observation
    sp2.plot(r_tau_ctotal_pr.data,
             lats_obs.points,
             color='k',
             lw=1.1,
             label=diag_config['obs_info']['source_label'])
    sp2.fill_betweenx(lats_obs.points,
                      r_tau_ctotal_pr_5.data,
                      r_tau_ctotal_pr_95.data,
                      facecolor='grey',
                      alpha=0.40)

    # PLOTTING for models

    # loop over models and plot zonal correlations
    for row_m in range(nmodels):
        row_mod = models[row_m]
        r_mod = zonal_correlation_mod[row_mod]['data']
        lats_mod = zonal_correlation_mod[row_mod]['latitude']
        r_tau_tas_c_pr_mod = r_mod[:, 0]
        r_tau_pr_c_tas_mod = r_mod[:, 1]
        if row_mod in ['MultiModelMedian', 'MultiModelMean']:
            sp1.plot(np.ma.masked_equal(r_tau_tas_c_pr_mod, np.nan),
                     lats_mod.points,
                     lw=1.1,
                     color='blue',
                     label=row_mod)
            sp2.plot(np.ma.masked_equal(r_tau_pr_c_tas_mod, np.nan),
                     lats_mod.points,
                     lw=1.1,
                     color='blue',
                     label=row_mod)
        else:
            sp1.plot(np.ma.masked_equal(r_tau_tas_c_pr_mod, np.nan),
                     lats_mod.points,
                     lw=0.3,
                     label=row_mod)
            sp2.plot(np.ma.masked_equal(r_tau_pr_c_tas_mod, np.nan),
                     lats_mod.points,
                     lw=0.3,
                     label=row_mod)

    # normalized mean correlations from model

    # remove the multimodel estimates
    models = list(zonal_correlation_mod.keys())
    for _mm in multimodel_stats:
        if _mm in models:
            models.remove(_mm)

    nmodels = len(models)

    r_tau_pr_c_tas_all = np.ones((len(lats_obs.points), nmodels)) * np.nan
    r_tau_tas_c_pr_all = np.ones((len(lats_obs.points), nmodels)) * np.nan
    for row_m in range(nmodels):
        row_mod = models[row_m]
        r_mod = zonal_correlation_mod[row_mod]['data']
        lats_mod = zonal_correlation_mod[row_mod]['latitude']
        r_tau_tas_c_pr_all[:, row_m] = r_mod[:, 0]
        r_tau_pr_c_tas_all[:, row_m] = r_mod[:, 1]

    r_mmod, r_mmod_std_low, r_mmod_std_hi = _get_multimodel_stats(
        r_tau_tas_c_pr_all)
    sp1.plot(np.ma.masked_equal(r_mmod, np.nan),
             lats_mod.points,
             color='red',
             ls='--',
             lw=1,
             label='Norm. Mean r')
    sp1.fill_betweenx(lats_mod.points,
                      np.ma.masked_equal(r_mmod_std_low, np.nan),
                      np.ma.masked_equal(r_mmod_std_hi, np.nan),
                      facecolor='#42d4f4',
                      alpha=0.25)

    r_mmod, r_mmod_std_low, r_mmod_std_hi = _get_multimodel_stats(
        r_tau_pr_c_tas_all)

    sp2.plot(np.ma.masked_equal(r_mmod, np.nan),
             lats_mod.points,
             color='red',
             ls='--',
             lw=1,
             label='Norm. Mean r')
    sp2.fill_betweenx(lats_mod.points,
                      np.ma.masked_equal(r_mmod_std_low, np.nan),
                      np.ma.masked_equal(r_mmod_std_hi, np.nan),
                      facecolor='#42d4f4',
                      alpha=0.25)

    plt.gca().yaxis.set_label_position("right")

    # draw the legend
    leg = plut.draw_line_legend(ax_fs=fig_config['ax_fs'])

    plut.save_figure(plot_path, _extr_art=[leg])
    plt.close()


def main(diag_config):
    """
    Diagnostic to evaluate zonal correlation between turnover time and climate.

    Argument:
    --------
        diag_config - nested dictionary of metadata
    """
    model_data_dict = group_metadata(diag_config['input_data'].values(),
                                     'dataset')
    fig_config = _get_fig_config(diag_config)
    zonal_correlation_mod = {}
    for model_name, model_dataset in model_data_dict.items():
        zonal_correlation_mod[model_name] = {}
        mod_coords = {}
        ctotal = _load_variable(model_dataset, 'ctotal')
        gpp = _load_variable(model_dataset, 'gpp')
        precip = _load_variable(model_dataset, 'pr')
        tas = _load_variable(model_dataset, 'tas')
        tau_ctotal = (ctotal / gpp)
        tau_ctotal.convert_units('yr')
        # set the attributes
        tau_ctotal.var_name = 'tau_ctotal'
        for coord in gpp.coords():
            mod_coords[coord.name()] = coord

        _tau_dat = _remove_invalid(tau_ctotal.data, fill_value=np.nan)
        _precip_dat = _remove_invalid(precip.data, fill_value=np.nan)
        _tas_dat = _remove_invalid(tas.data, fill_value=np.nan)
        zon_corr = _calc_zonal_correlation(_tau_dat, _precip_dat, _tas_dat,
                                           mod_coords['latitude'].points,
                                           fig_config)
        zonal_correlation_mod[model_name]['data'] = zon_corr
        zonal_correlation_mod[model_name]['latitude'] = mod_coords['latitude']
    zonal_correlation_obs = _get_obs_data_zonal(diag_config)

    base_name = '{title}_{corr}_{source_label}_{grid_label}z'.format(
        title='r_tau_ctotal_climate',
        corr=fig_config['correlation_method'],
        source_label=diag_config['obs_info']['source_label'],
        grid_label=diag_config['obs_info']['grid_label'])

    provenance_record = _get_provenance_record(
        "Comparison of latitudinal (zonal) variations of pearson"
        " correlation between turnover time and climate: turnover"
        " time and precipitation, controlled for temperature"
        " (left) and vice-versa (right). Reproduces figures 2c"
        " and 2d in Carvalhais et al. (2014).", ['corr', 'perc'], ['zonal'],
        _get_ancestor_files(diag_config, 'tau_ctotal'))

    if diag_config['write_netcdf']:
        model_cubes = [
            c for c in zonal_correlation_mod.values()
            if isinstance(c, iris.cube.Cube)
        ]
        obs_cubes = [
            c for c in zonal_correlation_obs.values()
            if isinstance(c, iris.cube.Cube)
        ]
        netcdf_path = get_diagnostic_filename(base_name, diag_config)
        save_cubes = iris.cube.CubeList(model_cubes + obs_cubes)
        iris.save(save_cubes, netcdf_path)

        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(netcdf_path, provenance_record)

    if diag_config['write_plots']:
        plot_path = get_plot_filename(base_name, diag_config)
        _plot_zonal_correlation(plot_path, zonal_correlation_mod,
                                zonal_correlation_obs, diag_config)
        provenance_record['plot_file'] = plot_path

        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
