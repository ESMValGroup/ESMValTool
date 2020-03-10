"""
calculates and compares the correlation between the turnover time of carbon and climate defined as the partial correlations with precipitation and temperature
"""
# place your module imports here:
import extraUtils as xu

# operating system manipulations (e.g. path constructions)
import os
import sys

# to manipulate iris cubes
import iris
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

# internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, select_metadata, run_diagnostic


def _load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube


## Classes and settings


class dot_dict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def _get_fig_config(__cfg):

    fcfg = {}
    fcfg['fill_val'] = np.nan
    fcfg['multimodel'] = False
    fcfg['correlation_method'] = 'pearson'

    # define the data and information for plotting ratios
    fcfg['ax_fs'] = 7.1
    fcfg['valrange_x'] = (-1, 1)
    fcfg['valrange_y'] = (-70, 90)
    fcfg['minPoints'] = 3
    fcfg['bandsize'] = 9.5
    fcfg['obs_label'] = 'Carvalhais2014'
    fcfg['gpp_threshold'] = 0  # gC m-2 yr -1

    fcfgL = list(fcfg.keys())
    for _fc in fcfgL:
        if __cfg.get(_fc) is not None:
            fcfg[_fc] = __cfg.get(_fc)
    _figSet = dot_dict(fcfg)
    return _figSet


# data and calculations


def _apply_common_mask(_dat1, _dat2, _dat3):
    '''
    apply a common mask to three arrays so that they have the same locations of all valid and invalid (non numeric) grid cells
    '''
    _dat1Mask = np.ma.getmask(np.ma.masked_invalid(_dat1))
    _dat2Mask = np.ma.getmask(np.ma.masked_invalid(_dat2))
    _dat3Mask = np.ma.getmask(np.ma.masked_invalid(_dat3))
    _valMaskA = 1 - (1 - _dat1Mask) * (1 - _dat2Mask) * (1 - _dat3Mask)
    _valMask = np.ma.nonzero(_valMaskA)
    _dat1[_valMask] = np.nan
    _dat2[_valMask] = np.nan
    _dat3[_valMask] = np.nan
    _dat1 = np.ma.masked_invalid(_dat1)
    _dat2 = np.ma.masked_invalid(_dat2)
    _dat3 = np.ma.masked_invalid(_dat3)
    return _dat1, _dat2, _dat3

def _apply_gpp_threshold(gpp_dat, fig_config):
    '''
    returns the input array with values below the threshold of gpp set to nan
    '''    
    # converting gC m-2 yr-1 to kgC m-2 s-1
    gpp_thres = fig_config.gpp_threshold / (
        86400.0 * 365 * 1000.)  
    gpp_dat = np.ma.masked_less(gpp_dat, gpp_thres).filled(fig_config.fill_value)
    return gpp_dat


def _get_obs_data(cfg):
    """
    Get and handle the observations of turnover time from Carvalhais 2014.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        dictionary with observation data with different variables as keys
    """
    if not cfg.get('obs_files'):
        raise ValueError('The observation files needs to be specified in the '
                         'recipe (see recipe description for details)')
    else:
        input_files = [
            os.path.join(cfg['auxiliary_data_dir'], obs_file)
            for obs_file in cfg.get('obs_files')
        ]
    all_data = {}
    varList = cfg.get('obs_variables')
    for _var in varList:
        variable_constraint = iris.Constraint(
            cube_func=(lambda c: c.var_name == _var))
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data[_var] = cube.data
        for coord in cube.coords():
            print(coord.name())
            all_data[coord.name()] = coord.points
    return all_data


def _get_zonal_correlation(cfg):
    """
    A diagnostic function to calculate the zonal correlation between ecosystem 
    carbon turnover time and climate.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic
    """
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    fcfg = _get_fig_config(cfg)
    allModdat = {}
    for key, value in my_files_dict.items():
        allModdat[key] = {}
        mod_coords = {}
        c_total = _load_variable(value, 'ctotal')
        gpp = _load_variable(value, 'gpp')
        pr = _load_variable(value, 'pr')
        tas = _load_variable(value, 'tas')
        _gpp = gpp.data
        _gpp = _apply_gpp_threshold(_gpp, fcfg)
        ctau = (c_total / _gpp) / (86400 * 365)
        for coord in gpp.coords():
            mod_coords[coord.name()] = coord.points
        zon_corr = _zonal_correlation(ctau.data, pr.data, tas.data,
                                      mod_coords['latitude'], fcfg)
        allModdat[key]['data'] = zon_corr
        allModdat[key]['latitude'] = mod_coords['latitude']
    allObsdat = _get_obs_data(cfg)
    _plot_zonal_correlation(allModdat, allObsdat, cfg)
    return 'zonal correlation diagnostic is complete'


def _partialCorr(C, _fcfg):
    d1 = C[:, 0]
    d2 = C[:, 1]
    d3 = C[:, 2]
    if _fcfg.correlation_method == 'pearson':
        r12 = stats.pearsonr(d1, d2)[0]
        r13 = stats.pearsonr(d1, d3)[0]
        r23 = stats.pearsonr(d2, d3)[0]
    elif _fcfg.correlation_method == 'spearman':
        r12 = stats.spearmanr(d1, d2)[0]
        r13 = stats.spearmanr(d1, d3)[0]
        r23 = stats.spearmanr(d2, d3)[0]
    else:
        sys.exit('set a valid correlation_method [pearson/spearman]')
    r123 = (r12 - r13 * r23) / np.sqrt((1 - r13**2) * (1 - r23**2))
    return r123


def _zonal_correlation(_tau, _pr, _tas, _lats, _fcfg):
    _latint = abs(_lats[1] - _lats[0])
    windowSize = int(_fcfg.bandsize / (_latint * 2))
    __dat = np.ones((np.shape(_tau)[0], 2)) * np.nan
    _tau, _pr, _tas = _apply_common_mask(_tau, _pr, _tas)
    minPoints = np.shape(_tau)[1] / 8
    for lat_index in range(len(__dat)):
        istart = max(0, lat_index - windowSize)
        iend = min(np.size(_lats), lat_index + windowSize + 1)
        _tauZone = _tau[istart:iend, :]
        _prZone = _pr[istart:iend, :]
        _tasZone = _tas[istart:iend, :]
        d1 = np.ma.masked_invalid(_tauZone).compressed().flatten()
        d2 = np.ma.masked_invalid(_prZone).compressed().flatten()
        d3 = np.ma.masked_invalid(_tasZone).compressed().flatten()
        nValids = sum(~np.isnan(d1 + d2 + d3))
        if nValids > minPoints:
            __dat[lat_index, 0] = _partialCorr(np.vstack((d1, d2, d3)).T, _fcfg)
            __dat[lat_index, 1] = _partialCorr(np.vstack((d1, d3, d2)).T, _fcfg)
    return __dat


def _get_multimodel_stats(r_multimodel):
    """
    returns the mean, low and high correlations of all models using the fisher's z transformation

    Arguments:
        r_multimodel - zonal correlation from the models in the column dimensions

    Returns:
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

    # get the lower bound of correlation using inverse of fisher's z transformation
    z_low = zmm_ens - zmm_ens_std
    r_low = (np.exp(2 * z_low) - 1) / (np.exp(2 * z_low) + 1)

    # get the upper bound of correlation using inverse of fisher's z transformation
    z_high = zmm_ens + zmm_ens_std
    r_hi = (np.exp(2 * z_high) - 1) / (np.exp(2 * z_high) + 1)
    return r_mean, r_low, r_hi


# Plotting functions


def _fix_axis(x_lab, fig_config, ax_fs=8, axlw=0.4, rem_list=['top', 'right']):
    """
    fixes the axis limits, labels and lines

    Arguments:
        x_lab - axis labels
        fig_config - figure configurations
        ax_fs - fontsize for axis and tick labels
        ax_lw - linewidth of axis lines
        rem_list - list of axis lines to remove

    Returns:
    """
    plt.xlim(fig_config.valrange_x[0], fig_config.valrange_x[1])
    plt.ylim(fig_config.valrange_y[0], fig_config.valrange_y[1])
    plt.axhline(y=0, lw=0.48, color='grey')
    plt.axvline(x=0, lw=0.48, color='grey')
    plt.xlabel(x_lab, fontsize=fig_config.ax_fs)
    ax = plt.gca()
    for loc, spine in ax.spines.items():
        if loc in rem_list:
            spine.set_position(('outward', 0))
            spine.set_linewidth(0.)
        else:
            spine.set_linewidth(axlw)
    return


def _plot_zonal_correlation(all_mod_dat, all_obs_dat, cfg):
    """
    makes the line plots of zonal correlations from all models 

    Arguments:
        cfg - nested dictionary of metadata
        all_mod_dat - dictionary of correlations from all models
        all_obs_dat - dictionary of correlations and ranges from observation

    Returns:
        string; makes some time-series plots
    """
    _fcfg = _get_fig_config(cfg)
    models = list(all_mod_dat.keys())
    nmodels = len(models)
    models = sorted(models, key=str.casefold)

    plt.figure(figsize=(5, 4))
    # tau-tas correlations
    sp1 = plt.subplot(1, 2, 1)
    x_lab = '$r_{\\tau-tas,pr}$'
    _fix_axis(x_lab, _fcfg)
    plt.ylabel('Latitude ($^\\circ N$)', fontsize=_fcfg.ax_fs, ma='center')
    # get the observations out of the dictionary
    lats_obs = all_obs_dat['latitude']
    r_tau_ctotal_tas = all_obs_dat['r_tau_ctotal_tas']
    r_tau_ctotal_tas_5 = all_obs_dat['r_tau_ctotal_tas_5']
    r_tau_ctotal_tas_95 = all_obs_dat['r_tau_ctotal_tas_95']
    # plot the correlations from observation
    sp1.plot(r_tau_ctotal_tas,
             lats_obs,
             color='k',
             lw=1.1,
             label='Observation')
    sp1.fill_betweenx(lats_obs,
                      r_tau_ctotal_tas_5,
                      r_tau_ctotal_tas_95,
                      facecolor='grey',
                      alpha=0.40)

    # tau-pr correlations
    sp2 = plt.subplot(1, 2, 2)
    x_lab = '$r_{\\tau-pr,tas}$'
    _fix_axis(x_lab, _fcfg)

    # get the observations out of the dictionary
    r_tau_ctotal_pr = all_obs_dat['r_tau_ctotal_pr']
    r_tau_ctotal_pr_5 = all_obs_dat['r_tau_ctotal_pr_5']
    r_tau_ctotal_pr_95 = all_obs_dat['r_tau_ctotal_pr_95']

    # plot the correlations from observation
    sp2.plot(r_tau_ctotal_pr,
             lats_obs,
             color='k',
             lw=1.1,
             label=_fcfg.obs_label)
    sp2.fill_betweenx(lats_obs,
                      r_tau_ctotal_pr_5,
                      r_tau_ctotal_pr_95,
                      facecolor='grey',
                      alpha=0.40)

    #### PLOTTING for models

    # define arrays to store the zonal correlation of each model

    r_tau_pr_c_tas_all = np.ones((len(lats_obs), nmodels)) * np.nan
    r_tau_tas_c_pr_all = np.ones((len(lats_obs), nmodels)) * np.nan

    # loop over models and plot zonal correlations
    for row_m in range(nmodels):
        row_mod = models[row_m]
        r_mod = all_mod_dat[row_mod]['data']
        lats_mod = all_mod_dat[row_mod]['latitude']
        r_tau_tas_c_pr_mod = r_mod[:, 1]
        r_tau_tas_c_pr_all[:, row_m] = r_tau_tas_c_pr_mod
        sp1.plot(np.ma.masked_equal(r_tau_tas_c_pr_mod, np.nan),
                 lats_mod,
                 lw=0.3,
                 label=row_mod)
        r_tau_pr_c_tas_mod = r_mod[:, 0]
        r_tau_pr_c_tas_all[:, row_m] = r_tau_pr_c_tas_mod
        sp2.plot(np.ma.masked_equal(r_tau_pr_c_tas_mod, np.nan),
                 lats_mod,
                 lw=0.3,
                 label=row_mod)
    r_mmod, r_mmod_std_low, r_mmod_std_hi = _get_multimodel_stats(
        r_tau_tas_c_pr_all)
    sp1.plot(np.ma.masked_equal(r_mmod, np.nan),
             lats_mod,
             color='blue',
             ls='--',
             lw=1,
             label='Norm. Mean r')
    sp1.fill_betweenx(lats_mod,
                      np.ma.masked_equal(r_mmod_std_low, np.nan),
                      np.ma.masked_equal(r_mmod_std_hi, np.nan),
                      facecolor='#42d4f4',
                      alpha=0.25)

    r_mmod, r_mmod_std_low, r_mmod_std_hi = _get_multimodel_stats(
        r_tau_pr_c_tas_all)

    sp2.plot(np.ma.masked_equal(r_mmod, np.nan),
             lats_mod,
             color='blue',
             ls='--',
             lw=1,
             label='Norm. Mean r')
    sp2.fill_betweenx(lats_mod,
                      np.ma.masked_equal(r_mmod_std_low, np.nan),
                      np.ma.masked_equal(r_mmod_std_hi, np.nan),
                      facecolor='#42d4f4',
                      alpha=0.25)

    plt.gca().yaxis.set_label_position("right")
    # draw the legend
    leg = xu.draw_line_legend(ax_fs=_fcfg.ax_fs)
    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)

    # save and close the figure
    local_path = cfg['plot_dir']
    png_name = 'comparison_zonal_' + _fcfg.correlation_method + '_correlation_turnovertime_climate_' + _fcfg.obs_label + '.png'
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x, leg],
                dpi=450)
    plt.close()

    return 'Plotting of zonal correlation is complete'


if __name__ == '__main__':
    with run_diagnostic() as config:
        _get_zonal_correlation(config)
