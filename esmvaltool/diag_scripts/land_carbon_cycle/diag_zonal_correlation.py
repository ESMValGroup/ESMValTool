"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared.supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
# place your module imports here:
import extraUtils as xu

# operating system manipulations (e.g. path constructions)
import os

# to manipulate iris cubes
import iris
import matplotlib.pyplot as plt
import numpy as np

# internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename, group_metadata, select_metadata, run_diagnostic, extract_variables
from esmvalcore.preprocessor import area_statistics


def _load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube


### user-defined functions for calculating partial correlations

## Classes and settings


class dotDict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def _get_fig_config(__cfg):

    fcfg = {}
    fcfg['fill_val'] = np.nan
    fcfg['multimodel'] = False
    fcfg['correlation_method'] = 'spearman'

    # define the data and information for plotting ratios
    fcfg['ax_fs'] = 7.1
    fcfg['valrange_x'] = (-1, 1)
    fcfg['valrange_y'] = (-70, 90)
    fcfg['minPoints'] = 3
    fcfg['bandsize'] = 9.5
    fcfg['outlierPerc'] = 0
    fcfg['obs_label'] = 'Carvalhais2014'
    fcfg['gpp_threshold'] = 0  # gC m-2 yr -1
    fcfgL = list(fcfg.keys())
    for _fc in fcfgL:
        if __cfg.get(_fc) != None:
            fcfg[_fc] = __cfg.get(_fc)
    _figSet = dotDict(fcfg)
    return _figSet


# data and calculations


def _apply_common_mask(_dat1, _dat2, _dat3):
    '''
    returns a mask array with 1 where all three data arrays have valid numeric
    values and zero elsewhere
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


def _apply_gpp_threshold(_gppDat, _fcfg):
    gpp_thres = _fcfg.gpp_threshold / (
        86400.0 * 365 * 1000.)  # converting gC m-2 yr-1 to kgC m-2 s-1
    _gppDat = np.ma.masked_less(_gppDat, gpp_thres).filled(_fcfg.fill_val)
    return _gppDat


def _fisher_z(_rdat):
    _zdat = 0.5 * (np.log(1 + _rdat) - np.log(1 - _rdat))
    return _zdat


def _get_obs_data(cfg):
    """Get and handle the observations of turnover time from Carvalhais 2014.

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


def _inverse_fisher_z(_zdat):
    _rdat = (np.exp(2 * _zdat) - 1) / (np.exp(2 * _zdat) + 1)
    return _rdat


def _partialCorr(C, _fcfg):
    d1 = C[:, 0]
    d2 = C[:, 1]
    d3 = C[:, 2]
    if d1.size > _fcfg.minPoints:
        d1out = _percentile_based_outlier(d1, threshold=_fcfg.outlierPerc)
        d1[d1out] = np.nan
        d2out = _percentile_based_outlier(d2, threshold=_fcfg.outlierPerc)
        d2[d2out] = np.nan
        d3out = _percentile_based_outlier(d3, threshold=_fcfg.outlierPerc)
        d3[d3out] = np.nan
        d1, d2, d3 = _apply_common_mask(d1, d2, d3)
        d1 = np.ma.masked_invalid(d1).compressed().flatten()
        d2 = np.ma.masked_invalid(d2).compressed().flatten()
        d3 = np.ma.masked_invalid(d3).compressed().flatten()
        if _fcfg.correlation_method == 'pearson':
            r12, p = xu.calc_pearson_r(d1, d2, outlierPerc=_fcfg.outlierPerc)
            r13, p = xu.calc_pearson_r(d1, d3, outlierPerc=_fcfg.outlierPerc)
            r23, p = xu.calc_pearson_r(d2, d3, outlierPerc=_fcfg.outlierPerc)
        elif _fcfg.correlation_method == 'spearman':
            r12, p = xu.calc_spearman_r(d1, d2, outlierPerc=_fcfg.outlierPerc)
            r13, p = xu.calc_spearman_r(d1, d3, outlierPerc=_fcfg.outlierPerc)
            r23, p = xu.calc_spearman_r(d2, d3, outlierPerc=_fcfg.outlierPerc)
        else:
            print('set a valid correLation_method [pearson/spearman]')
            exit
        r123 = (r12 - r13 * r23) / np.sqrt((1 - r13**2) * (1 - r23**2))
    else:
        r123 = np.nan
    return r123


def _percentile_based_outlier(data, threshold=2):
    diff = threshold / 2.0
    minval, maxval = np.percentile(data, [diff, 100 - diff])
    return (data < minval) | (data > maxval)


def _zonal_correlation(_dat, _pr, _tas, _lats, _fcfg):
    _latint = abs(_lats[1] - _lats[0])
    windowSize = int(_fcfg.bandsize / (_latint * 2))
    __dat = np.zeros((np.shape(_dat)[0], 2))
    _dat, _pr, _tas = _apply_common_mask(_dat, _pr, _tas)
    for li in range(len(__dat)):
        istart = max(0, li - windowSize)
        iend = min(np.size(_lats), li + windowSize + 1)
        _datZone = _dat[istart:iend, :]
        _prZone = _pr[istart:iend, :]
        _tasZone = _tas[istart:iend, :]
        _datZoneC = np.ma.masked_invalid(_datZone).compressed().flatten()
        _prZoneC = np.ma.masked_invalid(_prZone).compressed().flatten()
        _tasZoneC = np.ma.masked_invalid(_tasZone).compressed().flatten()
        pc_vpt = _partialCorr(
            np.column_stack((_datZoneC, _prZoneC, _tasZoneC)), _fcfg)
        pc_vtp = _partialCorr(
            np.column_stack((_datZoneC, _tasZoneC, _prZoneC)), _fcfg)
        __dat[li, 0] = pc_vpt
        __dat[li, 1] = pc_vtp
    return __dat


# Plotting functions


def _fix_axis(x_lab, _fcfg, ax_fs=8, axlw=0.4, rem_list=['top', 'right']):
    plt.xlim(_fcfg.valrange_x[0], _fcfg.valrange_x[1])
    plt.ylim(_fcfg.valrange_y[0], _fcfg.valrange_y[1])
    plt.axhline(y=0, lw=0.48, color='grey')
    plt.axvline(x=0, lw=0.48, color='grey')
    plt.xlabel(x_lab, fontsize=_fcfg.ax_fs)
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

    rAll_pt = np.zeros((len(lats_obs), nmodels))
    rAll_tp = np.zeros((len(lats_obs), nmodels))
    # loop over models and plot zonal correlations
    for row_m in range(nmodels):
        row_mod = models[row_m]
        mod_dat_row = all_mod_dat[row_mod]['data']
        lats_mod = all_mod_dat[row_mod]['latitude']
        r_mod_vtp = mod_dat_row[:, 1]
        rAll_tp[:, row_m] = r_mod_vtp
        sp1.plot(np.ma.masked_equal(r_mod_vtp, np.nan),
                 lats_mod,
                 lw=0.3,
                 label=row_mod)
        r_mod_vpt = mod_dat_row[:, 0]
        rAll_pt[:, row_m] = r_mod_vpt
        sp2.plot(np.ma.masked_equal(r_mod_vpt, np.nan),
                 lats_mod,
                 lw=0.3,
                 label=row_mod)

    # get the normalized mean zonal correlation of all models for tau-tas,pr
    zAll_tp = _fisher_z(rAll_tp)  # do the fisher's transformation of r
    zAll_tp[np.isinf(zAll_tp)] = np.nan
    # mean of fisher's z
    zmm_ens = np.nanmean(zAll_tp, axis=1)
    zmm_ens_std = np.nanstd(zAll_tp, axis=1)
    r_mmod = inverse__fisher_z(zmm_ens)
    # get the uncertainty ranges
    r_mmod_std_low = inverse__fisher_z(zmm_ens - zmm_ens_std)
    r_mmod_std_hi = inverse__fisher_z(zmm_ens + zmm_ens_std)
    # plot the normalized mean and uncertainty
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

    # get the normalized mean zonal correlation of all models for tau-pr,tas
    zAll_pt = _fisher_z(rAll_pt)
    zAll_pt[np.isinf(zAll_pt)] = np.nan
    zmm_ens = np.nanmean(zAll_pt, axis=1)
    zmm_ens_std = np.nanstd(zAll_pt, axis=1)
    r_mmod = inverse__fisher_z(zmm_ens)
    r_mmod_std_low = inverse__fisher_z(zmm_ens - zmm_ens_std)
    r_mmod_std_hi = inverse__fisher_z(zmm_ens + zmm_ens_std)

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
    # generate output path and save

    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)
    local_path = cfg['plot_dir']
    png_name = 'comparison_zonal_' + _fcfg.correlation_method + '_correlation_turnovertime_climate_' + _fcfg.obs_label + '.png'
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x, leg],
                dpi=450)
    plt.close()

    return 'Plotting complete'


if __name__ == '__main__':
    with run_diagnostic() as config:
        _get_zonal_correlation(config)
