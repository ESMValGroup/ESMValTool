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

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename, group_metadata, select_metadata, run_diagnostic, extract_variables
from esmvalcore.preprocessor import area_statistics

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
    fcfg['correlation_method'] = 'pearson'
    fcfg['ax_fs'] = 7.1
    fcfg['valrange_x'] = (2, 1000)
    fcfg['valrange_y'] = (-70, 90)
    fcfg['minPoints'] = 0
    fcfg['bandsize'] = 1.986
    fcfg['obs_label'] = 'Carvalhais2014'
    fcfg['gpp_threshold'] = 10  #gC m-2 yr -1
    fcfgL = list(fcfg.keys())
    for _fc in fcfgL:
        if __cfg.get(_fc) != None:
            fcfg[_fc] = __cfg.get(_fc)
    _figSet = dotDict(fcfg)
    _figSet.ax_fs = 9
    return _figSet


# calculation and data functions
def _apply_common_mask(_dat1, _dat2):
    '''
    returns a mask array with 1 where all three data arrays have valid numeric values and zero elsewhere
    '''
    _dat1Mask = np.ma.getmask(np.ma.masked_invalid(_dat1))
    _dat2Mask = np.ma.getmask(np.ma.masked_invalid(_dat2))
    _valMaskA = 1 - (1 - _dat1Mask) * (1 - _dat2Mask)
    _valMask = np.ma.nonzero(_valMaskA)
    _dat1[_valMask] = np.nan
    _dat2[_valMask] = np.nan
    _dat1 = np.ma.masked_invalid(_dat1)
    _dat2 = np.ma.masked_invalid(_dat2)
    return _dat1, _dat2


def _apply_gpp_threshold(_gppDat, _fcfg):
    gpp_thres = _fcfg.gpp_threshold / (
        86400.0 * 365 * 1000.)  # converting gC m-2 yr-1 to kgC m-2 s-1
    _gppDat = np.ma.masked_less(_gppDat, gpp_thres).filled(_fcfg.fill_val)
    return _gppDat


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
    print(varList)
    for _var in varList:
        variable_constraint = iris.Constraint(
            cube_func=(lambda c: c.var_name == _var))
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data[_var] = cube.data
    for coord in cube.coords():
        all_data[coord.name()] = coord.points
    return all_data


def _get_zonal_tau(cfg):
    """
    A diagnostic function to calculate the total carbon stock from the list of variables
    passed to the diagnostic script.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    """

    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')

    fcfg = _get_fig_config(cfg)
    zonal_tau_mod = {}
    global_gpp_mod = {}
    global_ctotal_mod = {}
    for key, value in my_files_dict.items():
        mod_coords = {}
        zonal_tau_mod[key] = {}
        ctotal = _load_variable(value, 'ctotal')
        gpp = _load_variable(value, 'gpp')
        gppDat = xu.remove_invalid(gpp.data, _fill_val=fcfg._fill_val)
        ctotalDat = xu.remove_invalid(ctotal.data, _fill_val=fcfg._fill_val)
        for coord in gpp.coords():
            mod_coords[coord.name()] = coord.points
        zonal_tau_mod[key]['data'] = _zonal_tau(gppDat, ctotalDat,
                                                mod_coords['latitude'], fcfg)
        zonal_tau_mod[key]['latitude'] = mod_coords['latitude']
        print(key, value)

        global_ctotal_mod[key] = ctotalDat
        global_gpp_mod[key] = gppDat

    # get the multimodel median GPP and ctotal and calculate zonal tau from multimodel median
    mm_ctotal = xu.remove_invalid(np.nanmedian(np.array(
        [_ctotal for _ctotal in global_ctotal_mod.values()]),
                                               axis=0),
                                  _fill_val=fcfg.fill_val)
    mm_gpp = xu.remove_invalid(np.nanmedian(np.array(
        [_gpp for _gpp in global_gpp_mod.values()]),
                                            axis=0),
                               _fill_val=fcfg.fill_val)
    zonal_tau_mod['zmultimodel'] = {}
    zonal_tau_mod['zmultimodel']['data'] = _zonal_tau(mm_gpp, mm_ctotal,
                                                      mod_coords['latitude'],
                                                      fcfg)
    zonal_tau_mod['zmultimodel']['latitude'] = mod_coords['latitude']
    # get the observation of zonal tau
    zonal_tau_obs = _get_obs_data(cfg)
    # plot the figure
    _plot_zonal_tau(zonal_tau_mod, zonal_tau_obs, cfg)
    return 'done plotting the zonal turnover time'


def _load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube


def _zonal_tau(_datgpp, _datcTotal, _lats, _fcfg):
    __dat = np.ones_like(_lats) * np.nan
    _latint = abs(_lats[1] - _lats[0])
    windowSize = int(_fcfg.bandsize / (_latint * 2))

    _datgpp = xu.remove_invalid(_datgpp, _fill_val=_fcfg.fill_val)
    _datcTotal = xu.remove_invalid(_datcTotal, _fill_val=_fcfg.fill_val)
    # _datgpp = _apply_gpp_threshold(_datgpp,_fcfg)
    _datgppC, _datcTotalC = _apply_common_mask(_datgpp, _datcTotal)
    for li in range(len(__dat)):
        istart = max(0, li - windowSize)
        iend = min(np.size(_lats), li + windowSize + 1)
        _datgppZone = _datgpp[istart:iend, :]  #* _areaZone
        _datcTotalZone = _datcTotal[istart:iend, :]  #* _areaZone
        nValids = np.nansum(1 -
                            np.ma.getmask(np.ma.masked_invalid(_datgppZone)))
        if nValids > _fcfg.minPoints:
            __dat[li] = (np.nansum(_datcTotalZone) /
                         np.nansum(_datgppZone)) / (86400 * 365)
        else:
            __dat[li] = np.nan
    return __dat


# Plotting functions


def _plot_zonal_tau(all_mod_dat, all_obs_dat, cfg):
    """
    makes the maps of variables 

    Arguments:
        cfg - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot

    Returns:
        string; makes some time-series plots

    Note: this function is private; remove the '_'
    so you can make it public.
    """
    _fcfg = _get_fig_config(cfg)
    models = list(all_mod_dat.keys())
    models = sorted(models, key=str.casefold)
    nmodels = len(models)

    lats_obs = all_obs_dat['latitude']
    tau_obs = all_obs_dat['tau_ctotal']
    tau_obs_5 = all_obs_dat['tau_ctotal_5']
    tau_obs_95 = all_obs_dat['tau_ctotal_95']

    plt.figure(figsize=(3, 5))

    sp0 = plt.subplot(1, 1, 1)
    sp0.plot(tau_obs, lats_obs, color='k', lw=1.5, label=_fcfg.obs_label)
    sp0.fill_betweenx(lats_obs,
                      tau_obs_5,
                      tau_obs_95,
                      facecolor='grey',
                      alpha=0.40)

    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_mod_tau = all_mod_dat[row_mod]['data']
        lats_mod = all_mod_dat[row_mod]['latitude']
        if row_mod == 'zmultimodel':
            sp0.plot(np.ma.masked_equal(dat_mod_tau, np.nan),
                     lats_mod,
                     lw=1.5,
                     color='blue',
                     label='Multimodel')
        else:
            sp0.plot(np.ma.masked_equal(dat_mod_tau, np.nan),
                     lats_mod,
                     lw=0.5,
                     label=row_mod)

    leg = xu.draw_line_legend(ax_fs=_fcfg.ax_fs)

    plt.gca().set_xscale('log')
    plt.xlim(_fcfg.valrange_x[0], _fcfg.valrange_x[1])
    plt.ylim(_fcfg.valrange_y[0], _fcfg.valrange_y[1])
    plt.axhline(y=0, lw=0.48, color='grey')
    x_lab = '$\\tau$'
    plt.xlabel(x_lab, fontsize=_fcfg.ax_fs)
    plt.ylabel('Latitude ($^\\circ N$)', fontsize=_fcfg.ax_fs, ma='center')
    xu.rem_axLine(['top', 'right'])

    local_path = cfg['plot_dir']
    png_name = 'comparison_zonal_turnovertime_' + _fcfg.obs_label + '.png'
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[leg],
                dpi=450)
    plt.close()

    return 'Plotting complete for zonal turnover'


# main
if __name__ == '__main__':
    with run_diagnostic() as config:
        _get_zonal_tau(config)
