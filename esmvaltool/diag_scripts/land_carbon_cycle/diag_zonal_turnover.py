'''
compares the zonal turnover time from the models with observation from
Carvalhais et al., 2014
'''

# operating system manipulations (e.g. path constructions)
import os

# to manipulate iris cubes
import math
import iris
import matplotlib.pyplot as plt
import numpy as np
# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata
from esmvaltool.diag_scripts.shared import run_diagnostic

# helper functions
import extraUtils as xu
from shared import _load_variable, _get_obs_data_zonal


# Classes and settings
def _get_fig_config(diag_config):
    '''
    get the default settings of the figure, and replace default with
    runtime settings from recipe

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        a dot dictionary of settings
    '''
    fig_config = {
        'fill_value': np.nan,
        'multimodel': False,
        'min_points_frac': 0.02,
        'ax_fs': 7.1,
        'valrange_x': (2, 1000),
        'valrange_y': (-70, 90),
        'bandsize': 1.86,
        'obs_label': 'Carvalhais2014',
        'gpp_threshold': 10,  # gC m-2 yr -1
    }
    fig_config.update(diag_config.get('fig_config'))
    return fig_config


# calculation and data functions
def _apply_common_mask(dat_1, dat_2):
    '''
    apply a common mask of valid grid cells across two input arrays
    '''
    dat_1_mask = np.ma.getmask(np.ma.masked_invalid(dat_1))
    dat_2_mask = np.ma.getmask(np.ma.masked_invalid(dat_2))
    _val_mask_a = 1 - (1 - dat_1_mask) * (1 - dat_2_mask)
    _val_mask = np.ma.nonzero(_val_mask_a)
    dat_1[_val_mask] = np.nan
    dat_2[_val_mask] = np.nan
    dat_1 = np.ma.masked_invalid(dat_1)
    dat_2 = np.ma.masked_invalid(dat_2)
    return dat_1, dat_2


def _get_zonal_tau(diag_config):
    '''
    A diagnostic function to calculate the total carbon stock from the list of
    variables passed to the diagnostic script.

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    '''

    my_files_dict = group_metadata(diag_config['input_data'].values(),
                                   'dataset')

    fig_config = _get_fig_config(diag_config)
    zonal_tau_mod = {}
    for key, value in my_files_dict.items():
        zonal_tau_mod[key] = {}
        ctotal = _load_variable(value, 'ctotal')
        gpp = _load_variable(value, 'gpp')
        zonal_tau_mod[key]['data'] = _calc_zonal_tau(gpp, ctotal, fig_config)

    zonal_tau_obs = _get_obs_data_zonal(diag_config)
    # plot the figure
    _plot_zonal_tau(zonal_tau_mod, zonal_tau_obs, diag_config)
    return 'done plotting the zonal turnover time'


def _calc_zonal_tau(gpp, ctotal, fig_config):
    '''
    calculate zonal turnover time

    Arguments:
        gpp - cube of global gpp
        ctotal - cube of total carbon content
        fig_config - figure/diagnostic configurations

    Returns:
        zonal_tau - zonal turnover time of carbon
    '''
    # get the interval of latitude and create array for partial correlation
    dat_lats = gpp.coord('latitude').points
    lat_int = abs(dat_lats[1] - dat_lats[0])
    # get the size of the sliding window based on the bandsize in degrees
    window_size = max(2, math.ceil(fig_config['bandsize'] / lat_int))

    gpp_zs = gpp.collapsed('longitude', iris.analysis.SUM)
    ctotal_zs = ctotal.collapsed('longitude', iris.analysis.SUM)
    gpp_z = gpp_zs.rolling_window('latitude',
                                  iris.analysis.SUM, window_size)
    ctotal_z = ctotal_zs.rolling_window('latitude',
                                        iris.analysis.SUM, window_size)

    zonal_tau = ctotal_z / gpp_z
    zonal_tau.data = zonal_tau.core_data() / (86400 * 365)

    return zonal_tau


# Plotting functions
def _plot_zonal_tau(all_mod_dat, all_obs_dat, diag_config):
    '''
    makes the maps of variables

    Arguments:
        diag_config - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot

    Returns:
        string; makes some time-series plots

    Note: this function is private; remove the '_'
    so you can make it public.
    '''
    fig_config = _get_fig_config(diag_config)
    models = list(all_mod_dat.keys())
    models = sorted(models, key=str.casefold)
    nmodels = len(models)

    lats_obs = all_obs_dat['latitude']
    tau_obs = all_obs_dat['tau_ctotal']
    tau_obs_5 = all_obs_dat['tau_ctotal_5']
    tau_obs_95 = all_obs_dat['tau_ctotal_95']

    plt.figure(figsize=(3, 5))

    sp0 = plt.subplot(1, 1, 1)
    sp0.plot(tau_obs, lats_obs, color='k', lw=1.5,
             label=fig_config['obs_label'])
    sp0.fill_betweenx(lats_obs,
                      tau_obs_5,
                      tau_obs_95,
                      facecolor='grey',
                      alpha=0.40)

    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_mod_tau = all_mod_dat[row_mod]['data']
        lats_mod = dat_mod_tau.coord('latitude')
        if row_mod == 'MultiModelMedian':
            sp0.plot(dat_mod_tau.data,
                     lats_mod.points,
                     lw=1.5,
                     color='blue',
                     label='Multimodel')
        else:
            sp0.plot(dat_mod_tau.data,
                     lats_mod.points,
                     lw=0.5,
                     label=row_mod)

    leg = xu.draw_line_legend(ax_fs=fig_config['ax_fs'])

    plt.gca().set_xscale('log')
    plt.xlim(fig_config['valrange_x'][0], fig_config['valrange_x'][1])
    plt.ylim(fig_config['valrange_y'][0], fig_config['valrange_y'][1])
    plt.axhline(y=0, lw=0.48, color='grey')
    x_lab = '$\\tau$'
    plt.xlabel(x_lab, fontsize=fig_config['ax_fs'])
    plt.ylabel('Latitude ($^\\circ N$)',
               fontsize=fig_config['ax_fs'],
               ma='center')
    xu.rem_axLine(['top', 'right'])

    local_path = diag_config['plot_dir']
    png_name = ('comparison_zonal_turnovertime_'
                + fig_config['obs_label'] + '.png')
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
