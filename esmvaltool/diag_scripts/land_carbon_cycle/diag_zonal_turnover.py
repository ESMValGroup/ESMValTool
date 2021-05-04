"""
evaluates the zonal distribution of turnover time.

compare the model simulations with observation from
Carvalhais et al., 2014
"""

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

import esmvaltool.diag_scripts.land_carbon_cycle.plot_utils as plut
from esmvaltool.diag_scripts.land_carbon_cycle.shared import (
    _get_obs_data_zonal,
    _load_variable,
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
        'ax_fs': 7.1,
        'valrange_x': (2, 1000),
        'valrange_y': (-70, 90),
        'bandsize': 2.5,
        'gpp_threshold': 0.01
    }
    fig_config.update(diag_config.get('fig_config'))
    return fig_config


def _calc_zonal_tau(gpp, ctotal, fig_config):
    """
    Calculate zonal turnover time.

    Argument:
    --------
        gpp - cube of global gpp
        ctotal - cube of total carbon content
        fig_config - figure/diagnostic configurations

    Return:
    ------
        zonal_tau - zonal turnover time of carbon
    """
    gpp_zs = gpp.collapsed('longitude', iris.analysis.SUM)
    ctotal_zs = ctotal.collapsed('longitude', iris.analysis.SUM)

    # get the size of the sliding window based on the bandsize in degrees
    if fig_config['bandsize'] is not None:
        # get the interval of latitude and create array for partial correlation
        dat_lats = gpp.coord('latitude').points
        lat_int = abs(dat_lats[1] - dat_lats[0])
        window_size = np.int(max(2,
                                 np.round(fig_config['bandsize'] / lat_int)))
        gpp_z = gpp_zs.rolling_window('latitude', iris.analysis.SUM,
                                      window_size)
        ctotal_z = ctotal_zs.rolling_window('latitude', iris.analysis.SUM,
                                            window_size)
    else:
        gpp_z = gpp_zs
        ctotal_z = ctotal_zs

    zonal_tau = ctotal_z / gpp_z
    zonal_tau.convert_units('yr')

    return zonal_tau


def _plot_zonal_tau(plot_path, all_mod_dat, all_obs_dat, diag_config):
    """
    Plot the zonal distribution of turnover time.

    Argument:
    --------
        diag_config - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot
    """
    fig_config = _get_fig_config(diag_config)
    models = list(all_mod_dat.keys())
    models = sorted(models, key=str.casefold)
    for _mm in ['MultiModelMedian', 'MultiModelMean']:
        if _mm in models:
            models.append(models.pop(models.index(_mm)))
    nmodels = len(models)

    lats_obs = all_obs_dat['latitude']
    obs_var = diag_config.get('obs_variable')[0]
    tau_obs = all_obs_dat[obs_var]
    tau_obs_5 = all_obs_dat[obs_var + '_5']
    tau_obs_95 = all_obs_dat[obs_var + '_95']

    plt.figure(figsize=(3, 5))

    sp0 = plt.subplot(1, 1, 1)
    sp0.plot(tau_obs.data,
             lats_obs.points,
             color='k',
             lw=1.5,
             label=diag_config['obs_info']['source_label'])
    sp0.fill_betweenx(lats_obs.points,
                      tau_obs_5.data,
                      tau_obs_95.data,
                      facecolor='grey',
                      alpha=0.40)

    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_mod_tau = all_mod_dat[row_mod]
        if row_mod in ['MultiModelMedian', 'MultiModelMean']:
            sp0.plot(dat_mod_tau.data,
                     dat_mod_tau.coord('latitude').points,
                     lw=1.5,
                     color='blue',
                     label=row_mod)
        else:
            sp0.plot(dat_mod_tau.data,
                     dat_mod_tau.coord('latitude').points,
                     lw=0.5,
                     label=row_mod)

    leg = plut.draw_line_legend(ax_fs=fig_config['ax_fs'])

    plt.gca().set_xscale('log')
    plt.xlim(fig_config['valrange_x'][0], fig_config['valrange_x'][1])
    plt.ylim(fig_config['valrange_y'][0], fig_config['valrange_y'][1])
    plt.axhline(y=0, lw=0.48, color='grey')
    plt.xlabel(f'{tau_obs.long_name} ({tau_obs.units})',
               fontsize=fig_config['ax_fs'],
               ma='center')
    plt.ylabel(f'{lats_obs.long_name} ({lats_obs.units})',
               fontsize=fig_config['ax_fs'],
               ma='center')
    plut.rem_ax_line(['top', 'right'])

    plut.save_figure(plot_path, _extr_art=[leg])
    plt.close()


def main(diag_config):
    """
    Diagnostic function to compare the zonal turnover time.

    Argument:
    --------
        diag_config - nested dictionary of metadata
    """
    model_data_dict = group_metadata(diag_config['input_data'].values(),
                                     'dataset')

    fig_config = _get_fig_config(diag_config)
    zonal_tau_mod = {}
    for model_name, model_dataset in model_data_dict.items():
        zonal_tau_mod[model_name] = {}
        ctotal = _load_variable(model_dataset, 'ctotal')
        gpp = _load_variable(model_dataset, 'gpp')
        zonal_tau_mod[model_name] = _calc_zonal_tau(gpp, ctotal, fig_config)

    zonal_tau_obs = _get_obs_data_zonal(diag_config)

    obs_var = diag_config.get('obs_variable')[0]
    tau_obs = zonal_tau_obs[obs_var]
    base_name = '{title}_{source_label}_{grid_label}z'.format(
        title=tau_obs.long_name,
        source_label=diag_config['obs_info']['source_label'],
        grid_label=diag_config['obs_info']['grid_label'])

    provenance_record = _get_provenance_record(
        "Comparison of latitudinal (zonal) variations of observation-based and"
        " modelled ecosystem carbon turnover time. The zonal turnover time is"
        " calculated as the ratio of zonal `ctotal` and `gpp`. Reproduces "
        " figure 2a and 2b in Carvalhais et al. (2014).", ['mean', 'perc'],
        ['zonal'], _get_ancestor_files(diag_config, obs_var))

    if diag_config['write_netcdf']:
        model_cubes = [
            c for c in zonal_tau_mod.values() if isinstance(c, iris.cube.Cube)
        ]
        obs_cubes = [
            c for c in zonal_tau_obs.values() if isinstance(c, iris.cube.Cube)
        ]
        netcdf_path = get_diagnostic_filename(base_name, diag_config)
        save_cubes = iris.cube.CubeList(model_cubes + obs_cubes)
        iris.save(save_cubes, netcdf_path)
        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(netcdf_path, provenance_record)

    if diag_config['write_plots']:
        plot_path = get_plot_filename(base_name, diag_config)
        _plot_zonal_tau(plot_path, zonal_tau_mod, zonal_tau_obs, diag_config)
        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
