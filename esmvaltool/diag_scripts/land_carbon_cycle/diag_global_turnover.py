"""Function to compare global distributions of turnover time."""

import os.path
import iris
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
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
    _load_variable,
    _remove_invalid,
    _var_name_constraint,
)
from esmvaltool.diag_scripts.land_carbon_cycle.provenance import (
    _get_ancestor_files,
    _get_provenance_record,
)

# set the properties of the lines used for hatching
mpl.rcParams['hatch.color'] = 'yellow'
mpl.rcParams['hatch.linewidth'] = 0.7


# Figure settings and colorbar info
def _get_diagonal_colorbar_info():
    """
    Get dictionary of colormap and colorbar information for diagonal maps.

    needed for plotting the maps along the diagonal, i.e., the maps of turnover
    time
    """
    cb_info_diagonal = {}
    cb_name = 'plasma_r'
    cb_info_diagonal['tickBounds'] = np.concatenate(
        ([1], np.linspace(8, 16, num=10)[:-1], np.linspace(16, 32,
                                                           num=10)[:-1],
         np.linspace(32, 64, num=10)[:-1], np.linspace(64, 128, num=10)[:-1],
         np.linspace(128, 256,
                     num=10)[:-1], np.linspace(256, 1000, num=2,
                                               endpoint=True)))
    cb_info_diagonal['ticksLoc'] = np.array([1, 8, 16, 32, 64, 128, 256])
    clist_ = plut.get_colomap(cb_name,
                              cb_info_diagonal['tickBounds'],
                              lowp=0.,
                              hip=1)
    cb_info_diagonal['colMap'] = mpl.colors.ListedColormap(clist_)
    return cb_info_diagonal


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
    nmodels = len(group_metadata(diag_config['input_data'].values(),
                                 'dataset')) + 1
    w_pl = 1. / nmodels
    h_pl = w_pl
    aspect_map = 0.5

    fig_config = {
        # generic settings
        'ax_fs': 7.1,
        'fill_value': np.nan,
        # settings of the figure and maps
        'x0': 0.02,
        'y0': 1.0,
        'wp': w_pl,
        'hp': h_pl,
        'xsp': 0.0,
        'ysp': -0.03,
        'aspect_map': aspect_map,
        # settings for the location of scatterplots
        'xsp_sca': w_pl / 3 * aspect_map,
        'ysp_sca': h_pl / 3 * aspect_map,
        # colorbar specific settings
        'hcolo': 0.0123,
        'wcolo': 0.25,
        'cb_off_y': 0.06158,
        'x_colo_d': 0.02,
        'x_colo_r': 0.76,
        'y_colo_single': 0.1086,
        # the correlation method for metric
        # given in the title of the scatterplot
        'correlation_method': 'spearman',
        'tx_y_corr': 1.075,
        # define the range of data and masks
        'valrange_sc': (2, 256),
        'obs_global': 23,
        'gpp_threshold': 0.01
    }
    # replace default values with those provided in recipe
    fig_config.update(diag_config.get('fig_config'))
    return fig_config


def _get_ratio_colorbar_info():
    """
    Get dictionary of colormap and colorbar information for off-diagonal maps.

    The maps of ratios above the diagonal.
    """
    cb_info_ratio = {}
    border = 0.9
    ncolo = 128
    num_gr = int(ncolo // 4)
    num_col = num_gr - 4
    # get the colormap
    cb_info_ratio['tickBounds'] = np.concatenate(
        (np.geomspace(0.2, 0.25,
                      num=num_col), np.geomspace(0.25, 0.33, num=num_col),
         np.geomspace(0.33, 0.5,
                      num=num_col), np.geomspace(0.5, border, num=num_col),
         np.linspace(border, 1 / border,
                     num=num_gr), np.geomspace(1 / border, 2, num=num_col),
         np.geomspace(2, 3, num=num_col), np.geomspace(3, 4, num=num_col),
         np.geomspace(4, 5, num=num_col)))
    colors1 = plt.cm.Blues(np.linspace(0.15, 0.998, (num_col) * 4))[::-1]
    colorsgr = np.tile(np.array([0.8, 0.8, 0.8, 1]),
                       num_gr).reshape(num_gr, -1)
    colors2 = plt.cm.Reds(np.linspace(0.15, 0.998, (num_col) * 4))

    # combine them and build a new colormap
    colors1g = np.vstack((colors1, colorsgr))
    colors = np.vstack((colors1g, colors2))
    cb_info_ratio['colMap'] = mpl.colors.LinearSegmentedColormap.from_list(
        'my_colormap', colors)
    cb_info_ratio['ticksLoc'] = [0.2, 0.25, 0.33, 0.5, 0.9, 1.1, 2, 3, 4, 5]
    cb_info_ratio['ticksLab'] = [
        '$\\dfrac{1}{5}$', '$\\dfrac{1}{4}$', '$\\dfrac{1}{3}$',
        '$\\dfrac{1}{2}$', '$\\dfrac{1}{1.1}$', '$1.1$', '$2$', '$3$', '$4$',
        '$5$'
    ]
    return cb_info_ratio


def _get_agreement_mask(mmdat, dat_5, dat_95, nmodel_reject=2):
    """
    Get mask of multimodel agreement.

    Finds regions where fewer than one quarter of the model
    simulations are outside the range of observational uncertainty.
    """
    _maskf = np.zeros_like(mmdat)
    _maskf[(mmdat < dat_95) & (mmdat > dat_5)] = 1
    num_count = _maskf.sum(0)
    agreement_mask = np.zeros_like(num_count)
    agreement_mask[num_count < nmodel_reject] = 1
    wnan = np.ma.masked_invalid(dat_5).mask
    agreement_mask[wnan] = 0.
    return agreement_mask


def _get_hex_data(dat_1, dat_2, fig_config):
    """
    Get data for density plots.

    Requires that both the arrays have the same mask with regards to valid data
    points
    """
    dat_1[(dat_1 < fig_config['valrange_sc'][0] * 0.5)] = np.nan
    dat_1[(dat_1 > fig_config['valrange_sc'][1] * 1.5)] = np.nan
    dat_2[(dat_2 < fig_config['valrange_sc'][0] * 0.5)] = np.nan
    dat_2[(dat_2 > fig_config['valrange_sc'][1] * 1.5)] = np.nan
    dat_1, dat_2 = _apply_common_mask(dat_1, dat_2)
    dat_1mc = np.ma.masked_equal(dat_1, np.nan).compressed()
    dat_2mc = np.ma.masked_equal(dat_2, np.nan).compressed()
    return dat_1mc, dat_2mc


def _get_obs_data(diag_config):
    """
    Get and handle the observations of turnover time from Carvalhais 2014.

    Argument:
    --------
        diag_config - nested dictionary of metadata

    Return:
    ------
        dictionary with observation data with different variables as keys
    """
    if not diag_config.get('obs_variable'):
        raise ValueError('The observation variable needs to be specified in '
                         'the recipe (see recipe description for details)')

    obs_dir = os.path.join(diag_config['auxiliary_data_dir'],
                           diag_config['obs_info']['obs_data_subdir'])

    all_data = {}
    all_data['global'] = {}
    all_data['grid'] = {}
    fig_config = _get_fig_config(diag_config)
    var_list = diag_config.get('obs_variable')

    input_files = []
    for _var in var_list:
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var,
                                                               perc=5))
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var,
                                                               perc=95))
        obs_filename = (f'{_var}_{{frequency}}_{{source_label}}_'
                        f'{{variant_label}}_{{grid_label}}.nc'.format(
                            **diag_config['obs_info']))
        input_files = np.append(input_files,
                                os.path.join(obs_dir, obs_filename))
    nvars = len(var_list)
    for v_ind in range(nvars):
        var_obs = var_list[v_ind]
        all_data['coords'] = {}
        variable_constraint = _var_name_constraint(var_obs)
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data['grid'][var_obs] = cube
        all_data['global'][var_obs] = fig_config['obs_global']
    for coord in cube.coords():
        all_data['coords'][coord.name()] = coord.points

    all_data['input_files'] = input_files
    return all_data


def _calc_turnover(ctotal, gpp, _model):
    """
    Calculate the turnover time from ctotal and gpp.

    Argument:
    --------
        ctotal- iris cube of total carbon stock
        gpp - iris cube of gross primary productivity

    Return:
    ------
        tau_ctotal - iris cube of turnover time in years
    """
    # calculate turnover and convert units to yr
    tau_ctotal = (ctotal / gpp)
    tau_ctotal.convert_units('yr')

    # set the attributes
    tau_ctotal.var_name = 'tau_ctotal'
    tau_ctotal.standard_name = None
    tau_ctotal.long_name = 'ecosystem_carbon_turnover_time'
    tau_ctotal.units = 'yr'

    return tau_ctotal


def _fix_map(axis_obj):
    """
    Beautify map object.

    Clean boundaries, coast lines, and removes the outline box/circle.
    """
    axis_obj.set_global()
    axis_obj.coastlines(linewidth=0.4, color='grey')
    plt.gca().outline_patch.set_visible(False)
    return axis_obj


def _get_data_to_plot(_data):
    """
    Get data to plot on map.

    Correct for the rotations of latitude and longitude.
    """
    xroll = _data.shape[1] / 2
    _data = np.roll(np.flipud(_data), int(xroll), axis=1)
    return _data


def _get_matrix_map_axes(_row_m, _col_m, _fig_config):
    """
    Get the axes object for matrix maps.

    Argument:
    --------
        _row_m - row location in the matrix
        _col_m - column location in the matrix
        _fig_config - figure settings

    Return:
    ------
        _ax - an axes object
    """
    if _row_m == _col_m:
        _ax = plt.axes([
            _fig_config['x0'] + _row_m * _fig_config['wp'] +
            _row_m * _fig_config['xsp'], _fig_config['y0'] -
            (_col_m * _fig_config['hp'] + _col_m * _fig_config['ysp']),
            _fig_config['wp'], _fig_config['hp']
        ], projection=ccrs.Robinson(central_longitude=0), frameon=False)
    if _row_m < _col_m:
        _ax = plt.axes([
            _fig_config['x0'] + _row_m * _fig_config['wp'] +
            _row_m * _fig_config['xsp'] + _fig_config['xsp_sca'],
            _fig_config['y0'] -
            (_col_m * _fig_config['hp'] + _col_m * _fig_config['ysp']) +
            _fig_config['ysp_sca'],
            _fig_config['wp'] * _fig_config['aspect_map'],
            _fig_config['hp'] * _fig_config['aspect_map']
        ])

    if _row_m > _col_m:
        _ax = plt.axes([
            _fig_config['x0'] + _row_m * _fig_config['wp'] +
            _row_m * _fig_config['xsp'], _fig_config['y0'] -
            (_col_m * _fig_config['hp'] + _col_m * _fig_config['ysp']),
            _fig_config['wp'], _fig_config['hp']
        ], projection=ccrs.Robinson(central_longitude=0), frameon=False)
    return _ax


def _fix_matrix_axes(row_m, col_m, models, nmodels, diag_config, fig_config):
    """Fix the axes lines and titles in matrix maps."""
    row_mod = models[row_m]
    col_mod = models[col_m]
    if row_m != 0 and col_m != nmodels - 1:
        plut.ax_clr()
        plut.rotate_labels(which_ax='x', axfs=fig_config['ax_fs'], rot=90)
    elif row_m == 0 and col_m != nmodels - 1:
        plut.ax_clr_x(axfs=fig_config['ax_fs'])
        plut.rotate_labels(which_ax='x', axfs=fig_config['ax_fs'], rot=90)
    elif col_m == nmodels - 1 and row_m != 0:
        plut.ax_clr_y(axfs=fig_config['ax_fs'])
        plut.rotate_labels(which_ax='x', axfs=fig_config['ax_fs'], rot=90)
    if row_m == 0 and col_m == nmodels - 1:
        plut.ax_orig(axfs=fig_config['ax_fs'])
        plut.rotate_labels(which_ax='x', axfs=fig_config['ax_fs'], rot=90)
        plt.ylabel('$model_{column}$', fontsize=fig_config['ax_fs'])
        plt.xlabel('$model_{row}$', fontsize=fig_config['ax_fs'])
    if col_m == 0:
        if row_mod == 'obs':
            _title_sp = diag_config['obs_info']['source_label']
        else:
            _title_sp = row_mod
        plt.title(str(row_m + 1) + ': ' + _title_sp,
                  fontsize=0.809 * fig_config['ax_fs'])
    if row_m == nmodels - 1:
        if col_mod == 'obs':
            _title_sp = diag_config['obs_info']['source_label']
        else:
            _title_sp = col_mod
        _title_sp = str(col_m + 1)
        t_x = plt.gca().text(1.1,
                             0.5,
                             _title_sp,
                             fontsize=0.809 * fig_config['ax_fs'],
                             va='center',
                             ha='center',
                             transform=plt.gca().transAxes)
    else:
        t_x = ''

    return t_x


def _draw_121_line():
    """Draw 1:1 line on the current axis."""
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    plt.plot((xmin, xmax), (ymin, ymax), 'k', lw=0.1)


def _plot_matrix_map(plot_path_matrix, global_tau_mod, global_tau_obs,
                     diag_config):
    """
    Plot the matrix of maps model-observation full factorial comparison.

    Argument:
    --------
        diag_config - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot
    """
    fig_config = _get_fig_config(diag_config)
    models = list(global_tau_mod['grid'].keys())
    models = sorted(models, key=str.casefold)
    multimodel_stats = 'MultiModelMedian MultiModelMean'.split()
    for _mm in multimodel_stats:
        if _mm in models:
            models.append(models.pop(models.index(_mm)))
    models.insert(0, 'obs')

    global_tau_mod['grid']['obs'] = global_tau_obs['grid']['tau_ctotal']
    global_tau_mod['global']['obs'] = global_tau_obs['global']['tau_ctotal']
    nmodels = len(models)

    # define the data and information for plotting ratios
    cb_info_ratio = _get_ratio_colorbar_info()

    # get the colormap for diagonal maps
    cb_info_diagonal = _get_diagonal_colorbar_info()

    plt.figure(figsize=(9, 6))
    for row_m in range(nmodels):
        dat_row = global_tau_mod['grid'][models[row_m]].data
        for col_m in range(nmodels):
            dat_col = global_tau_mod['grid'][models[col_m]].data
            _ax = _get_matrix_map_axes(row_m, col_m, fig_config)
            # plot the maps along the diagonal
            if row_m == col_m:
                plt.imshow(_get_data_to_plot(dat_row),
                           norm=mpl.colors.BoundaryNorm(
                               cb_info_diagonal['tickBounds'],
                               len(cb_info_diagonal['tickBounds'])),
                           cmap=cb_info_diagonal['colMap'],
                           origin='upper',
                           vmin=cb_info_diagonal['tickBounds'][0],
                           vmax=cb_info_diagonal['tickBounds'][-1],
                           transform=ccrs.PlateCarree())
                _fix_map(_ax)

            # plot the scatterplot/density plot below the diagonal
            if row_m < col_m:
                dat1h, dat2h = _get_hex_data(dat_col, dat_row, fig_config)
                _ax.hexbin(dat1h,
                           dat2h,
                           bins='log',
                           mincnt=3,
                           gridsize=40,
                           cmap='viridis_r',
                           linewidths=0)
                plt.ylim(fig_config['valrange_sc'][0],
                         fig_config['valrange_sc'][1] * 1.05)
                plt.xlim(fig_config['valrange_sc'][0],
                         fig_config['valrange_sc'][1] * 1.05)
                _draw_121_line()
                if fig_config['correlation_method'] == 'pearson':
                    corr = (stats.pearsonr(dat1h, dat2h)[0])**2
                else:
                    corr = (stats.spearmanr(dat1h, dat2h)[0])**2
                plt.title('$R^2$={corr:.2f}'.format(corr=corr),
                          fontsize=fig_config['ax_fs'] * 0.953,
                          ma='left',
                          y=fig_config['tx_y_corr'],
                          va="top")

            # plot the maps of ratio of models and observation above the
            # diagonal
            if row_m > col_m:
                plot_dat = _remove_invalid(dat_row / dat_col,
                                           fill_value=fig_config['fill_value'])
                _ax.imshow(_get_data_to_plot(plot_dat),
                           norm=mpl.colors.BoundaryNorm(
                               cb_info_ratio['tickBounds'],
                               len(cb_info_ratio['tickBounds'])),
                           interpolation='none',
                           vmin=cb_info_ratio['tickBounds'][0],
                           vmax=cb_info_ratio['tickBounds'][-1],
                           cmap=cb_info_ratio['colMap'],
                           origin='upper',
                           transform=ccrs.PlateCarree())
                _fix_map(_ax)
            t_x = _fix_matrix_axes(row_m, col_m, models, nmodels, diag_config,
                                   fig_config)

    # plot the colorbar for maps along the diagonal
    y_colo = fig_config['y0'] + fig_config['hp'] + fig_config['cb_off_y']
    _axcol_dia = [
        fig_config['x_colo_d'], y_colo, fig_config['wcolo'],
        fig_config['hcolo']
    ]
    cb_tit_d = '{name} ({unit})'.format(
        name=global_tau_mod['grid'][models[col_m]].long_name,
        unit=global_tau_mod['grid'][models[col_m]].units)
    col_bar = plut.mk_colo_tau(_axcol_dia,
                               cb_info_diagonal['tickBounds'],
                               cb_info_diagonal['colMap'],
                               tick_locs=cb_info_diagonal['ticksLoc'],
                               cbfs=0.86 * fig_config['ax_fs'],
                               cbtitle=cb_tit_d,
                               cbrt=90)

    # plot the colorbar for maps above the diagonal
    y_colo = fig_config['y0'] + fig_config['hp'] + fig_config['cb_off_y']
    _axcol_rat = [
        fig_config['x_colo_r'], y_colo, fig_config['wcolo'],
        fig_config['hcolo']
    ]
    col_bar = plut.mk_colo_cont(
        _axcol_rat,
        cb_info_ratio['tickBounds'],
        cb_info_ratio['colMap'],
        cbfs=0.7 * fig_config['ax_fs'],
        cbrt=90,
        col_scale='log',
        cbtitle='ratio ($model_{column}$/$model_{row}$)',
        tick_locs=cb_info_ratio['ticksLoc'])
    col_bar.ax.set_xticklabels(cb_info_ratio['ticksLab'],
                               fontsize=0.86 * fig_config['ax_fs'],
                               ha='center',
                               rotation=0)

    # save and close figure
    plut.save_figure(plot_path_matrix, _extr_art=[t_x])
    plt.close()


def _plot_multimodel_agreement(plot_path_multimodel, global_tau_mod,
                               global_tau_obs, diag_config):
    """
    Plot map of multimodel bias and multimodel agreement.

    Argument:
    --------
        global_tau_mod - dictionary of all model data
        global_tau_obs - dictionary of observed data
        diag_config - nested dictionary of metadata
    """
    # get the settings for plotting figure
    fig_config = _get_fig_config(diag_config)

    # get the observation data needed to calculate the bias and multimodel
    # agreement
    obs_var = diag_config.get('obs_variable')[0]
    tau_obs = global_tau_obs['grid'][obs_var].data
    tau_obs_5 = global_tau_obs['grid'][obs_var + '_5'].data
    tau_obs_95 = global_tau_obs['grid'][obs_var + '_95'].data

    # set the information of the colormap used for plotting bias
    cb_info = _get_ratio_colorbar_info()

    # calculate the bias of multimodel median turnover time
    models = list(global_tau_mod['grid'].keys())

    # remove multimodel estimates from the list of models
    multimodel_stats = 'MultiModelMedian MultiModelMean'.split()
    for _mm in multimodel_stats:
        if _mm in models:
            models.remove(_mm)

    nmodels = len(models)
    dat_tau_full = np.ones((nmodels, np.shape(tau_obs)[0],
                            np.shape(tau_obs)[1])) * fig_config['fill_value']
    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_tau = global_tau_mod['grid'][row_mod]
        dat_tau_full[row_m] = _remove_invalid(
            dat_tau.data, fill_value=fig_config['fill_value'])

    mm_tau = _remove_invalid(np.nanmedian(dat_tau_full, axis=0),
                             fill_value=fig_config['fill_value'])
    mm_bias_tau = mm_tau / tau_obs
    mm_bias_tau = _remove_invalid(mm_bias_tau,
                                  fill_value=fig_config['fill_value'])

    # define figure and main axis to plot the map
    plt.figure(figsize=(5, 3))
    _ax = plt.axes([0.1, 0.1, 0.9, 0.9],
                   projection=ccrs.Robinson(central_longitude=0),
                   frameon=False)

    # plot the data of multimodel bias (=bias of multimodel median turnover
    # time)
    _ax.imshow(_get_data_to_plot(mm_bias_tau),
               norm=mpl.colors.BoundaryNorm(cb_info['tickBounds'],
                                            len(cb_info['tickBounds'])),
               interpolation='none',
               vmin=cb_info['tickBounds'][0],
               vmax=cb_info['tickBounds'][-1],
               cmap=cb_info['colMap'],
               origin='upper',
               transform=ccrs.PlateCarree())
    _fix_map(_ax)

    # get the model agreement mask (less than quarter of the model within the
    # observational uncertainty)
    agreement_mask_tau = _get_agreement_mask(dat_tau_full,
                                             tau_obs_5,
                                             tau_obs_95,
                                             nmodel_reject=int(nmodels / 4))

    # plot the hatches for uncertainty/multimodel agreement
    lats = global_tau_obs['coords']['latitude']
    lons = global_tau_obs['coords']['longitude']
    latint = abs(lats[1] - lats[0])
    lonint = abs(lons[1] - lons[0])
    x_lat, y_lon = np.meshgrid(lons - lonint / 2, lats - latint / 2)

    _ax.contourf(x_lat,
                 y_lon,
                 agreement_mask_tau,
                 levels=[0, 0.5, 1],
                 alpha=0.,
                 hatches=['', '//////'],
                 linewidth=0.2,
                 transform=ccrs.PlateCarree())

    title_str = ('multimodel bias and agreement (-)\n{title}'.format(
        title=global_tau_obs['grid']['tau_ctotal'].long_name))
    plt.title(title_str, fontsize=0.98 * fig_config['ax_fs'])

    # plot colorbar using extraUtils
    _axcol_rat = [0.254, fig_config['y_colo_single'], 0.6, 0.035]

    col_bar = plut.mk_colo_cont(_axcol_rat,
                                cb_info['tickBounds'],
                                cb_info['colMap'],
                                cbfs=0.8 * fig_config['ax_fs'],
                                cbrt=90,
                                col_scale='log',
                                cbtitle='',
                                tick_locs=cb_info['ticksLoc'])
    col_bar.ax.set_xticklabels(cb_info['ticksLab'],
                               fontsize=0.9586 * fig_config['ax_fs'],
                               ha='center',
                               rotation=0)

    # save and close figure
    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)
    plut.save_figure(plot_path_multimodel, _extr_art=[t_x])
    plt.close()


def _plot_single_map(plot_path, _dat, _datglobal, _name, provenance_record,
                     diag_config):
    """
    Plot a map for a given variable.

    Argument:
    --------
        diag_config - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot
    """
    # figure configuration
    fig_config = _get_fig_config(diag_config)

    # colormap configuration
    cb_info = _get_diagonal_colorbar_info()

    # define the figure and axis
    plt.figure(figsize=(5, 3))
    _ax = plt.axes([0.1, 0.1, 0.9, 0.9],
                   projection=ccrs.Robinson(central_longitude=0),
                   frameon=False)
    # plot data over the map
    plt.imshow(_get_data_to_plot(_dat.data),
               norm=mpl.colors.BoundaryNorm(cb_info['tickBounds'],
                                            len(cb_info['tickBounds'])),
               cmap=cb_info['colMap'],
               origin='upper',
               vmin=cb_info['tickBounds'][0],
               vmax=cb_info['tickBounds'][-1],
               transform=ccrs.PlateCarree())
    _fix_map(_ax)

    # get the data and set the title of the map

    _dat_median = np.nanmedian(
        _remove_invalid(_dat.data, fill_value=fig_config['fill_value']))
    title_str = (f'{_dat.long_name} ({_dat.units}), {_name},\n'
                 f'global = {_datglobal:.2f}, median = {_dat_median:.2f}')

    plt.title(title_str, fontsize=0.98 * fig_config['ax_fs'])

    # draw the colorbar
    _axcol_dia = [0.254, fig_config['y_colo_single'], 0.6, 0.035]
    plut.mk_colo_tau(_axcol_dia,
                     cb_info['tickBounds'],
                     cb_info['colMap'],
                     tick_locs=cb_info['ticksLoc'],
                     cbfs=0.86 * fig_config['ax_fs'],
                     cbtitle='',
                     cbrt=90)

    # save and close figure
    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)
    plut.save_figure(plot_path, _extr_art=[t_x])
    plt.close()
    with ProvenanceLogger(diag_config) as provenance_logger:
        provenance_logger.log(plot_path, provenance_record)


def main(diag_config):
    """
    Evaluate global distribution of ecosystem carbon turnover time.

    Argument:
    --------
        diag_config - nested dictionary of metadata
    """
    model_data_dict = group_metadata(diag_config['input_data'].values(),
                                     'dataset')

    # get the data from the observation
    global_tau_obs = _get_obs_data(diag_config)
    base_name = ('{title}_{source_label}_'
                 '{grid_label}'.format(
                     title=global_tau_obs['grid']['tau_ctotal'].long_name,
                     source_label=diag_config['obs_info']['source_label'],
                     grid_label=diag_config['obs_info']['grid_label']))

    global_tau_mod = {}
    global_tau_mod['grid'] = {}
    global_tau_mod['global'] = {}

    provenance_record_matrix = _get_provenance_record(
        "Matrix Comparison of global distributions of turnover time of carbon",
        ['mean', 'perc'], ['global'],
        _get_ancestor_files(diag_config, 'tau_ctotal'))

    provenance_record_multimodel = _get_provenance_record(
        "Multimodel bias and agreements of global distributions of turnover"
        "time of carbon. Reproduces figure 3 in Carvalhais et al. (2014).",
        ['mean', 'perc'], ['global'],
        _get_ancestor_files(diag_config, 'tau_ctotal'))

    for model_name, model_dataset in model_data_dict.items():
        global_tau_mod[model_name] = {}

        # load the data
        ctotal = _load_variable(model_dataset, 'ctotal')
        gpp = _load_variable(model_dataset, 'gpp')
        tau_ctotal = _calc_turnover(ctotal, gpp, model_name)
        global_tau_mod['grid'][model_name] = tau_ctotal

        # apply the GPP threshold and set the data in dictionary
        gpp_global = gpp.collapsed(['latitude', 'longitude'],
                                   iris.analysis.SUM)
        ctotal_global = ctotal.collapsed(['latitude', 'longitude'],
                                         iris.analysis.SUM)
        tau_global = ctotal_global / gpp_global
        tau_global.convert_units('yr')

        global_tau_mod['global'][model_name] = np.float(tau_global
                                                        .core_data())

        if diag_config['write_plots']:
            base_name_mod = (
                'global_{title}_{source_label}_'
                '{grid_label}'.format(
                    title=global_tau_obs['grid']['tau_ctotal'].long_name,
                    source_label=model_name,
                    grid_label=diag_config['obs_info']['grid_label']))
            plot_path_mod = get_plot_filename(base_name_mod, diag_config)
            # plot_path_list.append(plot_path_mod)
            provenance_record_mod = _get_provenance_record(
                "Map of global distribution of turnover time of carbon",
                ['mean', 'perc'],
                ['global'],
                {model_name: model_dataset})
            _plot_single_map(plot_path_mod, tau_ctotal,
                             global_tau_mod['global'][model_name],
                             model_name,
                             provenance_record_mod,
                             diag_config)
        if diag_config['write_netcdf']:
            model_cubes = [
                c for c in global_tau_mod['grid'].values()
                if isinstance(c, iris.cube.Cube)
            ]
            obs_cubes = [
                c for c in global_tau_obs['grid'].values()
                if isinstance(c, iris.cube.Cube)
            ]
            netcdf_path = get_diagnostic_filename(base_name_mod, diag_config)
            save_cubes = iris.cube.CubeList(model_cubes + obs_cubes)
            iris.save(save_cubes, netcdf_path)

        else:
            netcdf_path = None

        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(netcdf_path, provenance_record_mod)

    if diag_config['write_plots']:
        # multimodel agreement
        base_name_multimodel = '{prefix}_{base_name}'.format(
            prefix='global_multimodelAgreement', base_name=base_name)
        plot_path_multimodel = get_plot_filename(base_name_multimodel,
                                                 diag_config)
        _plot_multimodel_agreement(plot_path_multimodel, global_tau_mod,
                                   global_tau_obs, config)
        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(plot_path_multimodel,
                                  provenance_record_multimodel)

        # map of observation
        base_name_obs = '{prefix}_{base_name}'.format(prefix='global',
                                                      base_name=base_name)
        plot_path_obs = get_plot_filename(base_name_obs, diag_config)
        provenance_record_obs = _get_provenance_record(
            "Map of observed global distribution of turnover time of carbon",
            ['mean', 'perc'],
            ['global'],
            global_tau_obs['input_files'].tolist())

        _plot_single_map(plot_path_obs,
                         global_tau_obs['grid']['tau_ctotal'],
                         global_tau_obs['global']['tau_ctotal'],
                         config['obs_info']['source_label'],
                         provenance_record_obs,
                         diag_config)

        # matrix of maps
        base_name_matrix = '{prefix}_{base_name}'.format(
            prefix='global_matrix_map', base_name=base_name)
        plot_path_matrix = get_plot_filename(base_name_matrix, diag_config)
        _plot_matrix_map(plot_path_matrix, global_tau_mod, global_tau_obs,
                         config)

        with ProvenanceLogger(diag_config) as provenance_logger:
            provenance_logger.log(plot_path_matrix, provenance_record_matrix)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
