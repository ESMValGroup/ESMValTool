"""
A diagnostic function to compare the global distributions of
ecosystem carbon turnover time
"""

# operating system manipulations (e.g. path constructions)
import os
import os.path

# to manipulate iris cubes
import iris
# plotting functions
import matplotlib.pyplot as plt
import matplotlib as mpl

# map library
import cartopy.crs as ccrs

import numpy as np
import scipy.stats as stats

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename
from esmvaltool.diag_scripts.shared import group_metadata
from esmvaltool.diag_scripts.shared import run_diagnostic

# user-defined functions
import extraUtils as xu
from shared import _apply_gpp_threshold, _load_variable

# set the properties of the lines used for hatching
mpl.rcParams['hatch.color'] = 'yellow'
mpl.rcParams['hatch.linewidth'] = 0.7


# Figure settings and colorbar info
def _get_diagonal_colorbar_info():
    '''
    get the dictionary of colormap and colorbar information needed for plotting
    the maps along the diagonal, i.e., the maps of turnover time
    '''

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
    clist_ = xu.get_colomap(cb_name,
                            cb_info_diagonal['tickBounds'],
                            lowp=0.,
                            hip=1)
    cb_info_diagonal['colMap'] = mpl.colors.ListedColormap(clist_)
    return cb_info_diagonal


def _get_fig_config(diag_config):
    '''
    get the default settings of the figure, and replace default with
    runtime settings from recipe

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        a dot dictionary of settings
    '''

    fig_config = {}

    # generic settings
    fig_config['ax_fs'] = 7.1
    fig_config['fill_value'] = np.nan

    # settings of the figure and maps
    nmodels = len(
        list(
            group_metadata(diag_config['input_data'].values(),
                           'dataset').keys())) + 1
    fig_config['x0'] = 0.02
    fig_config['y0'] = 1.0
    fig_config['wp'] = 1. / nmodels
    fig_config['hp'] = fig_config['wp']
    fig_config['xsp'] = 0.0
    fig_config['ysp'] = -0.03
    fig_config['aspect_map'] = 0.5

    # settings for the location of scatterplots
    fig_config['xsp_sca'] = fig_config['wp'] / 3 * (fig_config['aspect_map'])
    fig_config['ysp_sca'] = fig_config['hp'] / 3 * (fig_config['aspect_map'])

    # colorbar specific settings
    fig_config['hcolo'] = 0.0123
    fig_config['wcolo'] = 0.25
    fig_config['cb_off_y'] = 0.06158
    fig_config['x_colo_d'] = 0.02
    fig_config['x_colo_r'] = 0.76
    fig_config['y_colo_single'] = 0.1086

    # the correlation method for metric given in the title of the scatterplot
    fig_config['correlation_method'] = 'spearman'
    fig_config['tx_y_corr'] = 1.075

    # define the range of data and masks
    fig_config['valrange_sc'] = (2, 256)
    fig_config['obs_label'] = 'Carvalhais2014'
    fig_config['obs_global'] = 23
    fig_config['gpp_threshold'] = 10  # gC m-2 yr -1

    # name of the variable and unit
    fig_config['varName'] = '$\\tau$'
    fig_config['varUnit'] = 'yr'

    # replace default values with those provided in recipe
    fig_config_list = list(fig_config.keys())
    for _fc in fig_config_list:
        if diag_config.get(_fc) is not None:
            fig_config[_fc] = diag_config.get(_fc)

    return fig_config


def _get_ratio_colorbar_info():
    '''
    get the dictionary of colormap and colorbar information needed for plotting
    the maps of ratios above the diagonal
    '''
    cb_info_ratio = {}
    border = 0.9
    ncolo = 128
    num = int(ncolo // 4)
    # get the colormap
    cb_info_ratio['tickBounds'] = np.concatenate(
        (np.geomspace(0.2, 0.25, num=num), np.geomspace(0.25, 0.33, num=num),
         np.geomspace(0.33, 0.5, num=num), np.geomspace(0.5, border, num=num),
         np.linspace(border, 1 / border, num=int(ncolo / 4)),
         np.geomspace(1 / border, 2, num=num), np.geomspace(2, 3, num=num),
         np.geomspace(3, 4, num=num), np.geomspace(4, 5, num=num)))
    colors1 = plt.cm.Blues(np.linspace(0.15, 0.998, ncolo))[::-1]
    colorsgr = np.tile(np.array([0.8, 0.8, 0.8, 1]),
                       int(ncolo / 4)).reshape(int(ncolo / 4), -1)
    colors2 = plt.cm.Reds(np.linspace(0.15, 0.998, ncolo))  # [::-1]

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


# data and calculation functions


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


def _get_agreement_mask(mmdat, dat_5, dat_95, nmodel_reject=2):
    '''
    get the mask of regiosn where fewer than one quarter of the model
    simulations are outside the range of observational uncertainty
    '''
    _maskf = np.zeros_like(mmdat)
    _maskf[(mmdat < dat_95) & (mmdat > dat_5)] = 1
    num_count = _maskf.sum(0)
    agreement_mask = np.zeros_like(num_count)
    agreement_mask[num_count < nmodel_reject] = 1
    wnan = np.isnan(dat_5)
    agreement_mask[wnan] = np.nan
    return agreement_mask


def _get_hex_data(dat_1, dat_2):
    '''
    get the data to be plotted as density plots, which requires that both the
    arrays have the same mask with regards to valid data points
    '''
    dat_1, dat_2 = _apply_common_mask(dat_1, dat_2)
    dat_1mc = np.ma.masked_equal(dat_1, np.nan).compressed()
    dat_2mc = np.ma.masked_equal(dat_2, np.nan).compressed()
    return dat_1mc, dat_2mc


def _get_obs_data(diag_config):
    '''
    Get and handle the observations of turnover time from Carvalhais 2014.

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        dictionary with observation data with different variables as keys
    '''
    if not diag_config.get('obs_files'):
        raise ValueError('The observation files needs to be specified in the '
                         'recipe (see recipe description for details)')
    else:
        input_files = [
            os.path.join(diag_config['auxiliary_data_dir'], obs_file)
            for obs_file in diag_config.get('obs_files')
        ]
    all_data = {}
    fig_config = _get_fig_config(diag_config)
    var_list = diag_config.get('obs_variables')
    for v_ind in range(len(var_list)):
        var_obs = var_list[v_ind]
        all_data[var_obs] = {}
        all_data['coords'] = {}
        variable_constraint = iris.Constraint(
            cube_func=(lambda c: c.var_name == var_obs))
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data[var_obs]['grid'] = cube.data
        all_data[var_obs]['global'] = fig_config['obs_global']
    for coord in cube.coords():
        all_data['coords'][coord.name()] = coord.points
    return all_data


def _get_turnover_data(diag_config):
    '''
    A function to calculate the modelled ecosystem carbon turnover time from
    total carbon stock and gpp.

    Arguments:
        diag_config - nested dictionary of metadata

    Returns:
        dictionaries for model simulations and observation data
    '''
    fig_config = _get_fig_config(diag_config)
    my_files_dict = group_metadata(diag_config['input_data'].values(),
                                   'dataset')
    _all_mod_dat = {}
    for key, value in my_files_dict.items():
        _all_mod_dat[key] = {}

        # load the data
        ctotal = _load_variable(value, 'ctotal')
        gpp = _load_variable(value, 'gpp')

        # calculate turnover and convert seconds to yr
        ctau = (ctotal / gpp) / (86400 * 365)

        # set the attributes and save the cube
        ctau.var_name = 'tau_ctotal'
        ctau.standard_name = None
        ctau.fill_value = fig_config['fill_value']
        ctau.long_name = 'ecosystem_carbon_turnover_time'
        ctau.units = 'yr'
        ofilename = get_diagnostic_filename(key + '_tau_ctotal', diag_config)
        iris.save(ctau, ofilename)

        # apply the GPP threshold and set the data in dictionary
        _gpp = _apply_gpp_threshold(gpp.data, fig_config)
        _ctotal = ctotal.data
        _ctau = (_ctotal / _gpp) / (86400 * 365)
        _all_mod_dat[key]['grid'] = xu.remove_invalid(
            _ctau, fill_value=fig_config['fill_value'])
        _all_mod_dat[key]['global'] = (
            (np.nansum(
                xu.remove_invalid(_ctotal,
                                  fill_value=fig_config['fill_value']))
             / np.nansum(
                 xu.remove_invalid(_gpp,
                                   fill_value=fig_config['fill_value'])))
            / (86400 * 365))

        # plot the map of the turnover time from the model
        _plot_single_map(_ctau, _all_mod_dat[key]['global'], key, diag_config)
    # get the data from the observation
    _all_obs_dat = _get_obs_data(diag_config)

    return _all_mod_dat, _all_obs_dat


# Plotting functions


def _fix_map(axis_obj):
    '''
    fixes the map boundaries, coast lines, and removes the outline box/circle
    '''
    axis_obj.set_global()
    axis_obj.coastlines(linewidth=0.4, color='grey')
    plt.gca().outline_patch.set_visible(False)
    return axis_obj


def _get_data_to_plot(_data):
    '''
    get the data to be plotted on the map corrected for the rotations of
    latitude and longitude
    '''
    xroll = _data.shape[1] / 2
    _data = np.roll(np.flipud(_data), int(xroll), axis=1)
    return _data


def _plot_matrix_map(all_mod_dat, all_obs_dat, diag_config):
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
    models.insert(0, 'obs')
    all_mod_dat['obs'] = all_obs_dat['tau_ctotal']
    nmodels = len(models)

    # define the data and information for plotting ratios
    cb_info_ratio = _get_ratio_colorbar_info()

    # get the colormap for diagonal maps
    cb_info_diagonal = _get_diagonal_colorbar_info()

    plt.figure(figsize=(9, 6))
    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_row = all_mod_dat[row_mod]['grid']
        for col_m in range(nmodels):
            col_mod = models[col_m]
            dat_col = all_mod_dat[col_mod]['grid']
            print('---' + row_mod + ' vs ' + col_mod + '---')

            # plot the maps along the diagonal
            if row_m == col_m:
                _ax = plt.axes(
                    [fig_config['x0']
                     + row_m * fig_config['wp'] + row_m * fig_config['xsp'],
                     fig_config['y0']
                     - (col_m * fig_config['hp'] + col_m * fig_config['ysp']),
                     fig_config['wp'], fig_config['hp']],
                    projection=ccrs.Robinson(central_longitude=0),
                    frameon=False)
                plot_dat = dat_row
                plt.imshow(_get_data_to_plot(plot_dat),
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
                _ax = plt.axes(
                    [(fig_config['x0'] + row_m * fig_config['wp']
                      + row_m * fig_config['xsp'] + fig_config['xsp_sca']),
                     (fig_config['y0']
                      - (col_m * fig_config['hp'] + col_m * fig_config['ysp'])
                      + fig_config['ysp_sca']),
                     fig_config['wp'] * fig_config['aspect_map'],
                     fig_config['hp'] * fig_config['aspect_map']])
                xdat, ydat = dat_col, dat_row
                dat1h, dat2h = _get_hex_data(xdat, ydat)
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
                ymin, ymax = plt.ylim()
                xmin, xmax = plt.xlim()
                plt.plot((xmin, xmax), (ymin, ymax), 'k', lw=0.1)
                if fig_config['correlation_method'] == 'pearson':
                    corr = stats.pearsonr(dat1h, dat2h)[0]
                else:
                    corr = stats.spearmanr(dat1h, dat2h)[0]
                title_str = "$R^2$=" + str(round(corr**2, 2))
                plt.title(title_str,
                          fontsize=fig_config['ax_fs'] * 0.953,
                          ma='left',
                          y=fig_config['tx_y_corr'],
                          va="top")
                if row_m != 0 and col_m != nmodels - 1:
                    xu.ax_clr(axfs=fig_config['ax_fs'])
                    xu.rotate_labels(which_ax='x',
                                     axfs=fig_config['ax_fs'], rot=90)
                elif row_m == 0 and col_m != nmodels - 1:
                    xu.ax_clrX(axfs=fig_config['ax_fs'])
                    xu.rotate_labels(which_ax='x',
                                     axfs=fig_config['ax_fs'], rot=90)
                elif col_m == nmodels - 1 and row_m != 0:
                    xu.ax_clrY(axfs=fig_config['ax_fs'])
                    xu.rotate_labels(which_ax='x',
                                     axfs=fig_config['ax_fs'], rot=90)
                if row_m == 0 and col_m == nmodels - 1:
                    xu.ax_orig(axfs=fig_config['ax_fs'])
                    xu.rotate_labels(which_ax='x',
                                     axfs=fig_config['ax_fs'], rot=90)
                    plt.ylabel('$model_{column}$',
                               fontsize=fig_config['ax_fs'])
                    plt.xlabel('$model_{row}$', fontsize=fig_config['ax_fs'])

            # plot the maps of ratio of models and observation above the
            # diagonal
            if row_m > col_m:
                _ax = plt.axes(
                    [fig_config['x0']
                     + row_m * fig_config['wp'] + row_m * fig_config['xsp'],
                     fig_config['y0']
                     - (col_m * fig_config['hp'] + col_m * fig_config['ysp']),
                     fig_config['wp'], fig_config['hp']],
                    projection=ccrs.Robinson(central_longitude=0),
                    frameon=False)
                plot_dat = xu.remove_invalid(
                    dat_row / dat_col,
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
            if col_m == 0:
                if row_mod == 'obs':
                    _title_sp = fig_config['obs_label']
                else:
                    _title_sp = row_mod
                plt.title(str(row_m + 1) + ': ' + _title_sp,
                          fontsize=0.809 * fig_config['ax_fs'])
            if row_m == nmodels - 1:
                if col_mod == 'obs':
                    _title_sp = fig_config['obs_label']
                else:
                    _title_sp = col_mod
                _title_sp = str(col_m + 1)
                t_x = _ax.text(1.1,
                               0.5,
                               _title_sp,
                               fontsize=0.809 * fig_config['ax_fs'],
                               va='center',
                               ha='center',
                               transform=_ax.transAxes)

    # plot the colorbar for maps along the diagonal
    y_colo = fig_config['y0'] + fig_config['hp'] + fig_config['cb_off_y']
    _axcol_dia = [fig_config['x_colo_d'], y_colo,
                  fig_config['wcolo'], fig_config['hcolo']]
    cb_tit_d = 'ecosystem_carbon_turnover_time (yr)'
    cb = xu.mk_colo_tau(_axcol_dia,
                        cb_info_diagonal['tickBounds'],
                        cb_info_diagonal['colMap'],
                        tick_locs=cb_info_diagonal['ticksLoc'],
                        cbfs=0.86 * fig_config['ax_fs'],
                        cbtitle=cb_tit_d,
                        cbrt=90)

    # plot the colorbar for maps above the diagonal
    y_colo = fig_config['y0'] + fig_config['hp'] + fig_config['cb_off_y']
    _axcol_rat = [fig_config['x_colo_r'], y_colo,
                  fig_config['wcolo'], fig_config['hcolo']]
    cb = xu.mk_colo_cont(_axcol_rat,
                         cb_info_ratio['tickBounds'],
                         cb_info_ratio['colMap'],
                         cbfs=0.7 * fig_config['ax_fs'],
                         cbrt=90,
                         col_scale='log',
                         cbtitle='ratio ($model_{column}$/$model_{row}$)',
                         tick_locs=cb_info_ratio['ticksLoc'])
    cb.ax.set_xticklabels(cb_info_ratio['ticksLab'],
                          fontsize=0.86 * fig_config['ax_fs'],
                          ha='center',
                          rotation=0)

    # save and close the figure
    local_path = diag_config['plot_dir']
    png_name = ('global_comparison_matrix_models_'
                + fig_config['obs_label'] + '.png')
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x],
                dpi=450)
    plt.close()

    return 'Plotted full factorial model observation comparison matrix'


def _plot_multimodel_agreement(all_mod_dat, all_obs_dat, diag_config):
    '''
    makes the maps of bias of multimodel median turnover time with
    multimodel agreement

    Arguments:
        all_mod_dat - dictionary of all model data
        all_obs_dat - dictionary of observed data
        diag_config - nested dictionary of metadata

    Returns:
        string- on completion of plotting and saving

    '''
    # get the settings for plotting figure
    fig_config = _get_fig_config(diag_config)

    # get the observation data needed to calculate the bias and multimodel
    # agreement
    tau_obs = all_obs_dat['tau_ctotal']['grid']
    tau_obs_5 = all_obs_dat['tau_ctotal_5']['grid']
    tau_obs_95 = all_obs_dat['tau_ctotal_95']['grid']

    # set the information of the colormap used for plotting bias
    cb_info = _get_ratio_colorbar_info()

    # calculate the bias of multimodel median turnover time
    models = list(all_mod_dat.keys())
    nmodels = len(models)
    dat_tau_full = np.ones((nmodels, np.shape(tau_obs)[0],
                            np.shape(tau_obs)[1])) * fig_config['fill_value']
    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_tau = all_mod_dat[row_mod]['grid']
        dat_tau_full[row_m] = dat_tau

    mm_tau = xu.remove_invalid(np.nanmedian(dat_tau_full, axis=0),
                               fill_value=fig_config['fill_value'])
    mm_bias_tau = mm_tau / tau_obs
    mm_bias_tau = xu.remove_invalid(mm_bias_tau,
                                    fill_value=fig_config['fill_value'])
    # mm_bias_tau[np.isnan(mm_bias_tau)] = 0

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
    lats = all_obs_dat['coords']['latitude']
    lons = all_obs_dat['coords']['longitude']
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

    plt.title('ecosystem_carbon_turnover_time (yr), ' + fig_config['obs_label']
              + ',\n' + 'multimodel bias and agreement',
              fontsize=0.98 * fig_config['ax_fs'])

    # plot colorbar using extraUtils
    _axcol_rat = [0.254, fig_config['y_colo_single'], 0.6, 0.035]

    col_bar = xu.mk_colo_cont(_axcol_rat,
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
    local_path = diag_config['plot_dir']
    png_str = 'global_multimodelAgreement_turnovertime_'
    png_name = png_str + fig_config['obs_label'] + '.png'

    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x],
                dpi=450)
    plt.close()

    return 'Plotted multimodel bias and agreement of turnover time'


def _plot_single_map(_dat, _datglobal, _name, diag_config):
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
    plt.imshow(_get_data_to_plot(_dat),
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
        xu.remove_invalid(_dat, fill_value=fig_config['fill_value']))
    title_str = fig_config['varName'] + ": global = " + str(round(
        _datglobal, 2)) + ", median = " + str(round(_dat_median, 2))

    plt.title('ecosystem_carbon_turnover_time (yr), ' + _name + ',\n'
              + title_str,
              fontsize=0.98 * fig_config['ax_fs'])

    # draw the colorbar
    _axcol_dia = [0.254, fig_config['y_colo_single'], 0.6, 0.035]
    xu.mk_colo_tau(_axcol_dia,
                   cb_info['tickBounds'],
                   cb_info['colMap'],
                   tick_locs=cb_info['ticksLoc'],
                   cbfs=0.86 * fig_config['ax_fs'],
                   cbtitle='',
                   cbrt=90)

    # save the figure
    local_path = diag_config['plot_dir']
    png_name = 'global_turnovertime_' + _name + '.png'
    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x],
                dpi=450)
    plt.close()

    return 'Plotted map of turnover time for the given model/observation'


if __name__ == '__main__':
    with run_diagnostic() as config:
        mod_dat_all, Obs_dat_all = _get_turnover_data(config)
        fig_config = _get_fig_config(config)
        _plot_multimodel_agreement(mod_dat_all, Obs_dat_all, config)
        _plot_matrix_map(mod_dat_all, Obs_dat_all, config)
        _plot_single_map(Obs_dat_all['tau_ctotal']['grid'],
                         Obs_dat_all['tau_ctotal']['global'],
                         fig_config['obs_label'],
                         config)
