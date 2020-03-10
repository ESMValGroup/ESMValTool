"""
A diagnostic function to compare the global distributions of
ecosystem carbon turnover time
"""
# place your module imports here:

# operating system manipulations (e.g. path constructions)
import sys, os, os.path

# to manipulate iris cubes
import iris
# plotting functions
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors

# map library
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

import numpy as np
# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename, group_metadata, select_metadata, run_diagnostic, extract_variables
from esmvalcore.preprocessor import area_statistics

# user-defined functions
import extraUtils as xu
mpl.rcParams['hatch.color'] = 'yellow'
mpl.rcParams['hatch.linewidth'] = 0.7

# 

class dotDict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


# Figure settings and colorbar info


def _get_dia_colorbarInfo():
    cbInfo_diagonal = {}
    cbName = 'plasma_r'
    cbInfo_diagonal['tickBounds'] = np.concatenate(
        ([1], np.linspace(8, 16, num=10)[:-1], np.linspace(16, 32,
                                                           num=10)[:-1],
         np.linspace(32, 64, num=10)[:-1], np.linspace(64, 128, num=10)[:-1],
         np.linspace(128, 256,
                     num=10)[:-1], np.linspace(256, 1000, num=2,
                                               endpoint=True)))
    cbInfo_diagonal['ticksLoc'] = np.array([1, 8, 16, 32, 64, 128, 256])
    clist_ = xu.get_colomap(cbName, cbInfo_diagonal['tickBounds'], lowp=0., hip=1)
    cbInfo_diagonal['colMap'] = mpl.colors.ListedColormap(clist_)
    return cbInfo_diagonal


def _get_fig_config(__cfg):
    """
    get the default settings of the figure, and replace default with 
    runtime settings from recipe

    Arguments:
        __cfg - nested dictionary of metadata

    Returns:
        a dot dictionary of settings
    """

    fcfg = {}

    # generic settings
    fcfg['ax_fs'] = 7.1
    fcfg['fill_val'] = np.nan

    # settings of the figure and maps
    nmodels = len(
        list(group_metadata(__cfg['input_data'].values(),
                            'dataset').keys())) + 1
    fcfg['x0'] = 0.02
    fcfg['y0'] = 1.0
    fcfg['wp'] = 1. / nmodels
    fcfg['hp'] = fcfg['wp']
    fcfg['xsp'] = 0.0
    fcfg['ysp'] = -0.03
    fcfg['aspect_map'] = 0.5

    # settings for the location of scatterplots
    fcfg['xsp_sca'] = fcfg['wp'] / 3 * (fcfg['aspect_map'])
    fcfg['ysp_sca'] = fcfg['hp'] / 3 * (fcfg['aspect_map'])

    # colorbar specific settings
    fcfg['hcolo'] = 0.0123
    fcfg['wcolo'] = 0.25
    fcfg['cb_off_y'] = 0.06158
    fcfg['x_colo_d'] = 0.02
    fcfg['x_colo_r'] = 0.76
    fcfg['y_colo_single'] = 0.1086

    # the correlation method for metric given in the title of the scatterplot
    fcfg['correlation_method'] = 'spearman'
    fcfg['tx_y_corr'] = 1.075

    # define the range of data and masks
    fcfg['valrange_sc'] = (2, 256)
    fcfg['obs_label'] = 'Carvalhais2014'
    fcfg['obs_global'] = 23
    fcfg['gpp_threshold'] = 10  #gC m-2 yr -1

    # name of the variable and unit
    fcfg['varName'] = '$\\tau$'
    fcfg['varUnit'] = 'yr'

    # replace default values with those provided in recipe
    fcfgL = list(fcfg.keys())
    for _fc in fcfgL:
        if __cfg.get(_fc) != None:
            fcfg[_fc] = __cfg.get(_fc)

    # convert dictionary to dot notation
    _figSet = dotDict(fcfg)
    return _figSet


def _get_ratio_colorbarInfo():
    cbInfo_ratio = {}
    border = 0.9
    ncolo = 128
    # get the colormap
    cbInfo_ratio['tickBounds'] = np.concatenate(
        (np.geomspace(0.2, 0.25, num=int(ncolo) / 4),
         np.geomspace(0.25, 0.33, num=int(ncolo) / 4),
         np.geomspace(0.33, 0.5, num=int(ncolo) / 4),
         np.geomspace(0.5, border, num=int(ncolo) / 4),
         np.linspace(border, 1 / border, num=int(ncolo / 4)),
         np.geomspace(1 / border, 2, num=int(ncolo) / 4),
         np.geomspace(2, 3, num=int(ncolo) / 4),
         np.geomspace(3, 4, num=int(ncolo) / 4),
         np.geomspace(4, 5, num=int(ncolo) / 4)))
    colors1 = plt.cm.Blues(np.linspace(0.15, 0.998, ncolo))[::-1]
    colorsgr = np.tile(np.array([0.8, 0.8, 0.8, 1]),
                       int(ncolo / 4)).reshape(int(ncolo / 4), -1)
    colors2 = plt.cm.Reds(np.linspace(0.15, 0.998, ncolo))  #[::-1]

    # combine them and build a new colormap
    colors1g = np.vstack((colors1, colorsgr))
    colors = np.vstack((colors1g, colors2))
    cbInfo_ratio['colMap'] = mpl.colors.LinearSegmentedColormap.from_list(
        'my_colormap', colors)
    cbInfo_ratio['ticksLoc'] = [0.2, 0.25, 0.33, 0.5, 0.9, 1.1, 2, 3, 4, 5]
    cbInfo_ratio['ticksLab'] = [
        '$\\dfrac{1}{5}$', '$\\dfrac{1}{4}$', '$\\dfrac{1}{3}$',
        '$\\dfrac{1}{2}$', '$\\dfrac{1}{1.1}$', '$1.1$', '$2$', '$3$', '$4$',
        '$5$'
    ]
    return cbInfo_ratio


# data and calculation functions


def _apply_common_mask(_dat1, _dat2):
    '''
    returns the copies of two arrays with common mask applied
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
    # converting gC m-2 yr-1 to kgC m-2 s-1
    gpp_thres = _fcfg.gpp_threshold / (86400.0 * 365 * 1000.)
    _gppDat = np.ma.masked_less(_gppDat, gpp_thres).filled(_fcfg.fill_val)
    return _gppDat


def _get_agreement_mask(_mmdat,
                        _dat_5,
                        _dat_95,
                        _nmodels,
                        fill_val=-9999.,
                        nmodel_reject=2):
    _maskf = np.zeros_like(_mmdat)
    _maskf[(_mmdat < _dat_95) & (_mmdat > _dat_5)] = 1
    nCount = _maskf.sum(0)
    agreement_mask = np.zeros_like(nCount)
    agreement_mask[nCount < nmodel_reject] = 1
    wnan = np.isnan(_dat_5)
    agreement_mask[wnan] = np.nan
    return agreement_mask


def _get_hex_data(_dat1, _dat2, valrange_sc):
    _dat1, _dat2 = _apply_common_mask(_dat1, _dat2)
    _dat1mc = np.ma.masked_equal(_dat1, np.nan).compressed()
    _dat2mc = np.ma.masked_equal(_dat2, np.nan).compressed()
    return _dat1mc, _dat2mc


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
    fcfg = _get_fig_config(cfg)
    varList = cfg.get('obs_variables')
    # varList.append('latitude')
    print(varList)
    for _var in varList:
        all_data[_var] = {}
        all_data['coords'] = {}
        variable_constraint = iris.Constraint(
            cube_func=(lambda c: c.var_name == _var))
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data[_var]['grid'] = cube.data
        all_data[_var]['global'] = fcfg.obs_global
    for coord in cube.coords():
        print(coord.name())
        all_data['coords'][coord.name()] = coord.points
    return all_data


def _get_turnover_data(cfg):
    """
    A function to calculate the modelled ecosystem carbon turnover time from
    total carbon stock and gpp.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        dictionaries for model simulations and observation data
    """
    fcfg = _get_fig_config(cfg)
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    _allModdat = {}
    for key, value in my_files_dict.items():
        _allModdat[key] = {}

        # load the data
        ctotal = _load_variable(value, 'ctotal')
        gpp = _load_variable(value, 'gpp')

        # calculate turnover and convert seconds to yr
        ctau = (ctotal / gpp) / (86400 * 365)

        # set the attributes and save the cube
        ctau.var_name = 'tau_ctotal'
        ctau.standard_name = None
        ctau.fill_value = fcfg.fill_val
        ctau.long_name = 'ecosystem_carbon_turnover_time'
        ctau.units = 'yr'
        ofilename = get_diagnostic_filename(key + '_tau_ctotal', cfg)
        iris.save(ctau, ofilename)

        # apply the GPP threshold and set the data in dictionary
        _gpp = _apply_gpp_threshold(gpp.data, fcfg)
        _ctotal = ctotal.data
        _ctau = (_ctotal / _gpp) / (86400 * 365)
        _allModdat[key]['grid'] = xu.remove_invalid(_ctau,
                                                    _fill_val=fcfg.fill_val)
        _allModdat[key]['global'] = (
            np.nansum(xu.remove_invalid(_ctotal, _fill_val=fcfg.fill_val)) /
            np.nansum(xu.remove_invalid(_gpp, _fill_val=fcfg.fill_val))) / (
                86400 * 365)

        # plot the map of the turnover time from the model
        _plot_single_map(_ctau, _allModdat[key]['global'], key, cfg)
    # get the data from the observation
    _allObsdat = _get_obs_data(cfg)

    return _allModdat, _allObsdat


def _load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube


# Plotting functions


def _fix_map(__ax, __fcfg):
    __ax.set_global()
    __ax.coastlines(linewidth=0.4, color='grey')
    plt.gca().outline_patch.set_visible(False)
    return __ax


def _get_data_to_plot(_data, _fcfg):
    xroll = _data.shape[1] / 2
    _data = np.roll(np.flipud(_data), int(xroll), axis=1)
    return _data


def _plot_matrix_map(all_mod_dat, all_obs_dat, cfg):
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
    models.insert(0, 'obs')
    all_mod_dat['obs'] = all_obs_dat['tau_ctotal']
    nmodels = len(models)

    # define the data and information for plotting ratios
    cbInfo_ratio = _get_ratio_colorbarInfo()

    # get the colormap for diagonal maps
    cbInfo_diagonal = _get_dia_colorbarInfo()

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
                _ax = plt.axes([
                    _fcfg.x0 + row_m * _fcfg.wp + row_m * _fcfg.xsp, _fcfg.y0 -
                    (col_m * _fcfg.hp + col_m * _fcfg.ysp), _fcfg.wp, _fcfg.hp
                ],
                               projection=ccrs.Robinson(central_longitude=0),
                               frameon=False)
                plot_dat = dat_row
                plt.imshow(_get_data_to_plot(plot_dat, _fcfg),
                           norm=matplotlib.colors.BoundaryNorm(
                               cbInfo_diagonal['tickBounds'], len(cbInfo_diagonal['tickBounds'])),
                           cmap=cbInfo_diagonal['colMap'],
                           origin='upper',
                           vmin=cbInfo_diagonal['tickBounds'][0],
                           vmax=cbInfo_diagonal['tickBounds'][-1],
                           transform=ccrs.PlateCarree())
                _fix_map(_ax, _fcfg)

            # plot the scatterplot/density plot below the diagonal
            if row_m < col_m:
                _ax = plt.axes([
                    _fcfg.x0 + row_m * _fcfg.wp + row_m * _fcfg.xsp +
                    _fcfg.xsp_sca, _fcfg.y0 -
                    (col_m * _fcfg.hp + col_m * _fcfg.ysp) + _fcfg.ysp_sca,
                    _fcfg.wp * _fcfg.aspect_map, _fcfg.hp * _fcfg.aspect_map
                ])
                xdat, ydat = dat_col, dat_row
                dat1h, dat2h = _get_hex_data(xdat, ydat, _fcfg.valrange_sc)
                _ax.hexbin(dat1h,
                           dat2h,
                           bins='log',
                           mincnt=3,
                           gridsize=40,
                           cmap='viridis_r',
                           linewidths=0)
                plt.ylim(_fcfg.valrange_sc[0], _fcfg.valrange_sc[1] * 1.05)
                plt.xlim(_fcfg.valrange_sc[0], _fcfg.valrange_sc[1] * 1.05)
                ymin, ymax = plt.ylim()
                xmin, xmax = plt.xlim()
                plt.plot((xmin, xmax), (ymin, ymax), 'k', lw=0.1)
                if _fcfg.correlation_method == 'pearson':
                    r, p = xu.calc_pearson_r(xdat, ydat)
                else:
                    r, p = xu.calc_spearman_r(xdat, ydat)
                title_str = "$R^2$=" + str(round(r**2, 2))
                plt.title(title_str,
                          fontsize=_fcfg.ax_fs * 0.953,
                          ma='left',
                          y=_fcfg.tx_y_corr,
                          va="top")
                print(title_str)
                if row_m != 0 and col_m != nmodels - 1:
                    xu.ax_clr(axfs=_fcfg.ax_fs)
                    xu.rotate_labels(which_ax='x', axfs=_fcfg.ax_fs, rot=90)
                elif row_m == 0 and col_m != nmodels - 1:
                    xu.ax_clrX(axfs=_fcfg.ax_fs)
                    xu.rotate_labels(which_ax='x', axfs=_fcfg.ax_fs, rot=90)
                elif col_m == nmodels - 1 and row_m != 0:
                    xu.ax_clrY(axfs=_fcfg.ax_fs)
                    xu.rotate_labels(which_ax='x', axfs=_fcfg.ax_fs, rot=90)
                if row_m == 0 and col_m == nmodels - 1:
                    xu.ax_orig(axfs=_fcfg.ax_fs)
                    xu.rotate_labels(which_ax='x', axfs=_fcfg.ax_fs, rot=90)
                    plt.ylabel('$model_{column}$', fontsize=_fcfg.ax_fs)
                    plt.xlabel('$model_{row}$', fontsize=_fcfg.ax_fs)

            # plot the maps of ratio of models and observation above the diagonal
            if row_m > col_m:
                _ax = plt.axes([
                    _fcfg.x0 + row_m * _fcfg.wp + row_m * _fcfg.xsp, _fcfg.y0 -
                    (col_m * _fcfg.hp + col_m * _fcfg.ysp), _fcfg.wp, _fcfg.hp
                ],
                               projection=ccrs.Robinson(central_longitude=0),
                               frameon=False)
                plot_dat = xu.remove_invalid(dat_row / dat_col,
                                             _fill_val=_fcfg.fill_val)
                _ax.imshow(_get_data_to_plot(plot_dat, _fcfg),
                           norm=matplotlib.colors.BoundaryNorm(
                               cbInfo_ratio['tickBounds'], len(cbInfo_ratio['tickBounds'])),
                           interpolation='none',
                           vmin=cbInfo_ratio['tickBounds'][0],
                           vmax=cbInfo_ratio['tickBounds'][-1],
                           cmap=cbInfo_ratio['colMap'],
                           origin='upper',
                           transform=ccrs.PlateCarree())
                _fix_map(_ax, _fcfg)
            if col_m == 0:
                if row_mod == 'obs':
                    _title_sp = _fcfg.obs_label
                else:
                    _title_sp = row_mod
                plt.title(str(row_m + 1) + ': ' + _title_sp,
                          fontsize=0.809 * _fcfg.ax_fs)
            if row_m == nmodels - 1:
                if col_mod == 'obs':
                    _title_sp = _fcfg.obs_label
                else:
                    _title_sp = col_mod
                _title_sp = str(col_m + 1)
                t_x = _ax.text(1.1,
                               0.5,
                               _title_sp,
                               fontsize=0.809 * _fcfg.ax_fs,
                               va='center',
                               ha='center',
                               transform=_ax.transAxes)

    # plot the colorbar for maps along the diagonal
    y_colo = _fcfg.y0 + _fcfg.hp + _fcfg.cb_off_y
    _axcol_dia = [_fcfg.x_colo_d, y_colo, _fcfg.wcolo, _fcfg.hcolo]
    cb_tit_d = 'ecosystem_carbon_turnover_time (yr)'
    cb = xu.mk_colo_tau(_axcol_dia,
                        cbInfo_diagonal['tickBounds'],
                        cbInfo_diagonal['colMap'],
                        tick_locs=cbInfo_diagonal['ticksLoc'],
                        cbfs=0.86 * _fcfg.ax_fs,
                        cbtitle=cb_tit_d,
                        cbrt=90)

    # plot the colorbar for maps above the diagonal
    y_colo = _fcfg.y0 + _fcfg.hp + _fcfg.cb_off_y
    _axcol_rat = [_fcfg.x_colo_r, y_colo, _fcfg.wcolo, _fcfg.hcolo]
    cb = xu.mk_colo_cont(_axcol_rat,
                         cbInfo_ratio['tickBounds'],
                         cbInfo_ratio['colMap'],
                         cbfs=0.7 * _fcfg.ax_fs,
                         cbrt=90,
                         col_scale='log',
                         cbtitle='ratio ($model_{column}$/$model_{row}$)',
                         tick_locs=cbInfo_ratio['ticksLoc'])
    cb.ax.set_xticklabels(cbInfo_ratio['ticksLab'],
                          fontsize=0.86 * _fcfg.ax_fs,
                          ha='center',
                          rotation=0)

    # save and close the figure
    local_path = cfg['plot_dir']
    png_name = 'global_comparison_matrix_models_' + _fcfg.obs_label + '.png'
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x],
                dpi=450)
    plt.close()

    return 'Plotted full factorial model observation comparison matrix'


def _plot_multimodel_agreement(all_mod_dat, all_obs_dat, cfg):
    """
    makes the maps of bias of multimodel median turnover time with
    multimodel agreement

    Arguments:
        all_mod_dat - dictionary of all model data
        all_obs_dat - dictionary of observed data
        cfg - nested dictionary of metadata

    Returns:
        string- on completion of plotting and saving

    """
    # get the settings for plotting figure
    _fcfg = _get_fig_config(cfg)

    # get the observation data needed to calculate the bias and multimodel agreement
    tau_obs = all_obs_dat['tau_ctotal']['grid']
    tau_obs_5 = all_obs_dat['tau_ctotal_5']['grid']
    tau_obs_95 = all_obs_dat['tau_ctotal_95']['grid']

    # set the information of the colormap used for plotting bias
    cbInfo = _get_ratio_colorbarInfo()

    # calculate the bias of multimodel median turnover time
    models = list(all_mod_dat.keys())
    nmodels = len(models)
    dat_tau_full = np.ones(
        (nmodels, np.shape(tau_obs)[0], np.shape(tau_obs)[1])) * _fcfg.fill_val
    for row_m in range(nmodels):
        row_mod = models[row_m]
        dat_tau = all_mod_dat[row_mod]['grid']
        dat_tau_full[row_m] = dat_tau

    mm_tau = xu.remove_invalid(np.nanmedian(dat_tau_full, axis=0),
                               _fill_val=_fcfg.fill_val)
    mm_bias_tau = mm_tau / tau_obs
    mm_bias_tau = xu.remove_invalid(mm_bias_tau, _fill_val=_fcfg.fill_val)
    # mm_bias_tau[np.isnan(mm_bias_tau)] = 0

    # define figure and main axis to plot the map
    plt.figure(figsize=(5, 3))
    _ax = plt.axes([0.1, 0.1, 0.9, 0.9],
                   projection=ccrs.Robinson(central_longitude=0),
                   frameon=False)

    # plot the data of multimodel bias (=bias of multimodel median turnover time)
    _ax.imshow(_get_data_to_plot(mm_bias_tau, _fcfg),
               norm=matplotlib.colors.BoundaryNorm(cbInfo['tickBounds'],
                                                   len(cbInfo['tickBounds'])),
               interpolation='none',
               vmin=cbInfo['tickBounds'][0],
               vmax=cbInfo['tickBounds'][-1],
               cmap=cbInfo['colMap'],
               origin='upper',
               transform=ccrs.PlateCarree())
    _fix_map(_ax, _fcfg)

    # get the model agreement mask (less than quarter of the model within the
    # observational uncertainty)
    agreement_mask_tau = _get_agreement_mask(dat_tau_full,
                                             tau_obs_5,
                                             tau_obs_95,
                                             nmodels,
                                             nmodel_reject=int(nmodels / 4),
                                             fill_val=_fcfg.fill_val)

    # plot the hatches for uncertainty/multimodel agreement
    lats = all_obs_dat['coords']['latitude']
    lons = all_obs_dat['coords']['longitude']
    latint = abs(lats[1] - lats[0])
    lonint = abs(lons[1] - lons[0])
    x, y = np.meshgrid(lons - lonint / 2, lats - latint / 2)

    _ax.contourf(x,
                 y,
                 agreement_mask_tau,
                 levels=[0, 0.5, 1],
                 alpha=0.,
                 hatches=['', '//////'],
                 linewidth=0.2,
                 transform=ccrs.PlateCarree())

    plt.title('ecosystem_carbon_turnover_time (yr), ' + _fcfg.obs_label +
              ',\n' + 'multimodel bias and agreement',
              fontsize=0.98 * _fcfg.ax_fs)

    # plot colorbar using extraUtils
    _axcol_rat = [0.254, _fcfg['y_colo_single'], 0.6, 0.035]

    cb = xu.mk_colo_cont(_axcol_rat,
                         cbInfo['tickBounds'],
                         cbInfo['colMap'],
                         cbfs=0.8 * _fcfg.ax_fs,
                         cbrt=90,
                         col_scale='log',
                         cbtitle='',
                         tick_locs=cbInfo['ticksLoc'])
    cb.ax.set_xticklabels(cbInfo['ticksLab'],
                          fontsize=0.9586 * _fcfg.ax_fs,
                          ha='center',
                          rotation=0)

    # save and close figure
    local_path = cfg['plot_dir']
    png_name = 'global_multimodelAgreement_turnovertime_' + _fcfg.obs_label + '.png'
    t_x = plt.figtext(0.5, 0.5, ' ', transform=plt.gca().transAxes)
    plt.savefig(os.path.join(local_path, png_name),
                bbox_inches='tight',
                bbox_extra_artists=[t_x],
                dpi=450)
    plt.close()

    return 'Plotted multimodel bias and agreement of turnover time'


def _plot_single_map(_dat, _datglobal, _name, _cfg):
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
    # figure configuration
    _fcfg = _get_fig_config(_cfg)

    # colormap configuration
    cbInfo = _get_dia_colorbarInfo()

    # define the figure and axis
    plt.figure(figsize=(5, 3))
    _ax = plt.axes([0.1, 0.1, 0.9, 0.9],
                   projection=ccrs.Robinson(central_longitude=0),
                   frameon=False)
    # plot data over the map
    plt.imshow(_get_data_to_plot(_dat, _fcfg),
               norm=matplotlib.colors.BoundaryNorm(cbInfo['tickBounds'],
                                                   len(cbInfo['tickBounds'])),
               cmap=cbInfo['colMap'],
               origin='upper',
               vmin=cbInfo['tickBounds'][0],
               vmax=cbInfo['tickBounds'][-1],
               transform=ccrs.PlateCarree())
    _fix_map(_ax, _fcfg)

    # get the data and set the title of the map

    _datMedian = np.nanmedian(xu.remove_invalid(_dat,
                                                _fill_val=_fcfg.fill_val))
    title_str = _fcfg.varName + ": global = " + str(round(
        _datglobal, 2)) + ", median = " + str(round(_datMedian, 2))

    plt.title('ecosystem_carbon_turnover_time (yr), ' + _name + ',\n' +
              title_str,
              fontsize=0.98 * _fcfg.ax_fs)

    # draw the colorbar
    _axcol_dia = [0.254, _fcfg['y_colo_single'], 0.6, 0.035]
    cb = xu.mk_colo_tau(_axcol_dia,
                        cbInfo['tickBounds'],
                        cbInfo['colMap'],
                        tick_locs=cbInfo['ticksLoc'],
                        cbfs=0.86 * _fcfg.ax_fs,
                        cbtitle='',
                        cbrt=90)

    # save the figure
    local_path = _cfg['plot_dir']
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
        allModdat, allObsdat = _get_turnover_data(config)
        _fcfg = _get_fig_config(config)
        _plot_multimodel_agreement(allModdat, allObsdat, config)
        _plot_matrix_map(allModdat, allObsdat, config)
        _plot_single_map(allObsdat['tau_ctotal']['grid'],
                         allObsdat['tau_ctotal']['global'], _fcfg.obs_label,
                         config)
