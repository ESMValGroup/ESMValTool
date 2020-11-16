"""Provide shared functions for land carbon cycle diagnostic."""
import os
import iris
import numpy as np

from esmvaltool.diag_scripts.shared import select_metadata


def _apply_common_mask(*args):
    """
    Apply common mask to all arrays passed as argument.

    Argument:
    --------
        arrays

    Return:
    ------
        an array with size of nargs x common size of all input arrays
    """
    nargs = len(args)
    for arg_ in range(nargs):
        _dat = args[arg_]
        vars()['dat_mask' + str(arg_)] = np.ones((np.shape(_dat)))
        vars()['dat_mask_inv' + str(arg_)] = np.ma.masked_invalid(_dat).mask
        vars()['dat_mask' + str(arg_)][vars()['dat_mask_inv' + str(arg_)]] = 0
    dat_mask = vars()['dat_mask0']
    for arg_ in range(nargs):
        dat_mask = dat_mask * vars()['dat_mask' + str(arg_)]
    mask_where = np.ma.getmask(np.ma.masked_less(dat_mask, 1.))
    odat = []
    for arg_ in range(nargs):
        _dat = args[arg_].astype(np.float)
        _dat[mask_where] = np.nan
        odat = np.append(odat, np.ma.masked_invalid(_dat))
    odat = odat.reshape(nargs, _dat.shape[0], _dat.shape[1])
    return odat


def _apply_gpp_threshold(gpp_dat, fig_config):
    """Mask the gpp array below threshold."""
    # converting gC m-2 yr-1 to kgC m-2 s-1
    gpp_thres = fig_config["gpp_threshold"] / (86400.0 * 365 * 1000.)
    gpp_dat = np.ma.masked_less(gpp_dat,
                                gpp_thres).filled(fig_config["fill_value"])
    return gpp_dat


def _get_obs_data_zonal(diag_config):
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
    var_list = diag_config.get('obs_variable')

    input_files = []
    for _var in var_list:
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var,
                                                               perc=5))
        var_list = np.append(var_list, '{var}_{perc:d}'.format(var=_var,
                                                               perc=95))
        obs_filename = (f'{_var}_{{frequency}}_{{source_label}}_'
                        f'{{variant_label}}_{{grid_label}}z.nc'.format(
                            **diag_config['obs_info']))
        input_files = np.append(input_files,
                                os.path.join(obs_dir, obs_filename))

    nvars = len(var_list)
    for v_ind in range(nvars):
        var_obs = var_list[v_ind]
        variable_constraint = _var_name_constraint(var_obs)
        cube = iris.load_cube(input_files, constraint=variable_constraint)
        all_data[var_obs] = cube
    for coord in cube.coords():
        all_data[coord.name()] = coord
    return all_data


def _load_variable(metadata, var_name):
    """
    Load data for the variable listed in metadata of the diagnostic variable.

    Argument:
    --------
        metadata - nested dictionary of metadata

    Return:
    ------
        iris cube of the data
    """
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube


def _remove_invalid(tmp, fill_value=-9999.):
    """
    Remove the invalid non-numeric values from the input array.

    Fill it with fill_value.
    Remove all large and small values with magnitude
    beyond 1e15
    """
    tmp = np.ma.masked_outside(tmp, -1e15, 1e15).filled(fill_value)
    where_nan = np.isnan(tmp)
    tmp[where_nan] = fill_value
    where_inf = np.isinf(tmp)
    tmp[where_inf] = fill_value
    return tmp


def _var_name_constraint(var_name):
    """:mod:`iris.Constraint` using `var_name` of a :mod:`iris.cube.Cube`."""
    return iris.Constraint(cube_func=lambda c: c.var_name == var_name)
