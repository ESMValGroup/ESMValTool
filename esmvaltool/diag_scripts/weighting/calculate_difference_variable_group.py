"""Calculating the difference between two preprocessed files."""

import logging
import os
import xarray as xr

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    run_diagnostic,
)

from esmvaltool.diag_scripts.weighting.climwip.io_functions import (
    read_metadata,
    read_model_data,
    read_observation_data,
    log_provenance,
)

logger = logging.getLogger(os.path.basename(__file__))


def calculate_diff(data_clim: list, data_glob: list,
                   observations=False) -> tuple:
    """
    Read data and calculate differences.

    Return differences and ancestor files.
    """

    errmsg = '{}_{} not found but needed for anomaly calculation!'
    if not data_clim:
        raise ValueError(errmsg.format(data_clim[0]['short_name'], 'CLIM'))
    if not data_glob:
        raise ValueError(errmsg.format(data_glob[0]['short_name'], 'GLOBAL'))

    if observations:
        data_clim_read, data_files_clim = read_observation_data(data_clim)
        data_glob_read, data_files_glob = read_observation_data(data_glob)
    else:
        data_clim_read, data_files_clim = read_model_data(data_clim)
        data_glob_read, data_files_glob = read_model_data(data_glob)
    data_files_clim.extend(data_files_glob)

    diff = data_clim_read - data_glob_read

    diff.attrs['short_name'] = data_clim[0]['short_name']
    diff.attrs['units'] = data_clim[0]['units']

    return diff, data_files_clim


def _save_data(data: 'xr.DataArray', name: str, cfg: dict,
               ancestors: list):
    """Save data to netcdf for further use."""
    varn_new = f'{data.short_name}_ANOM'
    filename_data = get_diagnostic_filename(
        '%s%s' % (name, varn_new), cfg, extension='nc')
    data.to_dataset(name=varn_new).to_netcdf(filename_data)

    caption = '%s%s' % (name, varn_new)
    log_provenance(caption, filename_data, cfg, ancestors)


def main(cfg):
    """Compute the difference between the two given variable groups."""
    models, observations = read_metadata(cfg)
    variable_groups = list(models)
    varnames = list()
    for variable_group in variable_groups:
        varname = variable_group.split("_")[0]
        varnames.append(varname)
    short_names = set(varnames)

    for short_name in short_names:
        datasets_mod_clim = models[short_name + '_CLIM']
        datasets_mod_global = models[short_name + '_GLOBAL']
        diff, data_files_models = calculate_diff(
            datasets_mod_clim, datasets_mod_global)
        _save_data(diff, 'MODELS_', cfg, ancestors=data_files_models)

        obs_clim = observations[short_name + '_CLIM']
        obs_glob = observations[short_name + '_GLOBAL']
        diff_obs, data_files_obs = calculate_diff(
            obs_clim, obs_glob, observations=True)
        _save_data(diff_obs, 'OBS_', cfg, ancestors=data_files_obs)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
