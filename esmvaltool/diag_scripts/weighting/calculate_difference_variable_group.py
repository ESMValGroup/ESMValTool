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


def calculate_diff(ds1: 'xr.DataSet', ds2: 'xr.DataSet',
                   observations=False) -> tuple:
    """
    Read data and calculate differences.

    Return differences and ancestor files.
    """
    if observations:
        data1, data_files1 = read_observation_data(ds1)
        data2, data_files2 = read_observation_data(ds2)
    else:
        data1, data_files1 = read_model_data(ds1)
        data2, data_files2 = read_model_data(ds2)
    data_files1.extend(data_files2)

    diff = data1 - data2
    diff.attrs['short_name'] = ds1[0]['short_name']
    diff.attrs['units'] = ds1[0]['units']

    return diff, data_files1


def _save_data(data: 'xr.DataArray', name: str, cfg: dict,
               ancestors: list):
    """Save data to netcdf for further use."""
    filename_data = get_diagnostic_filename(
        '%s%s_ANOM' % (name, data.short_name), cfg, extension='nc')
    data.to_netcdf(filename_data)

    caption = '%s%s_ANOM' % (name, data.short_name)
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
        ds1 = models[short_name + '_CLIM']
        ds2 = models[short_name + '_GLOBAL']
        diff, data_files_models = calculate_diff(ds1, ds2)
        _save_data(diff, 'MODELS_', cfg, ancestors=data_files_models)

        obs1 = observations[short_name + '_CLIM']
        obs2 = observations[short_name + '_GLOBAL']
        diff_obs, data_files_obs = calculate_diff(
            obs1, obs2, observations=True)
        _save_data(diff_obs, 'OBS_', cfg, ancestors=data_files_obs)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
