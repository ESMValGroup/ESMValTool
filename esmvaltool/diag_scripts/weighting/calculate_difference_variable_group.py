"""Calculating the difference between two preprocessed files
 to be further used as other variable groups

"""
import logging
import os

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
)

from esmvaltool.diag_scripts.weighting.climwip.io_functions import (
    read_metadata,
    read_model_data,
    read_observation_data,
    log_provenance,
)

logger = logging.getLogger(os.path.basename(__file__))


def _save_data(data: 'xr.DataArray', name: str, cfg: dict,
              ancestors: list):
    """ Save data to netcdf for further use"""

    filename_data = get_diagnostic_filename(
        '%s%s_ANOM' %(name, data.short_name), cfg, extension='nc')
    data.to_netcdf(filename_data)

    caption = '%s%s_ANOM' %(name, data.short_name)
    log_provenance(caption, filename_data, cfg, ancestors)

def main(cfg):
    """Compute the difference between the two given variable groups."""

    models, observations = read_metadata(cfg)
    variable_group = list(models)

    dataset1_model = models[variable_group[0]]
    model_data1, model_data_files1 = read_model_data(dataset1_model)

    dataset2_model = models[variable_group[1]]
    model_data2, model_data_files2 = read_model_data(dataset2_model)
    model_data_files1.extend(model_data_files2)

    diff_model_data = model_data1 - model_data2
    diff_model_data.attrs['short_name'] = dataset1_model[0]['short_name']
    #import ipdb; ipdb.set_trace()
    _save_data(diff_model_data, 'MODELS_', cfg, ancestors=model_data_files1)

    dataset1_obs = observations[variable_group[0]]
    obs_data1, obs_data_files1 = read_observation_data(dataset1_obs)

    dataset2_obs = observations[variable_group[1]]
    obs_data2, obs_data_files2 = read_observation_data(dataset2_obs)
    obs_data_files1.extend(obs_data_files2)

    diff_obs_data = obs_data1 - obs_data2
    diff_obs_data.attrs['short_name'] = dataset1_obs[0]['short_name']

    _save_data(diff_obs_data, 'OBS_', cfg, ancestors=obs_data_files1)

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
