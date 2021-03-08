"""Python example diagnostic."""
import logging
from datetime import datetime
from pathlib import Path

# import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def make_standard_calendar(xrda: 'xr.DataArray'):
    """Make sure time coordinate uses the default calendar.

    Workaround for incompatible calendars 'standard' and 'no-leap'.
    Assumes yearly data.
    """
    try:
        years = xrda.time.dt.year.values
        xrda['time'] = [datetime(year, 7, 1) for year in years]
    except TypeError:
        # Time dimension is 0-d array
        pass
    except AttributeError:
        # Time dimension does not exist
        pass


def load_data(metadata: list):
    """Load all files from metadata into an Xarray dataset.

    ``metadata`` is a list of dictionaries with dataset descriptors.
    """

    data_arrays = []
    identifiers = []
    ancestors = []

    for infodict in metadata:
        alias = infodict['alias']
        input_file = infodict['filename']
        short_name = infodict['short_name']

        xrds = xr.open_dataset(input_file)
        xrda = xrds[short_name]

        # Make sure datasets can be combined
        make_standard_calendar(xrda)
        redundant_dims = np.setdiff1d(xrda.coords, xrda.dims)
        xrda = xrda.drop(redundant_dims)

        data_arrays.append(xrda)
        identifiers.append(alias)
        ancestors.append(input_file)

    # Combine along a new dimension
    data_array = xr.concat(data_arrays, dim='dataset')
    data_array['dataset'] = identifiers

    return data_array, ancestors


def area_weighted_mean(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate area mean weighted by the latitude."""
    weights_lat = np.cos(np.radians(data_array.lat))
    means = data_array.weighted(weights_lat).mean(dim=['lat', 'lon', 'time'])

    return means


def calculate_bias(model_data: 'xr.DataArray',
                   obs_data: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate area weighted RMSD with respect to (mean of) observations."""
    if len(obs_data['dataset']) > 1:
        obs_data = obs_data.mean(dim='dataset')
    else:
        obs_data = obs_data.squeeze()

    diff = model_data - obs_data
    bias = area_weighted_mean(diff**2)**0.5

    bias.attrs = model_data.attrs
    return bias


def main(cfg):
    metadata = cfg['input_data'].values()
    grouped_metadata = group_metadata(metadata, 'variable_group')

    biases = {}
    changes = {}
    for group, metadata in grouped_metadata.items():

        model_metadata = select_metadata(metadata, tag='model')
        model_data, model_ancestors = load_data(model_metadata)
        variable = model_data.name

        if group.endswith('bias'):
            obs_metadata = select_metadata(metadata, tag='observations')
            obs_data, obs_ancestors = load_data(obs_metadata)

            bias = calculate_bias(model_data, obs_data)
            biases[variable] = bias

        elif group.endswith('change'):
            changes[variable] = model_data.mean('time')

        else:
            logger.warning(
                "Got input for variable group %s"
                " but I don't know what to do with it.", group)

    # Combine all variables
    bias = xr.Dataset(biases)
    change = xr.Dataset(changes)
    combined = xr.concat([bias, change], dim='metric')
    combined['metric'] = ['bias', 'change']
    dataframe = combined.to_dataframe()
    tidy_df = dataframe.reset_index()

    # TODO: add time average for change to the diagnostic, and add reference
    # period for

    # We'll add some extra code to make sure the output is saved and tracked
    # provenance = {'caption': f"Climate model performance and change info for
    # impact modellers", 'domains': ['reg'], 'authors': ['kalverla_peter'],
    # 'references': ['acknow_project'], 'ancestors': [input_file],
    # }

    return tidy_df


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
