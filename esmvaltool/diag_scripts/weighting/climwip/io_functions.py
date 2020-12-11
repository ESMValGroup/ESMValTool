"""A collection of input-output functions."""
import logging
import os
from collections import defaultdict
from datetime import datetime

import natsort
import numpy as np
import xarray as xr

from esmvaltool.diag_scripts.shared import ProvenanceLogger, group_metadata

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption: str, ancestors: list):
    """Create a provenance record describing the diagnostic data and plots."""
    record = {
        'caption':
        caption,
        'domains': ['reg'],
        'authors': [
            'kalverla_peter',
            'smeets_stef',
            'brunner_lukas',
            'camphuijsen_jaro',
        ],
        'references': [
            'brunner2019',
            'lorenz2018',
            'knutti2017',
        ],
        'ancestors':
        ancestors,
    }
    return record


def log_provenance(caption: str, filename: str, cfg: dict, ancestors: list):
    """Log provenance info."""
    provenance_record = get_provenance_record(caption, ancestors=ancestors)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info('Output stored as %s', filename)


def read_metadata(cfg: dict) -> tuple:
    """Read the metadata from the config file.

    Returns a two dicts, one for the model data and one for the
    observational data. They are split based on the value of the
    'obs_data' variable in the recipe. The dictionaries are sorted by
    the variable.
    """
    obs_ids = cfg['obs_data']
    if isinstance(obs_ids, str):
        obs_ids = [obs_ids]

    input_data = cfg['input_data'].values()
    projects = group_metadata(input_data, attribute='project')

    observations = defaultdict(list)
    models = defaultdict(list)

    for project, metadata in projects.items():

        for item in metadata:
            variable_group = item['variable_group']

            if project in obs_ids:
                observations[variable_group].append(item)
            else:
                models[variable_group].append(item)

    return models, observations


def make_standard_calendar(xrds: 'xr.Dataset'):
    """Make sure time coordinate uses the default calendar.

    Workaround for incompatible calendars 'standard' and 'no-leap'.
    Assumes yearly data.
    """
    try:
        years = xrds.time.dt.year.values
        xrds['time'] = [datetime(year, 7, 1) for year in years]
    except TypeError:
        # Time dimension is 0-d array
        pass
    except AttributeError:
        # Time dimension does not exist
        pass


def read_input_data(metadata: list,
                    dim: str = 'data_ensemble',
                    identifier_fmt: str = '{dataset}') -> tuple:
    """Load data from metadata.

    Read the input data from the list of given data sets. `metadata` is
    a list of metadata containing the filenames to load. Only returns
    the given `variable`. The datasets are stacked along the `dim`
    dimension. Returns an xarray.DataArray.
    """
    data_arrays = []
    identifiers = []
    input_files = []
    for info in metadata:
        filename = info['filename']
        short_name = info['short_name']
        variable_group = info['variable_group']

        xrds = xr.open_dataset(filename)
        make_standard_calendar(xrds)
        xrda = xrds[short_name]
        xrda = xrda.rename(variable_group)
        data_arrays.append(xrda)

        identifier = identifier_fmt.format(**info)
        identifiers.append(identifier)
        input_files.append(filename)

    diagnostic = xr.concat(data_arrays, dim=dim)
    diagnostic[dim] = identifiers

    # Clean up unnecessary coordinate info
    redundant_dims = np.setdiff1d(diagnostic.coords, diagnostic.dims)
    diagnostic = diagnostic.drop(redundant_dims)

    # Use natural sorting order
    sorting = natsort.natsorted(identifiers, alg=natsort.IC)
    diagnostic = diagnostic.sel(indexers={dim: sorting})

    return diagnostic, input_files


def read_model_data(datasets: list) -> tuple:
    """Load model data from list of metadata."""
    return read_input_data(datasets,
                           dim='model_ensemble',
                           identifier_fmt='{dataset}_{ensemble}_{exp}')
