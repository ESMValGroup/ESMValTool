"""Calculate and plot bias and change for each model."""
import logging
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def log_provenance(filename, ancestors, caption, cfg):
    """Create a provenance record for the output file."""
    provenance = {
        'caption': caption,
        'domains': ['reg'],
        'authors': ['kalverla_peter'],
        'projects': ['isenes3'],
        'ancestors': ancestors,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)


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
        if infodict.get('ensemble') is not None:
            alias = "{project}_{dataset}_{ensemble}".format(**infodict)
        else:
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


def plot_scatter(tidy_df, ancestors, cfg):
    """Plot bias on one axis and change on the other."""
    grid = sns.relplot(
        data=tidy_df,
        x="Bias (RMSD of all gridpoints)",
        y="Mean change (Future - Reference)",
        hue="dataset",
        col="variable",
        facet_kws=dict(sharex=False, sharey=False),
        kind='scatter',
    )

    filename = get_plot_filename('bias_vs_change', cfg)
    grid.fig.savefig(filename, bbox_inches='tight')

    caption = "Bias and change for each variable"
    log_provenance(filename, ancestors, caption, cfg)


def plot_table(dataframe, ancestors, cfg):
    """Render pandas table as a matplotlib figure."""
    fig, axes = plt.subplots()
    pd.plotting.table(axes, dataframe.reset_index().round(2))
    axes.set_axis_off()

    filename = get_plot_filename('table', cfg)
    fig.savefig(filename, bbox_inches='tight')

    caption = "Bias and change for each variable"
    log_provenance(filename, ancestors, caption, cfg)


def plot_htmltable(dataframe, ancestors, cfg):
    """Render pandas table as html output.

    # https://pandas.pydata.org/pandas-docs/stable/user_guide/style.html
    """
    styles = [
        {
            "selector": ".index_name",
            "props": [("text-align", "right")]
        },
        {
            "selector": ".row_heading",
            "props": [("text-align", "right")]
        },
        {
            "selector": "td",
            "props": [("padding", "3px 25px")]
        },
    ]

    styled_table = dataframe\
        .unstack('variable')\
        .style\
        .set_table_styles(styles)\
        .background_gradient(cmap='RdYlGn', low=0, high=1, axis=0)\
        .format("{:.2e}", na_rep="-")\
        .render()

    filename = get_diagnostic_filename('bias_vs_change', cfg, extension='html')
    with open(filename, 'w') as htmloutput:
        htmloutput.write(styled_table)

    caption = "Bias and change for each variable"
    log_provenance(filename, ancestors, caption, cfg)


def main(cfg):
    """Calculate, visualize and save the bias and change for each model."""
    metadata = cfg['input_data'].values()
    grouped_metadata = group_metadata(metadata, 'variable_group')

    biases = {}
    changes = {}
    ancestors = []
    for group, metadata in grouped_metadata.items():

        model_metadata = select_metadata(metadata, tag='model')
        model_data, model_ancestors = load_data(model_metadata)
        ancestors.extend(model_ancestors)

        variable = model_data.name

        if group.endswith('bias'):
            obs_metadata = select_metadata(metadata, tag='observations')
            obs_data, obs_ancestors = load_data(obs_metadata)
            ancestors.extend(obs_ancestors)

            bias = calculate_bias(model_data, obs_data)
            biases[variable] = bias

        elif group.endswith('change'):
            changes[variable] = model_data

        else:
            logger.warning(
                "Got input for variable group %s"
                " but I don't know what to do with it.", group)

    # Combine all variables
    bias = xr.Dataset(biases)
    change = xr.Dataset(changes)
    combined = xr.concat([bias, change], dim='metric')
    combined['metric'] = [
        'Bias (RMSD of all gridpoints)', 'Mean change (Future - Reference)'
    ]

    dataframe = combined.rename(
        tas='Temperature (K)',
        pr='Precipitation (kg/m2/s)',
    ).to_dataframe()
    dataframe.columns.name = 'variable'
    tidy_df = dataframe.stack('variable').unstack('metric')

    plot_scatter(tidy_df, ancestors, cfg)
    plot_table(tidy_df, ancestors, cfg)
    plot_htmltable(tidy_df, ancestors, cfg)

    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
