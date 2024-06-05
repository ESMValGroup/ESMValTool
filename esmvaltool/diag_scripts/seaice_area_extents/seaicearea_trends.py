"""Diagnostic script to plot minima and maxima trends.

based on code from Anton Steketee's COSIMA cookbook notebook
https://cosima-recipes.readthedocs.io/en/latest/DocumentedExamples
    /SeaIce_Obs_Model_Compare.html
"""

import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure
from esmvaltool.diag_scripts.shared._base import get_plot_filename


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def sea_ice_area(sic, area, coordls):
    """Percent sic is in 0-100. Mulitply by portion so divide by 100."""
    sic = sic / 100
    # valid sic between 0.15 and 1
    return (sic * area).where((sic >= 0.15) * (sic <= 1)).sum(coordls)


def sea_ice_area_obs(xdataset):
    """Compute sea ice area for obs dataset."""
    sic = xdataset.siconc
    area_km2 = xdataset.areacello / 1e6
    result = sea_ice_area(sic, area_km2,
                          ['x', 'y']).to_dataset(name='cdr_area')

    # Theres a couple of data gaps which should be nan
    result.loc[{'time': '1988-01-01'}] = np.nan
    result.loc[{'time': '1987-12'}] = np.nan

    return result.sel(time=slice('1979', '2018'))


def sea_ice_area_model_sh(xdataset):
    """Compute sea ice area for model dataset."""
    sic = xdataset.siconc.where(xdataset.siconc.lat < -20, drop=True)

    area_km2 = xdataset.areacello / 1e6  # area convert to km2

    return sea_ice_area(sic, area_km2, ['i', 'j']).to_dataset(name='si_area')


def min_and_max(dataset):
    """Compute min and max for dataset."""
    def min_and_max_year(yeardata):
        result = xr.Dataset()
        result['min'] = yeardata.min()
        result['max'] = yeardata.max()
        return result
    annual_min_max_ds = dataset.si_area.groupby('time.year'
                                                ).apply(min_and_max_year)
    return annual_min_max_ds


def plot_trend(model_min_max, obs_a, minmax):
    """
    function to plot min or max trend

    Parameters
    ----------
    model_min_max: dictionary of model label and xarray ds with min and max.
    obs_a: xarray of observations area
    minmax: 'min' or 'max'
    """

    figure, _axes = plt.subplots()

    # add note for years change
    for mod_label, model_min_max_dt in model_min_max.items():
        model_min_max_dt[minmax].plot(label=mod_label)

    if minmax == 'max':
        obs_a.cdr_area.groupby('time.year').max().plot(label='Obs CDR')
        plt.title('Trends in Sea-Ice Maxima')
    elif minmax == 'min':
        obs_a.cdr_area.groupby('time.year').min().plot(label='Obs CDR')
        plt.title('Trends in Sea-Ice Minima')

    plt.ylabel('Sea-Ice Area (km2)')

    _ = plt.legend()
    return figure


def main(cfg):
    """Compute sea ice area for each input dataset."""
    data = []

    for dataset in cfg['input_data'].values():
        # data values to iterate
        logger.info("dataset: %s", dataset['long_name'])
        data.append([dataset['filename'], dataset['short_name'],
                    dataset['dataset']])

    inputfiles_df = pd.DataFrame(data, columns=['filename', 'short_name',
                                                'dataset'])
    # sort to ensure order of reading
    inputfiles_df.sort_values(['dataset', 'short_name'], inplace=True)
    logger.info(inputfiles_df[['short_name', 'dataset']])

    min_max = {}

    for filepath, short, data_name in inputfiles_df.itertuples(index=False):
        if data_name == 'NSIDC-G02202-sh':
            if short == 'areacello':
                area_obs = xr.open_dataset(filepath)
            else:
                obs_si = xr.open_dataset(filepath)
        else:  # other models
            if short == 'areacello':
                area_mod = xr.open_dataset(filepath)
                dt_label = data_name
            else:
                mod_si = xr.open_dataset(filepath)

                # make sure correct area with model
                if dt_label == data_name:
                    mod_si['areacello'] = area_mod['areacello']
                    model_area_dt = sea_ice_area_model_sh(mod_si)
                    model_min_max_dt = min_and_max(model_area_dt)
                    # make years in ACCESS model compariable
                    model_min_max_dt['year'] = model_min_max_dt.year + 1652
                    min_max[data_name] = model_min_max_dt
                else:
                    logger.warning("..%s missing a variable?", data_name)

    obs_si['areacello'] = area_obs['areacello']

    provenance = get_provenance_record(inputfiles_df['filename'].to_list())
    for trend_type in ['min', 'max']:
        fig = plot_trend(min_max, sea_ice_area_obs(obs_si), trend_type)
        # Save output
        save_figure(get_plot_filename(f'{trend_type}_trend', cfg),
                    provenance, cfg, figure=fig)


def get_provenance_record(ancestor_files):
    """Build provenance record."""
    record = {
        'ancestors': ancestor_files,
        'authors': [
            'chun_felicity',
        ],
        'caption': 'added 1652 years to model years for comparability',
        'domains': ['shpolar'],
        'plot_types': ['times'],
        'references': [],
        'statistics': ['mean'],
    }
    return record


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
