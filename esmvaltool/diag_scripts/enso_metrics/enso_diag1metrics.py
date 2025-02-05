"""Diagnostic script to plot ENSO metrics lifecycle and seasonality."""

import os
import iris
import logging
from pprint import pformat
import numpy as np
import sacpy as scp
import xarray as xr

import matplotlib.pyplot as plt

from esmvalcore.preprocessor import (extract_month, extract_season,
                                     climate_statistics)

from esmvaltool.diag_scripts.shared import (run_diagnostic,
                                            save_figure,
                                            get_diagnostic_filename,
                                            group_metadata,
                                            select_metadata,
                                            )

# stdout
logger = logging.getLogger(os.path.basename(__file__))


def plot_level1(input_data, metricval, y_label, title, dtls): 
    """Plots ENSO metric data, input data- model and obs."""
    figure = plt.figure(figsize=(10, 6), dpi=300)

    if title in ['ENSO lifecycle']:
        # model first
        plt.plot(*input_data[0], label=dtls[0])
        plt.plot(*input_data[1], label=f'ref: {dtls[1]}', color='black')
        plt.text(0.5, 0.95, f"RMSE: {metricval:.2f}", fontsize=12,
                 ha='center', transform=plt.gca().transAxes,
                 bbox={'facecolor'='white', 'alpha'=0.8, 'edgecolor'='none'})

    else:
        plt.scatter(range(len(input_data)), input_data,
                    c=['black','blue'], marker='D')
        # obs first
        plt.xlim(-0.5,2)
        plt.xticks([])
        plt.text(0.75,0.95, f'* {dtls[0]}', color='blue',
                 transform=plt.gca().transAxes)
        plt.text(0.75,0.9, f'* ref: {dtls[1]}', color='black',
                 transform=plt.gca().transAxes)
        plt.text(0.75, 0.8, f"metric(%): {metricval:.2f}", fontsize=12,
                 transform=plt.gca().transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    plt.title(title)  # metric name
    plt.legend()
    plt.grid(linestyle='--')
    plt.ylabel(y_label)

    if title == 'ENSO lifecycle':
        plt.axhline(y=0, color='black', linewidth=2)
        xticks = np.arange(1, 73, 6) - 36  # Adjust for lead/lag months
        xtick_labels = ['Jan', 'Jul'] * (len(xticks) // 2)
        plt.xticks(xticks, xtick_labels)
        plt.yticks(np.arange(-2, 2.5, step=1))

    logger.info("%s : metric:%.2f", dtls[0], metricval)

    return figure


def sst_regressed(n34_cube):  # update to matrix method
    """Regression function for sst_time_series on sst_enso."""
    n34_dec = extract_month(n34_cube, 12)
    n34_dec = xr.DataArray.from_iris(n34_dec)
    n34 = xr.DataArray.from_iris(n34_cube)
    leadlagyr = 3  # rolling window cut off, not include first year
    n34_dec_ct = n34_dec[leadlagyr:-leadlagyr]
    event_years = n34_dec_ct.time.dt.year
    # Ensure that the selected years are not the last 
    # or second last year in the n34 dataset
    years_of_interest_array = []
    
    # Fill the array with the years of interest for each event year 
    for i, year in enumerate(event_years):      
        years_of_interest_array.append([year.values - 2, year.values - 1,
                                        year.values, year.values + 1,
                                        year.values + 2, year.values + 3])
    
    n34_selected = []
    # create sst_time_series
    for i, epoch_years in enumerate(years_of_interest_array):
        # Select the data for the current year and append it to n34_selected
        n34_selected.append(n34.sel(time=n34['time.year'].isin(epoch_years)))

    # 1) linear regression of sst_time_series on sst_enso
    slope = scp.LinReg(n34_dec_ct.values, n34_selected).slope
    return slope


def compute_enso_metrics(input_pair, dt_ls, var_group, metric):
    """Compute the ENSO metric for the collected preprocessed data."""
    val, fig = None, None
    # input_pair: obs first
    if metric == 'lifecycle':
        model = sst_regressed(input_pair[1][var_group[0]])
        obs = sst_regressed(input_pair[0][var_group[0]])
        val = np.sqrt(np.mean((obs - model) ** 2))
        months = np.arange(1, 73) - 36
        # plot function, xticks, labels as dict/ls
        fig = plot_level1([(months,model),(months,obs)], val, 
                          'Degree C / C', f'ENSO {metric}', dt_ls)

    elif metric == 'seasonality':
        data_values = []
        for dset in input_pair:  # obs 0, mod 1
            preproc = {}
            for seas in ['NDJ','MAM']:
                cube = extract_season(dset[var_group[0]], seas)
                cube = climate_statistics(cube, operator="std_dev",
                                          period="full")
                preproc[seas] = cube.data

            data_values.append(preproc['NDJ']/preproc['MAM'])

        val = compute(data_values[1], data_values[0])
        fig = plot_level1(data_values, val, 'SSTA std (NDJ/MAM)(°C/°C)',
                          f'ENSO {metric}', dt_ls)

    return val, fig


def compute(obs, mod):
    """Compute reference metric."""
    return abs((mod-obs)/obs)*100

def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""

    record = {
        'caption': caption,
        'statistics': ['anomaly'],
        'domains': ['eq'],
        'plot_types': ['line'],
        'authors': [
            'chun_felicity',
            # 'sullivan_arnold',
        ],
        'references': [
            'access-nri',
        ],
        'ancestors': ancestor_files,
    }
    return record

def main(cfg):
    """Run ENSO metrics."""

    input_data = cfg['input_data'].values()

    # iterate through each metric and get variable group preprocessed data
    metrics = {'lifecycle':['tos_lifdur1'],  # 'tos_lifdurdiv2' for diag lvl2
               'seasonality': ['tos_seas_asym'],}

    # select twice with project to get obs, iterate through model selection
    for metric, var_preproc in metrics.items():
        logger.info("Metric: %s, group: %s", metric, var_preproc)
        obs, models = [], []
        for var_prep in var_preproc:
            obs += select_metadata(input_data, variable_group=var_prep,
                                   project='OBS')
            obs += select_metadata(input_data, variable_group=var_prep,
                                   project='OBS6')
            models += select_metadata(input_data, variable_group=var_prep,
                                      project='CMIP6')

        msg = f"{metric} : observation datasets {len(obs)}, models {pformat(models)}"
        logger.info(msg)
        
        # list dt_files
        dt_files = []
        for dset in models:
            dt_files.append(dset['filename'])
        prov_record = get_provenance_record(f'ENSO metrics {metric}', dt_files)
        # obs datasets for each model
        obs_datasets = {datas['variable_group']: 
                        iris.load_cube(datas['filename']) for datas in obs}

        # group models by dataset
        model_ds = group_metadata(models, 'dataset', sort='project')

        for dataset in model_ds:
            logger.info("%s, preprocessed cubes:%s, dataset:%s",
                        metric, len(model_ds), dataset)
            
            model_datasets = {attr['variable_group']:
                              iris.load_cube(attr['filename'])
                              for attr in model_ds[dataset]}
            input_pair = [obs_datasets, model_datasets]
            logger.info(pformat(model_datasets))
            # process function for each metric
            value, fig = compute_enso_metrics(input_pair,
                                              [dataset, obs[0]['dataset']],
                                              var_preproc, metric)

            # save metric for each pair, check not none
            if value:
                metricfile = get_diagnostic_filename('matrix', cfg,
                                                     extension='csv')
                with open(metricfile, 'a+') as f:
                    f.write(f"{dataset},{metric},{value}\n")

                save_figure(f'{dataset}_{metric}', prov_record, cfg,
                            figure=fig, dpi=300)
            #clear value,fig
            value = None


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
