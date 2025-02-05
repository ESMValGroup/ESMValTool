"""diagnostic script to plot ENSO metrics lifecycle and seasonality

"""

import matplotlib.pyplot as plt
import iris.quickplot as qplt

import iris
import os
import logging
from pprint import pformat
import numpy as np
import sacpy as scp
import xarray as xr

from esmvaltool.diag_scripts.shared import (run_diagnostic, 
                                            save_figure, 
                                            get_diagnostic_filename,
                                            group_metadata,
                                            select_metadata,
                                            )
from esmvalcore.preprocessor import (extract_month, extract_season,
                                     climate_statistics)


#stdout
logger = logging.getLogger(os.path.basename(__file__))

def plot_level1(input_data, metricval, y_label, title, dtls): #input data is 2 - model and obs

    figure = plt.figure(figsize=(10, 6), dpi=300)

    if title in ['ENSO pattern','ENSO lifecycle']:
        # model first 
        plt.plot(*input_data[0], label=dtls[0])
        plt.plot(*input_data[1], label=f'ref: {dtls[1]}', color='black')
        plt.text(0.5, 0.95, f"RMSE: {metricval:.2f}", fontsize=12, ha='center', transform=plt.gca().transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    else:
        plt.scatter(range(len(input_data)), input_data, c=['black','blue'], marker='D')
        # obs first
        plt.xlim(-0.5,2)#range(-1,3,1)) #['model','obs']
        plt.xticks([])
        plt.text(0.75,0.95, f'* {dtls[0]}', color='blue', transform=plt.gca().transAxes)
        plt.text(0.75,0.9, f'* ref: {dtls[1]}', color='black', transform=plt.gca().transAxes)
        plt.text(0.75, 0.8, f"metric(%): {metricval:.2f}", fontsize=12, transform=plt.gca().transAxes,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    plt.title(title) # metric name
    plt.legend()
    plt.grid(linestyle='--')
    plt.ylabel(y_label) #param
    
    if title == 'ENSO lifecycle':
        plt.axhline(y=0, color='black', linewidth=2)
        xticks = np.arange(1, 73, 6) - 36  # Adjust for lead/lag months
        xtick_labels = ['Jan', 'Jul'] * (len(xticks) // 2)
        plt.xticks(xticks, xtick_labels)
        plt.yticks(np.arange(-2,2.5, step=1))

    logger.info(f"{dtls[0]} : metric:{metricval}")

    return figure


def sst_regressed(n34_cube): #for lifecycle
    # params cubes, 
    n34_dec = extract_month(n34_cube, 12)
    n34_dec = xr.DataArray.from_iris(n34_dec)
    n34 = xr.DataArray.from_iris(n34_cube)
    leadlagyr = 3 #rolling window cut off, not include first year
    n34_dec_ct = n34_dec[leadlagyr:-leadlagyr]
    event_years = n34_dec_ct.time.dt.year #
    # Ensure that the selected years are not the last or second last year in the n34 dataset
    years_of_interest_array = []
    
    # Fill the array with the years of interest for each event year 
    for i, year in enumerate(event_years):# 2:-3 for dec        
        years_of_interest_array.append([year.values - 2, year.values - 1, year.values, year.values + 1, year.values + 2, year.values + 3])
    
    n34_selected = []
    for i in range(len(years_of_interest_array)): #creates sst_time_series
        # Select the data for the current year and append it to n34_selected #n34 is not dec month only
        n34_selected.append(n34.sel(time=n34['time.year'].isin(years_of_interest_array[i])))

    # 1) linear regression of sst_time_series on sst_enso
    slope = scp.LinReg(n34_dec_ct.values, n34_selected).slope
    return slope

def compute_enso_metrics(input_pair, dt_ls, var_group, metric):

    # input_pair: obs first
    if metric =='lifecycle':
        model = sst_regressed(input_pair[1][var_group[0]]) #n34_cube
        obs = sst_regressed(input_pair[0][var_group[0]])
        val = np.sqrt(np.mean((obs - model) ** 2))
        months = np.arange(1, 73) - 36 #build tuples?
        # plot function #need xticks, labels as dict/ls
        fig = plot_level1([ (months,model),(months,obs)], val, 'Degree C / C', f'ENSO {metric}', dt_ls)

    elif metric =='seasonality':
        data_values = []
        for ds in input_pair: #obs 0, mod 1
            preproc = {}
            for seas in ['NDJ','MAM']:
                cube = extract_season(ds[var_group[0]], seas)
                cube = climate_statistics(cube, operator="std_dev", period="full")
                preproc[seas] = cube.data

            data_values.append(preproc['NDJ']/preproc['MAM'])

        val = compute(data_values[1], data_values[0])
        fig = plot_level1(data_values, val, 'SSTA std (NDJ/MAM)(°C/°C)',f'ENSO {metric}', dt_ls)

    return val, fig 

def compute(obs, mod):
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

    # iterate through each metric and get variable group, select_metadata, map to function call
    metrics = { 'lifecycle':['tos_lifdur1'], #'tos_lifdurdiv2' for lev2
                'seasonality': ['tos_seas_asym'],}
    
    # select twice with project to get obs, iterate through model selection
    for metric, var_preproc in metrics.items(): #if empty or try
        logger.info(f"{metric},{var_preproc}")
        obs, models = [], []
        for var_prep in var_preproc: #enumerate 1 or 2 length? if 2 append,
            obs += select_metadata(input_data, variable_group=var_prep, project='OBS')
            obs += select_metadata(input_data, variable_group=var_prep, project='OBS6')
            models += select_metadata(input_data, variable_group=var_prep, project='CMIP6')

        # log
        msg = "{} : observation datasets {}, models {}".format(metric, len(obs), pformat(models))
        logger.info(msg)
        
        # list dt_files
        dt_files = []
        for ds in models: #and obs?
            dt_files.append(ds['filename'])
        prov_record = get_provenance_record(f'ENSO metrics {metric}', dt_files)
        # obs datasets for each model
        obs_datasets = {dataset['variable_group']: iris.load_cube(dataset['filename']) for dataset in obs}
        
        # group models by dataset
        model_ds = group_metadata(models, 'dataset', sort='project')        
        # dataset name
        
        for dataset in model_ds:
            logger.info(f"{metric}, preprocessed cubes:{len(model_ds)}, dataset:{dataset}")
            
            model_datasets = {attributes['variable_group']: iris.load_cube(attributes['filename']) 
                              for attributes in model_ds[dataset]}
            input_pair = [obs_datasets, model_datasets]
            logger.info(pformat(model_datasets))
            # process function for each metric
            value, fig = compute_enso_metrics(input_pair, [dataset, obs[0]['dataset']], var_preproc, metric)

            # save metric for each pair, check not none
            if value:
                metricfile = get_diagnostic_filename('matrix', cfg, extension='csv')
                with open(metricfile, 'a+') as f:
                    f.write(f"{dataset},{metric},{value}\n")

                save_figure(f'{dataset}_{metric}', prov_record, cfg, figure=fig, dpi=300)
            #clear value,fig
            value = None


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
