"""Diagnostic script to plot ENSO metrics lifecycle and seasonality."""

import logging
import os
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np
import iris.plot as iplt
import cartopy.crs as ccrs
from shapely import box
import shapely.vectorized as shp_vect

from esmvalcore.preprocessor import (
    climate_statistics,
    extract_month,
    extract_season,
    anomalies,
)

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# stdout
logger = logging.getLogger(os.path.basename(__file__))


def plot_level1(model_data, obs, metric_values, y_label, title, dtls):
    """Plots ENSO metric data, input data- model and obs."""
    figure = plt.figure(figsize=(10, 6), dpi=300)
    modls = [d[0] for d in metric_values]
    if title in ['ENSO lifecycle']:
        # model first :list
        for datamod, name in zip(model_data[1], modls):
            plt.plot(model_data[0], datamod, label=name)
        plt.plot(*obs, label=f'ref: {dtls[0]}', linestyle='dashdot', linewidth=3, color='black')

    else:
        plt.scatter(0, obs, c=['black'], marker='D')
        plt.scatter(np.ones(len(model_data), np.int8), model_data,
                    c=model_data, cmap='tab10', marker='D')
        plt.xlim(-0.5, 2)
        plt.xticks([])
        logger.info("number of models:%s", len(model_data))

        plt.text(0.75, 0.9, f'* ref: {dtls[0]}', color='black',
                 transform=plt.gca().transAxes)

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

    plt.tight_layout
    return figure

def plot_map1(input_data, rmse, title): #input data is 2 - model and obs

    figure = plt.figure(figsize=(20, 6), dpi=300)

    proj = ccrs.PlateCarree(central_longitude=180)
    figure.suptitle(title)
    i =121

    for label, cube in input_data.items():
        
        ax1 = plt.subplot(i,projection=proj)
        ax1.coastlines()
        cf1 = iplt.contourf(cube, levels=np.arange(-1,1,0.1), extend='both',cmap='RdBu_r')
        ax1.set_title(label)
        gl1 = ax1.gridlines(draw_labels=True, linestyle='--')
        gl1.top_labels = False
        gl1.right_labels = False
        i+=1

    plt.text(0.1, -0.3, f'RMSE: {rmse:.2f} ', fontsize=12, ha='left',
         transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))    
    # Add a single colorbar at the bottom
    cax = plt.axes([0.15,0.08,0.7,0.05])
    cbar = figure.colorbar(cf1, cax=cax, orientation='horizontal', extend='both', ticks=np.arange(-1,1.5,0.5))
    cbar.set_label('regression (°C/°C)')
    logger.info(f"{title}, {label} : metric:{rmse}")
    plt.tight_layout
    
    return figure


def sst_regressed(n34_cube):
    """Regression function for sst_time_series on sst_enso."""
    n34_dec = extract_month(n34_cube, 12)
    leadlagyr = 3
    n34_dec_years = [n34_dec.coord('time').units.num2date(time).year
                     for time in n34_dec.coord('time').points]
    event_years = n34_dec_years[leadlagyr:-leadlagyr]
    # Ensure that the selected years are not the last years
    years_of_interest = []
    for yr in event_years:
        years_of_interest.append([yr - 2, yr - 1, yr, yr + 1, yr + 2, yr + 3])

    n34_selected = []
    for enso_epoch in years_of_interest:
        # Select the data for the current year and append it to n34_selected
        year_enso = iris.Constraint(time=lambda cell: cell.point.year in
                                    enso_epoch)
        cube_2 = n34_cube.extract(year_enso)
        n34_selected.append(cube_2.data.data)

    event_constr = iris.Constraint(time=lambda cell: cell.point.year in
                                   event_years)
    n34_dec_ct = n34_dec.extract(event_constr)

    # 1) linear regression of sst_time_series on sst_enso
    a_data = np.array(n34_selected)
    b_data = n34_dec_ct.data
    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T
    coefs, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)

    return coefs[0]

def lin_regress_matrix(cubeA, cubeB): #array must not contain infs or NaNs
    """
    Calculate the linear regression of cubeA on cubeB using matrix operations.

    Parameters
    ----------
    cubeA: iris.cube.Cube
        The 2D input cube for which the regression is calculated.
    
    cubeB: iris.cube.Cube
        The cube used as the independent variable in the regression.

    Returns
    -------
    iris.cube.Cube
        A new cube containing the slope of the regression for each spatial point.
    """
    # Get data as flattened arrays
    A_data = cubeA.data.reshape(cubeA.shape[0], -1)  # Shape (time, spatial_points)
    B_data = cubeB.data.flatten()  # Shape (time,)
    logger.info("cubes: %s, %s", cubeA.name, cubeB.name)
    # Add intercept term by stacking a column of ones with cubeB
    B_with_intercept = np.vstack([B_data, np.ones_like(B_data)]).T

    # Solve the linear equations using least squares method
    coefs, _, _, _ = np.linalg.lstsq(B_with_intercept, A_data, rcond=None)
    logger.info("%s, %s",cubeA.coords(), cubeA.shape)
    # Extract slopes from coefficients #coefs 1
    slopes = coefs[0].reshape(cubeA.shape[1], cubeA.shape[2])

    # Create a new Iris Cube for the regression results
    result_cube = iris.cube.Cube(slopes, long_name='regression ENSO SSTA',
                                 dim_coords_and_dims=[(cubeA.coord('latitude'), 0),
                                                      (cubeA.coord('longitude'), 1)])

    return result_cube


def compute_enso_metrics(input_pair, dt_ls, var_group, metric):
    """Compute the ENSO metric for the collected preprocessed data."""
    fig = None
    metric_values = []
    model_plot = []
    # input_pair: obs first
    if metric == 'lifecycle':
        obs = sst_regressed(input_pair[0][var_group[0]])
        # model list, calc values
        for dataset in input_pair[1]:
            logger.info("%s, preprocessed cubes:%s, dataset:%s",
                        metric, len(input_pair[1]), dataset)

            model_datasets = {attr['variable_group']:
                              iris.load_cube(attr['filename'])
                              for attr in input_pair[1][dataset]}
            model = sst_regressed(model_datasets[var_group[0]])
            val = np.sqrt(np.mean((obs - model) ** 2))
            metric_values.append((dataset, val))
            model_plot.append(model)

        months = np.arange(1, 73) - 36
        # plot function, xticks, labels as dict/ls
        fig = plot_level1((months, model_plot), (months, obs), metric_values,
                          'Degree C / C', f'ENSO {metric}', dt_ls)

    elif metric == 'seasonality':
        obs = seasonality_calc(input_pair[0][var_group[0]])
        for dataset in input_pair[1]:
            logger.info("%s, preprocessed cubes:%s, dataset:%s",
                        metric, len(input_pair[1]), dataset)

            model_datasets = {attr['variable_group']:
                              iris.load_cube(attr['filename'])
                              for attr in input_pair[1][dataset]}
            model = seasonality_calc(model_datasets[var_group[0]])

            val = abs((model-obs)/obs)*100
            metric_values.append((dataset, val))
            model_plot.append(model)
        fig = plot_level1(model_plot, obs, metric_values,
                          'SSTA std (NDJ/MAM)(°C/°C)',
                          f'ENSO {metric}', dt_ls)

    return metric_values, fig


def compute_telecon_metrics(input_pair, var_group, metric):

    if metric =='pr_telecon':
        title = '{} PR Teleconnection' # both seasons
    elif metric == 'ts_telecon':
        title = '{} SST Teleconnection'

    val, fig = {}, {}
    for seas in ['DJF','JJA']:
        data_values = []
        cubes = {}
        for label, ds in input_pair.items(): #obs 0, mod 1
            preproc = {}
            for variable in var_group:
                cube = extract_season(ds[variable].copy(), seas)
                preproc[variable] = anomalies(cube, period="full")

            regcube = lin_regress_matrix(preproc[var_group[1]], preproc[var_group[0]])
            reg_masked = mask_pacific(regcube)

            data_values.append(reg_masked.data)
            cubes[label] = reg_masked

        val[seas] = np.sqrt(np.mean((data_values[0] - data_values[1]) ** 2))
        fig[seas] = plot_map1(cubes, val[seas], title.format(seas))

    return val, fig 


def seasonality_calc(dset_cube):
    preproc = {}
    for seas in ['NDJ', 'MAM']:
        cube = extract_season(dset_cube, seas)
        cube = climate_statistics(cube, operator="std_dev",
                                  period="full")
        preproc[seas] = cube.data

    return preproc['NDJ']/preproc['MAM']


def mask_pacific(cube):
    region = box(130.,-15.,270.,15) #remove land
    x_p, y_p = np.meshgrid(
        cube.coord(axis="X").points,
        cube.coord(axis="Y").points,
    )

    mask = shp_vect.contains(region, x_p, y_p)
    cube.data.mask = mask
    return cube


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""

    record = {
        'caption': caption,
        'statistics': ['anomaly'],
        'domains': ['eq'],
        'plot_types': ['line'], #map telecon
        'authors': [
            'chun_felicity',
            # 'beucher_romain',
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
    metrics = {'lifecycle': ['tos_lifdur1'],  # 'tos_lifdurdiv2' for diag lvl2
               'seasonality': ['tos_seas_asym'],
               'pr_telecon': ['tos_enso', 'pr_global'],
               'ts_telecon':['tos_enso','tos_global']}

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

        msg = f"{metric}: observation datasets {len(obs)},"
        msg += f"models {pformat(models)}"
        logger.info(msg)

        # obs datasets for each model
        obs_datasets = {datas['variable_group']:
                        iris.load_cube(datas['filename']) for datas in obs}

        # group models by dataset
        model_ds = group_metadata(models, 'dataset', sort='project')

        if metric.endswith('telecon'):
            for dataset in model_ds:
                logger.info(f"{metric}, preprocessed cubes:{len(model_ds)}, dataset:{dataset}")
                dt_files = [ds['filename'] 
                            for ds in obs] + [ds['filename'] 
                                            for ds in model_ds[dataset]]

                model_datasets = {attributes['variable_group']: iris.load_cube(attributes['filename']) 
                                for attributes in model_ds[dataset]}
                input_pair = {obs[0]['dataset']:obs_datasets, dataset:model_datasets}

                values, fig = compute_telecon_metrics(input_pair, var_preproc, metric)

                # save metric for each pair, check not none teleconnection metric value djf, jja
                for seas, val in values.items():
                    metricfile = get_diagnostic_filename('matrix', cfg, extension='csv')
                    with open(metricfile, 'a+') as f:
                        f.write(f"{dataset},{seas}_{metric},{val}\n")
                    
                    prov_record = get_provenance_record(f'ENSO metrics {seas} {metric}',
                                                        dt_files)
                    save_figure(f'{dataset}_{seas}_{metric}', prov_record, cfg,
                                figure=fig[seas], dpi=300)#
        
        else:
            input_pair = [obs_datasets, model_ds]
            
            dt_files = [datas['filename'] for datas in obs]
            for dset in models:
                dt_files.append(dset['filename'])
            prov_record = get_provenance_record(f'ENSO metrics {metric}', dt_files)

            # process function for each metric
            values, fig = compute_enso_metrics(input_pair, [obs[0]['dataset']],
                                            var_preproc, metric)
            save_figure(metric, prov_record, cfg,
                        figure=fig, dpi=300)

            # save metric for each pair, check not none
            metricfile = get_diagnostic_filename('matrix', cfg,
                                                extension='csv')
            with open(metricfile, 'a+') as f:
                for dataset, value in values:
                    f.write(f"{dataset},{metric},{value}\n")


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
