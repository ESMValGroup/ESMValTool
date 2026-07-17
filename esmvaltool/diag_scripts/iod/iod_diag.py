import logging
from pathlib import Path
from pprint import pformat
import matplotlib.pyplot as plt
import iris # type: ignore
import numpy as np
import numpy.ma as ma
from datetime import datetime
import calendar
import scipy # type: ignore
from scipy.stats import linregress # type: ignore
from scipy.signal import welch, detrend # type: ignore
import cartopy.crs as ccrs # type: ignore
import cartopy.feature as cfeature # type: ignore


from esmvaltool.diag_scripts.shared import ( # type: ignore
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot # type: ignore

from basic_functions import get_provenance_record, load_data, to_set, compute_cube_diff

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

def find_peak_years(cube, threshold, start_months):
    # To capture a 3-month running mean over a 4-month period
    # ASON for IOD, NDJF for ENSO
    # Add 'year' coordinate if not already present
    if 'year' not in [coord.name() for coord in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if 'month_number' not in [coord.name() for coord in cube.coords()]:
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')

    data = cube.data
    months = cube.coord('month_number').points
    years = cube.coord('year').points

    positive_years = []
    negative_years = []

    for i in range(len(data) - 2):
        window = data[i:i+3]

        if months[i] in start_months:  
            
            if all(val > abs(threshold) for val in window):
                positive_years.append(years[i])
    
            if all(val < -abs(threshold) for val in window):
                negative_years.append(years[i])

    
    return sorted(set(positive_years)), sorted(set(negative_years))


def plot_ts(cfg, dict, title, output_basename):
    """
    Plot all datasets in a single figure.
    """
    plt.figure(figsize=(10, 5))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(dict)))  # Generate distinct colors

    input_filenames = set()

    for i, (dataset, info) in enumerate(dict.items()):
        cube = info['cube']
        file = info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])

        if 'HadISST' in dataset:
            iris.plot.plot(cube, label=dataset, color='black', linewidth=2)
        else:
            iris.plot.plot(cube, label=dataset, color=colors[i])
    
    plt.xlabel("Time", fontsize=12)
    plt.ylabel("Temperature / $^\circ$C", fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(True)

    provenance_record = get_provenance_record(title, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg)
    logger.info(f"Timeseries plot saved: {output_basename}")
    plt.close()

def plot_std(cfg, dict, title, output_basename):
    """
    Plot all datasets in a single figure.
    """
    months = list(calendar.month_abbr)[1:]
    plt.figure(figsize=(10, 5))
    colors = plt.cm.viridis(np.linspace(0, 1, len(dict)))  # Generate distinct colors

    input_filenames = set()

    for i, (dataset, info) in enumerate(dict.items()):
        cube = info['cube']
        file = info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])

        if 'month_number' not in [coord.name() for coord in cube.coords()]:
            iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
        cube_std = cube.aggregated_by('month_number', iris.analysis.STD_DEV) 
        if 'HadISST' in dataset:
            plt.plot(months, cube_std.data, label=dataset, color='black', linewidth=2)
        else:
            plt.plot(months, cube_std.data, label=dataset, color=colors[i])
    
    plt.xlabel("Time", fontsize=12)
    plt.ylabel("STD / $^\circ$C", fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(True)

    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg)
    logger.info(f"Standard deviation plot saved: {output_basename}")
    plt.close()

def plot_spectra(cfg, dict, nperseg, title, output_basename):
    """
    Plot all datasets in a single figure.
    """

    # Compute the Power Spectral Density (PSD)
    fs = 1  # Monthly data
    nperseg = 12*nperseg  # Segment length (5 years, adjust as needed to capture both annual and decadal variability)

    plt.figure(figsize=(10, 5))
    colors = plt.cm.viridis(np.linspace(0, 1, len(dict)))  # Generate distinct colors

    input_filenames = set()
    
    for i, (dataset, info) in enumerate(dict.items()):

        cube = info['cube']
        file = info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])

        f, Pxx = welch(cube.data, fs=fs, nperseg=nperseg)

        # Convert frequency to period in years
        period_years = 1 / (f * 12)
        
        if 'HadISST' in dataset:
            plt.plot(period_years, Pxx, label=dataset, color='black', linewidth=2)
        else:
            plt.plot(period_years, Pxx, label=dataset, color=colors[i])
    
    plt.xlabel("Period / years", fontsize=12)
    plt.ylabel("Power", fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(True)

    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg)
    logger.info(f"Spectra plot saved: {output_basename}")
    plt.show()

def add_subplot(ax, cube, cb_label, cb_range, cmap, title):
    # Add map features
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True, linewidth=0.8, color='grey', linestyle='--', alpha=1)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 12}
    gl.ylabel_style = {'size': 12}

    # Contour plot
    contour = ax.contourf(
        cube.coord('longitude').points,
        cube.coord('latitude').points,
        cube.data,
        cmap=cmap,  # Change colormap here
        levels=cb_range, extend='both'
    )

    # Title for each subplot
    ax.set_title(title, fontsize=14)

    # Add individual colorbar
    cb = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.07, fraction=0.04)
    cb.set_label(cb_label, fontsize=12)
    cb.ax.tick_params(labelsize=10)

def plot_map(cfg, plot_dict, output_basename):

    input_filenames = set()
    (dataset, info), = plot_dict.items() 
    file = info['filename']
    input_filenames.update(file if isinstance(file, list) else [file])
    cube = info['cube']

    fig, axes = plt.subplots(1, figsize=(10,6), subplot_kw={'projection': ccrs.PlateCarree()})
    add_subplot(axes, cube, 'mm $\mathregular{day^{-1}}$ per std DMI', np.arange(-1.9,2.1,0.2), 'PuOr_r', f'Regression of DMI on precip - {dataset}')
    
    # Adjust layout
    fig.tight_layout()

    save_title = f'{output_basename}_{dataset}'
    provenance_record = get_provenance_record(save_title, list(input_filenames))
    save_figure(save_title, provenance_record, cfg)
    logger.info(f"Map plot saved: {save_title}")
    plt.close()


def plot_double_map(cfg, cube1, pos_count, pos_title, cube2, neg_count, neg_title, dataset, output_basename, prov_files):

    # cube1 is positive events
    # cube2 is negative events

    fig, axes = plt.subplots(2, figsize=(8, 10), subplot_kw={'projection': ccrs.PlateCarree()})

    add_subplot(axes.flat[0], cube1, "SST anomaly / $^\\circ$C", np.arange(-2.5,2.7,0.2), 'RdBu_r', pos_title)
    add_subplot(axes.flat[1], cube2, "SST anomaly / $^\\circ$C", np.arange(-2.5,2.7,0.2), 'RdBu_r', neg_title)

    fig.suptitle(f'Composites for {dataset}', fontsize=16, y=0.98)
    fig.text(0.83, 0.93, f'Event count = {pos_count}')
    fig.text(0.83, 0.46, f'Event count = {neg_count}')
    
    # Adjust layout
    fig.tight_layout()

    save_figure(output_basename, list(prov_files), cfg)
    logger.info(f"Composite plot saved: {output_basename}")
    # Show the figure
    plt.show()
        
def extract_composite_cube(years, start_month, monthly_anoms):
    # Start month of composite in zero-indexing
    # i.e. for IOD, plot ASON, August = 7
    # for ENSO, plot NDJF, Nov = 10

    if len(years) == 0:
        return create_blank_cube(monthly_anoms)

    if 'year' not in [coord.name() for coord in monthly_anoms.coords()]:
        iris.coord_categorisation.add_year(monthly_anoms, 'time', name='year')
    if 'month_number' not in [coord.name() for coord in monthly_anoms.coords()]:
        iris.coord_categorisation.add_month_number(monthly_anoms, 'time', name='month_number')
    year_one = monthly_anoms.coord('year').points[0]

    monthly_inds = [(y-year_one)*12+start_month for y in years]

    slices = [monthly_anoms[i:i+4] for i in monthly_inds]
    merged_cube = iris.cube.CubeList(slices).concatenate_cube()
    mean_cube = merged_cube.collapsed('time', iris.analysis.MEAN)

    return mean_cube

def create_blank_cube(reference_cube):
    """Create a blank cube with the same lat/lon as reference."""
    lat = reference_cube.coord('latitude').points
    lon = reference_cube.coord('longitude').points
    data = ma.masked_all((len(lat), len(lon)))
    
    blank_cube = iris.cube.Cube(
        data,
        long_name='Blank',
        units='',
        dim_coords_and_dims=[
            (reference_cube.coord('latitude'), 0),
            (reference_cube.coord('longitude'), 1)
        ]
    )
    return blank_cube

def get_event_sets(iod_yrs, enso_yrs):
    both_yrs = np.intersect1d(iod_yrs, enso_yrs)
    iod_only_yrs= np.setdiff1d(iod_yrs, both_yrs)
    enso_only_yrs = np.setdiff1d(enso_yrs, both_yrs)
    return both_yrs, iod_only_yrs, enso_only_yrs

def composite_map(cfg, dmi_dict, nino_dict, sst_anom_dict):

    common_datasets = set(dmi_dict.keys()) & set(nino_dict.keys())

    if not common_datasets:
        raise ValueError("No matching datasets found for DMI and Nino 3.4 anomalies.")

    # Ensure HadISST is first in an ordered list
    ordered_datasets = ['HadISST'] + [x for x in common_datasets if x != 'HadISST']

    for dataset in ordered_datasets:
        logger.info(f"Plotting composites for {dataset}")

        dmi_info = dmi_dict[dataset]
        nino_info = nino_dict[dataset]
        sst_info= sst_anom_dict[dataset]

        dmi_cube = dmi_info['cube']
        nino_cube = nino_info['cube']
        sst_cube = sst_info['cube']

        prov_files = to_set(dmi_info['filename']) | to_set(sst_info['filename']) | to_set(nino_info['filename'])

        # Compute standard deviations for thresholds
        dmi_std = dmi_cube.collapsed('time', iris.analysis.STD_DEV).data
        nino_std = nino_cube.collapsed('time', iris.analysis.STD_DEV).data

        # Find event years
        pIOD_yrs, nIOD_yrs = find_peak_years(dmi_cube, dmi_std, (8,9))
        EN_yrs, LN_yrs = find_peak_years(nino_cube, nino_std, (11,12))
        print('nIOD', nIOD_yrs)
        print('LN', LN_yrs)

        pIOD_EN_yrs, pIOD_only_yrs, EN_only_yrs = get_event_sets(pIOD_yrs, EN_yrs)
        nIOD_LN_yrs, nIOD_only_yrs, LN_only_yrs = get_event_sets(nIOD_yrs, LN_yrs)

        # --- Extract composites ---
        composites = {
            'pIOD_all': extract_composite_cube(pIOD_yrs, 7, sst_cube),
            'nIOD_all': extract_composite_cube(nIOD_yrs, 7, sst_cube),
            'EN_all': extract_composite_cube(EN_yrs, 10, sst_cube),
            'LN_all': extract_composite_cube(LN_yrs, 10, sst_cube),
            'pIOD_only': extract_composite_cube(pIOD_only_yrs, 7, sst_cube),
            'nIOD_only': extract_composite_cube(nIOD_only_yrs, 7, sst_cube),
            'EN_only': extract_composite_cube(EN_only_yrs, 10, sst_cube),
            'LN_only': extract_composite_cube(LN_only_yrs, 10, sst_cube),
            'pIOD_both': extract_composite_cube(pIOD_EN_yrs, 7, sst_cube),
            'EN_both': extract_composite_cube(pIOD_EN_yrs, 10, sst_cube),
            'nIOD_both': extract_composite_cube(nIOD_LN_yrs, 7, sst_cube),
            'LN_both': extract_composite_cube(nIOD_LN_yrs, 10, sst_cube),
        }

        # --- Plotting ---
        # IOD

        output_basename_iod = f'IOD_composite_maps_{dataset}'
        output_basename_enso = f'ENSO_composite_maps_{dataset}'

        plot_double_map(cfg, composites['pIOD_all'], len(pIOD_yrs), 'Positive IOD - ASON', 
                        composites['nIOD_all'], len(nIOD_yrs), 'Negative IOD - ASON', 
                        dataset, output_basename_iod, prov_files)

        plot_double_map(cfg, composites['pIOD_only'], len(pIOD_only_yrs), 'Only Positive IOD - ASON', 
                        composites['nIOD_only'], len(nIOD_only_yrs), 'Only Negative IOD - ASON', 
                        dataset, f'only_{output_basename_iod}', prov_files)
        
        plot_double_map(cfg, composites['pIOD_both'], len(pIOD_EN_yrs), 'Positive IOD - ASON - both event years', 
                        composites['nIOD_both'], len(nIOD_LN_yrs), 'Negative IOD - ASON - both event years', 
                        dataset, f'both_events_{output_basename_iod}', prov_files)

        # ENSO

        plot_double_map(cfg, composites['EN_all'], len(EN_yrs), 'El Niño - NDJF', 
                        composites['LN_all'], len(LN_yrs), 'La Niña - NDJF', 
                        dataset, output_basename_enso, prov_files)
        
        plot_double_map(cfg, composites['EN_only'], len(EN_only_yrs), 'Only El Niño - NDJF', 
                        composites['LN_only'], len(LN_only_yrs), 'Only La Niña - NDJF', 
                        dataset, f'only_{output_basename_enso}', prov_files)
                                                                             
        plot_double_map(cfg, composites['EN_both'], len(pIOD_EN_yrs), 'El Niño - NDJF - both events years', 
                        composites['LN_both'], len(nIOD_LN_yrs), 'La Niña - NDJF - both event years', 
                        dataset, f'both_events_{output_basename_enso}', prov_files)                  

def compute_regression(cfg, predictor_ts, predictand_map, output_basename):

    regress_results = {}

    (dataset_ts, info_ts), = predictor_ts.items()
    (dataset_map, info_map), = predictand_map.items()
    
    # Compute DMI for each dataset

    cube_ts_all = info_ts['cube']
    cube_map_all = info_map['cube']
    files = to_set(info_ts['filename']) | to_set(info_map['filename'])

    # Extract chosen months to regress on
    for cube in [cube_ts_all, cube_map_all]:
        if 'month_number' not in [coord.name() for coord in cube.coords()]:
            iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    
    month_constraint = iris.Constraint(month_number=[8,9,10,11])
    cube_ts = cube_ts_all.extract(month_constraint)
    cube_map = cube_map_all.extract(month_constraint)

    # Detrend DMI timeseries
    data_detrend = detrend(cube_ts.data)
    
    # Standardise timeseries (zero mean, unit variance)
    ts_mean = np.mean(data_detrend)
    ts_std = np.std(data_detrend)
    ts_standard = (data_detrend - ts_mean) / ts_std
    logger.info(f"Standardised timeseries for {dataset_ts}")

    cube_data = cube_map.data
    # Reshape data to (time, grid_points) for vectorized regression
    time, lat, lon = cube_data.shape
    grid_points = lat * lon
    data_reshaped = cube_data.reshape(time, grid_points)

    # Step 3: Calculate slope and intercept for each grid point in one go
    # Vectorized calculation
    x_mean = np.mean(ts_standard)
    y_mean = np.mean(data_reshaped, axis=0)
    xy_cov = np.mean((ts_standard[:, np.newaxis] * data_reshaped), axis=0) - x_mean * y_mean
    x_var = np.mean(ts_standard**2) - x_mean**2

    # Regression slope (beta)
    slope = xy_cov / x_var

    # Reshape back to lat, lon
    regression_map = slope.reshape((lat, lon))

    # Create a new cube for the regression map
    reg_cube = iris.cube.Cube(
        regression_map,
        long_name='Regression Index',
        var_name='regression',
        units='',
        dim_coords_and_dims=[(cube_map.coord('latitude'), 0), (cube_map.coord('longitude'), 1)]
    )

    dataset_reg = 'Obs' if dataset_ts == 'HadISST' else dataset_ts

    regress_results[dataset_reg] = {
        'cube': reg_cube,
        'filename': list(files)
    }

    # Save output
    save_string = f"{output_basename}_{dataset_reg}"
    provenance_record = get_provenance_record(save_string, list(files))
    save_data(save_string, provenance_record, cfg, reg_cube)

    return regress_results

def check_item_load(metadata, variable, dict_name):
    if (item := load_data(metadata, variable, get_filenames=True)):
            dict_name.update(item)


def main(cfg):
    """Compute the Dipole Mode Index and plot results for all datasets."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')
    dmi_results, west_data, east_data, nino_data, sst_data, west_anoms, east_anoms, precip_data = {}, {}, {}, {}, {}, {}, {}, {}
    for group_name, group_md in grouped_data.items():
        print(group_name)
        west_anoms = load_data(group_md, 'sst_anomalies_west', get_filenames=True)
        east_anoms = load_data(group_md, 'sst_anomalies_east', get_filenames=True)
        if west_anoms and east_anoms:
            dmi_res = compute_cube_diff(cfg, west_anoms, east_anoms, 'dmi')
            dmi_results.update(dmi_res)
        
        # Load west and east SST absolute data
        if (item := load_data(group_md, 'sst_abs_west', get_filenames=True)):
            west_data.update(item)
        if (item := load_data(group_md, 'sst_abs_east', get_filenames=True)):
            east_data.update(item)   
        if (item := load_data(group_md, 'sst_anomalies_nino', get_filenames=True)):
            nino_data.update(item)  
        if (item := load_data(group_md, 'sst_anomalies_global', get_filenames=True)):
            sst_data.update(item)                
        
        if (item := load_data(group_md, 'precip_anomalies_global', get_filenames=True)):
            precip_data.update(item)


    for dataset, dmi_info in dmi_results.items():

        dmi_single ={dataset: dmi_info}
        if dataset == 'HadISST':
            precip_single = {'NCEP': precip_data['NCEP']}
        else:
            precip_single = {dataset: precip_data[dataset]}

        regress_res = compute_regression(cfg, dmi_single, precip_single, 'dmi_precip_regression_ASON')

        # Plot regression maps
        plot_map(cfg, regress_res, 'regression_map_dmi_precip_ASON')

    composite_map(cfg, dmi_results, nino_data, sst_data)
    
    # Plot results for all datasets
    plot_ts(cfg, west_data, 'West SSTs for all datasets', 'west_sst')
    plot_ts(cfg, east_data, 'East SSTs for all datasets', 'east_sst')
    plot_ts(cfg, dmi_results, 'DMI for all datasets', 'dmi_plot')
    plot_ts(cfg, nino_data, 'Nino 3.4 for all datasets', 'nino_plot')

    plot_std(cfg, west_data, 'Standard deviation of west SSTs', 'west_sst_std')
    plot_std(cfg, east_data, 'Standard deviation of east SSTs', 'east_sst_std')
    plot_std(cfg, dmi_results, 'Standard deviation of DMI', 'dmi_std')
    plot_std(cfg, nino_data, 'Standard deviation of Nino 3.4', 'nino_std')

    plot_spectra(cfg, west_data, 5, 'Power spectra of west SSTs', 'west_sst_spectra')
    plot_spectra(cfg, east_data, 5, 'Power spectra of east SSTs', 'east_sst_spectra')
    plot_spectra(cfg, dmi_results, 12, 'Power spectra of DMI', 'dmi_spectra')
    plot_spectra(cfg, nino_data, 12, 'Power spectra of Nino 3.4', 'nino_spectra')

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

