import logging
from pathlib import Path
from pprint import pformat
import matplotlib.pyplot as plt
import iris # type: ignore
import numpy as np
import numpy.ma as ma

from basic_functions import (get_provenance_record, iso_depth_4d, load_and_update_dict, 
                             load_data)

from esmvaltool.diag_scripts.shared import ( # type: ignore
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot # type: ignore

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)



def extract_eq_IO(cube):
    """
    Extract the IO equatorial region. given an global equatorial band.
    Fold in the NEMO grid at ~73E makes direct extraction not possible.
    """
    
    # Average over the extracted +-2 degrees region
    t20d_mean = np.nanmean(cube.data, axis=1)

    # Extract t20d along equatorial IO (45-100E)
    # Split in grid in ORCA at 73E therefore need to rearrange values
    lon_points = cube[0].coord('longitude').points[0]
    print(cube.coord('longitude').points.shape)
    print(lon_points)
    io_inds = np.squeeze(np.where((lon_points >= 45) & (lon_points <= 100)))
    io_lons = lon_points[io_inds]
    io_t20d = t20d_mean[:,io_inds]

    sorted_t20d = []

    for i in range(4):
        sorted_lons, sorted_season = zip(*sorted(zip(io_lons, io_t20d[i])))
        sorted_t20d.append(sorted_season)

    return sorted_t20d    

def get_iso_data(cfg, group_md, variable):
    """
    Compute T20D from 4D temperature fields and return dataset-mean cube dictionary.
    """
    logger.info(f"Processing iso data for variable: {variable}")
    iso_results = {}
    ISO_LEVEL = 20

# Load temperatures using load_data
    temp_data = load_data(group_md, variable, get_filenames=True)

    if not temp_data or not list(temp_data.keys()):
        logger.warning("No temperature data found. Skipping.")
        return None  

    (dataset, info), = temp_data.items()
    cube = info['cube']
    input_filenames = set()
    file = info['filename']
    input_filenames.update(file if isinstance(file, list) else [file])
    logger.debug(f"Selected dataset: {dataset}")

    t20d_cube = iso_depth_4d(cube, ISO_LEVEL, 'season_number')
    t20d_cube.data = np.ma.masked_invalid(t20d_cube.data)
    lon_coord = t20d_cube.coord('longitude')
    # if lon_coord.ndim == 2:
    #     print(f"{dataset} has 2D longitude coordinates, applying extract_eq_IO.")
    #     eq_iso_data = extract_eq_IO(t20d_cube)
    # else:
    #     print(f"{dataset} has 1D longitude coordinates, applying longitude constraint.")
    #     # Extract longitude range
    #     lon_constraint = iris.Constraint(longitude=lambda lon: 45 <= lon <= 100)
    #     extracted = t20d_cube.extract(lon_constraint)
    # Average over latitude to match 2D behavior
    eq_iso_mean = t20d_cube.collapsed('latitude', iris.analysis.MEAN)
    # Convert to list format matching extract_eq_IO output
    eq_iso_data = [eq_iso_mean[i].data for i in range(len(eq_iso_mean.coord('season_number').points))]

    iso_results[dataset] = {
                'cube': eq_iso_data,
                'longitude': eq_iso_mean.coord('longitude').points,
                'filename': file
            }

    # Save output
    output_basename = f"{dataset}_{variable}_t20d"
    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_data(output_basename, provenance_record, cfg, t20d_cube)

    logger.info(f"T20D processed and saved for {dataset}")
    return iso_results


def prepare_wind_data(wind_dict):
    """Collapse raw eq_winds cubes (season_number, latitude, longitude) to
    per-season longitude profiles, matching the format used for iso_results.
    """
    processed = {}
    for dataset, info in wind_dict.items():
        cube = info['cube']
        file = info['filename']
        mean_cube = cube.collapsed('latitude', iris.analysis.MEAN)
        n_seasons = len(mean_cube.coord('season_number').points)
        processed[dataset] = {
            'cube': [mean_cube[i].data for i in range(n_seasons)],
            'longitude': mean_cube.coord('longitude').points,
            'filename': file,
        }
    return processed


def _plot_panel(ax, cfg, plot_dict, season_idx, ylabel, obs_name,
                invert_yaxis=False):
    """Plot one variable's seasonal longitude profile onto a shared axis.

    Returns the common longitude array and the set of input filenames used,
    for axis formatting and provenance tracking by the caller.
    """
    highlight_colors = [
        'tab:orange',
        'tab:red',
        'tab:green',
        'tab:brown',
        'tab:pink',
        'tab:olive',
        'tab:gray',
    ]
    highlight_idx = 0
    multimodel_lons = None
    multimodel_profiles = []
    input_filenames = set()

    for dataset, dict_info in plot_dict.items():
        cube = dict_info['cube']
        file = dict_info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])
        lons = np.asarray(dict_info['longitude'], dtype=float)

        model_vals = np.ma.filled(
            np.ma.masked_invalid(cube[season_idx]).astype(float), np.nan)

        if multimodel_lons is None:
            multimodel_lons = lons
        elif not np.array_equal(lons, multimodel_lons):
            raise ValueError(
                f"{dataset} is not on the common 1x1 degree grid "
                "expected by plot_ts")

        if dataset != obs_name:
            multimodel_profiles.append(model_vals)

        if dataset == obs_name:
            color = 'black'
            linewidth = 2.5
            alpha = 1.0
        elif dataset in cfg.get("highlight_datasets", []):
            color = highlight_colors[highlight_idx % len(highlight_colors)]
            highlight_idx += 1
            linewidth = 1.2
            alpha = 0.6
        else:
            # Skip plotting non-highlight model lines while still using them for multimodel stats.
            continue

        ax.plot(
            lons,
            model_vals,
            label=dataset,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
        )

    if multimodel_profiles:
        profiles = np.asarray(multimodel_profiles, dtype=float)
        multimodel_median = np.nanmedian(profiles, axis=0)
        upper = np.nanpercentile(profiles, 75, axis=0)
        lower = np.nanpercentile(profiles, 25, axis=0)

        ax.fill_between(
            multimodel_lons,
            lower,
            upper,
            color='tab:blue',
            alpha=0.2,
            linewidth=0,
            label='Multimodel IQR',
        )

        ax.plot(
            multimodel_lons,
            multimodel_median,
            label='Multimodel median',
            color='tab:blue',
            linewidth=3,
        )

    if invert_yaxis:
        ax.invert_yaxis()
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(True)

    return multimodel_lons, input_filenames


def plot_ts(cfg, thermocline_dict, wind_dict, title, output_basename):
    """
    Plot zonal wind (top) and thermocline depth (bottom) sharing an x-axis.

    All datasets are expected on a common 1x1 degree grid with the
    equatorial Indian Ocean region already extracted.
    """

    seasons = ['DJF','MAM','JJA','SON']

    for n, season in enumerate(seasons):
        logger.info(f"Plotting thermocline and winds for {season}")

        fig, (ax_wind, ax_therm) = plt.subplots(
            2, 1, figsize=(12, 9), sharex=True,
            gridspec_kw={'height_ratios': [1, 1.3]})

        wind_lons, wind_files = _plot_panel(
            ax_wind, cfg, wind_dict, n,
            "1000 hPa zonal wind (ua) / m s$^{-1}$", obs_name='NCEP')
        ax_wind.axhline(0, color='k', linewidth=1, linestyle='--')
        therm_lons, therm_files = _plot_panel(
            ax_therm, cfg, thermocline_dict, n,
            "20 $^\\circ$C isotherm depth / m", obs_name='EN4',
            invert_yaxis=True)

        multimodel_lons = therm_lons if therm_lons is not None else wind_lons
        if multimodel_lons is not None:
            lon_start = np.floor(multimodel_lons[0] / 5.0) * 5.0
            lon_end = np.ceil(multimodel_lons[-1] / 5.0) * 5.0
            ticks = np.arange(lon_start, lon_end + 5.0, 5.0)
            ticks = ticks[(ticks >= multimodel_lons[0])
                          & (ticks <= multimodel_lons[-1])]
            ax_therm.set_xlim(multimodel_lons[0], multimodel_lons[-1])
            ax_therm.set_xticks(ticks)
            ax_therm.set_xticklabels([f"{int(t)}$^\\circ$E" for t in ticks])
        ax_therm.set_xlabel("Longitude", fontsize=12)

        handles, labels = ax_therm.get_legend_handles_labels()
        fig.legend(handles, labels, bbox_to_anchor=(1.0, 0.95), loc='upper right',
                   borderaxespad=0.5)
        fig.suptitle(f'{title} - {season}', fontsize=14)
        fig.tight_layout(rect=(0, 0, 0.8, 1))

        save_title = f'{output_basename}_{season}'
        input_filenames = wind_files | therm_files
        provenance_record = get_provenance_record(save_title, list(input_filenames))
        save_figure(save_title, provenance_record, cfg)
        logger.info(f"Scatter plot saved: {save_title}")
        plt.close(fig)


def _plot_bias_panel(ax, cfg, plot_dict, season_idx, obs_name, base_color,
                      ylabel, linestyle='-'):
    """Plot model-minus-obs bias (highlighted models, MMM median, IQR) for
    one variable onto `ax`, which may be a twin axis sharing the x-axis with
    another variable's bias panel. Returns the common longitude array.
    """
    if obs_name not in plot_dict:
        raise ValueError(f"Observational dataset {obs_name} not found in plot_dict")

    highlight_colors = [
        'tab:orange',
        'tab:green',
        'tab:brown',
        'tab:pink',
        'tab:olive',
        'tab:gray',
        'tab:purple',
    ]

    obs_info = plot_dict[obs_name]
    lons = np.asarray(obs_info['longitude'], dtype=float)
    obs_vals = np.ma.filled(
        np.ma.masked_invalid(obs_info['cube'][season_idx]).astype(float), np.nan)

    highlight_idx = 0
    bias_profiles = []
    var_label = ylabel.split(' bias')[0]

    for dataset, dict_info in plot_dict.items():
        if dataset == obs_name:
            continue
        model_lons = np.asarray(dict_info['longitude'], dtype=float)
        if not np.array_equal(model_lons, lons):
            raise ValueError(
                f"{dataset} is not on the common 1x1 degree grid "
                "expected by plot_bias_ts")

        model_vals = np.ma.filled(
            np.ma.masked_invalid(dict_info['cube'][season_idx]).astype(float),
            np.nan)
        bias = model_vals - obs_vals
        bias_profiles.append(bias)

        if dataset in cfg.get("highlight_datasets", []):
            color = highlight_colors[highlight_idx % len(highlight_colors)]
            highlight_idx += 1
            ax.plot(
                lons,
                bias,
                label=f'{dataset}',
                color=color,
                linewidth=1.2,
                alpha=0.7,
                linestyle=linestyle,
            )

    if bias_profiles:
        profiles = np.asarray(bias_profiles, dtype=float)
        median_bias = np.nanmedian(profiles, axis=0)
        upper = np.nanpercentile(profiles, 75, axis=0)
        lower = np.nanpercentile(profiles, 25, axis=0)

        ax.fill_between(
            lons,
            lower,
            upper,
            color=base_color,
            alpha=0.15,
            linewidth=0,
            label=f'{var_label} IQR',
        )
        ax.plot(
            lons,
            median_bias,
            color=base_color,
            linewidth=3,
            linestyle=linestyle,
            label=f'{var_label} MMM bias',
        )

    ax.set_ylabel(ylabel, fontsize=12)
    ax.axhline(0, color='k', linewidth=1.5, linestyle=':')
    ax.grid(True)

    return lons


def plot_bias_ts(cfg, thermocline_dict, wind_dict, title, output_basename):
    """
    Plot model-minus-obs bias for thermocline depth (top) and zonal wind
    (bottom) sharing an x-axis, showing multimodel median bias, IQR,
    and highlighted models for each variable.
    """

    seasons = ['DJF', 'MAM', 'JJA', 'SON']

    for n, season in enumerate(seasons):
        logger.info(f"Plotting thermocline and wind bias for {season}")

        fig, (ax_therm, ax_wind) = plt.subplots(
            2, 1, figsize=(12, 9), sharex=True,
            gridspec_kw={'height_ratios': [1.3, 1]})

        wind_lons = _plot_bias_panel(
            ax_wind, cfg, wind_dict, n, obs_name='NCEP',
            base_color='tab:red',
            ylabel="ua bias (model - NCEP) / m s$^{-1}$") 
        therm_lons = _plot_bias_panel(
            ax_therm, cfg, thermocline_dict, n, obs_name='EN4',
            base_color='tab:blue', ylabel="T20D bias (model - EN4) / m")

        multimodel_lons = therm_lons if therm_lons is not None else wind_lons
        if multimodel_lons is not None:
            lon_start = np.floor(multimodel_lons[0] / 5.0) * 5.0
            lon_end = np.ceil(multimodel_lons[-1] / 5.0) * 5.0
            ticks = np.arange(lon_start, lon_end + 5.0, 5.0)
            ticks = ticks[(ticks >= multimodel_lons[0])
                          & (ticks <= multimodel_lons[-1])]
            ax_wind.set_xlim(multimodel_lons[0], multimodel_lons[-1])
            ax_wind.set_xticks(ticks)
            ax_wind.set_xticklabels([f"{int(t)}$^\\circ$E" for t in ticks])
        ax_wind.set_xlabel("Longitude", fontsize=12)

        handles1, labels1 = ax_therm.get_legend_handles_labels()
        handles2, labels2 = ax_wind.get_legend_handles_labels()
        
        # Deduplicate labels while preserving order
        seen = set()
        all_handles, all_labels = [], []
        for h, l in zip(handles1 + handles2, labels1 + labels2):
            if l not in seen:
                all_handles.append(h)
                all_labels.append(l)
                seen.add(l)
        
        if all_handles:
            fig.legend(all_handles, all_labels,
                       bbox_to_anchor=(0.98, 1), loc='upper right',
                       fontsize=10)
        fig.suptitle(f'{title} - {season}', fontsize=14)
        fig.tight_layout(rect=(0, 0, 0.78, 0.98))

        save_title = f'{output_basename}_{season}'
        input_filenames = set()
        for plot_dict in (thermocline_dict, wind_dict):
            for info in plot_dict.values():
                file = info['filename']
                input_filenames.update(file if isinstance(file, list) else [file])
        provenance_record = get_provenance_record(save_title, list(input_filenames))
        save_figure(save_title, provenance_record, cfg)
        logger.info(f"Bias plot saved: {save_title}")
        plt.close(fig)


def main(cfg):
    """Compute the 20degree isotherm along the equatorial Indian Ocean and plot results for all datasets."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')
    iso_results, eq_winds = {}, {}
    for group_name, group_md in grouped_data.items():
        iso_res = get_iso_data(cfg, group_md, 'eq_temps')
        if iso_res:
            iso_results.update(iso_res)
        load_and_update_dict(group_md, "eq_winds", eq_winds)
    
    logger.info("Thermocline calculated, now plotting.")
    wind_results = prepare_wind_data(eq_winds)
    # Plot results for all datasets
    plot_ts(cfg, iso_results, wind_results, 'Variation along the equatorial Indian Ocean', 'eq_IO_therm_winds')
    plot_bias_ts(cfg, iso_results, wind_results, 'Model bias relative to observations', 'eq_IO_therm_winds_bias')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

