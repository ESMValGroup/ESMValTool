"""Surface concentration of trace gases diagnostics.

The diagnostics rely on comparing model data to ground-based NOAA GML surface
observations of trace gases.

Some classes and functions largely adapted from the AOD-AERONET diagnostic
at /diag_scripts/aerosols/aod_aeronet_assess.py.
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import iris.analysis
import iris.coords
import iris.plot as iplt
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy
from cf_units import Unit
from esmvalcore.preprocessor import climate_statistics
from matplotlib import colors, gridspec
from numpy import ma

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)
from esmvaltool.diag_scripts.surface_trace_gases.utils_surface_trace_gases import (
    TRACE_GASES_FACTOR,
    TRACE_GASES_UNITS,
    _aggregate_model_stats,
    _colocate_obs_model,
    _extract_pt,
    _latitude_weights,
    _quick_fix_cube,
    _setup_growth_cube,
)

logger = logging.getLogger(Path(__file__).stem)
fontsizedict = {"title": 25, "axis": 20, "legend": 18, "ticklabel": 18}

PLOT_PARAM = {
    "ch4": {
        "cmap_step": 5,
        "cmap_width": 40,
        "scatter_step": 20,
    },
    "co2": {
        "cmap_step": 1,
        "cmap_width": 6,
        "scatter_step": 50,
    },
    "n2o": {
        "cmap_step": 10,
        "cmap_width": 10,
        "scatter_step": 50,
    },
}
COLORS_MARKERS = ["cornflowerblue", "lightsteelblue", "royalblue"]
MARKERS = ["o", "^", "s"]
LATITUDE_TITLES = [
    r"Latitudes 60$^\circ$N - 90$^\circ$N",
    r"Latitudes 30$^\circ$N - 60$^\circ$N",
    r"Latitudes 30$^\circ$S - 30$^\circ$N",
    r"Latitudes 90$^\circ$S - 30$^\circ$S",
]
LATITUDE_RANGES = {
    "60N - 90N": (60, 90),
    "30N - 60N": (30, 60),
    "30S - 30N": (-30, 30),
    "90S - 30S": (-90, -30),
}


def get_provenance_record(
    ancestors: str | list,
    model_name: str | list,
    project: str | list,
    experiment: str | list,
    timerange: str,
    trace_gas: str,
    obs_name: str,
    var: str,
) -> dict:
    """Create a provenance record describing the diagnostic data and plot.

    Parameters
    ----------
    ancestors : str, list, dict
        List of ancestor files.
    project : str
        Project facet.
    model_name : str, list
        Model facet.
    experiment : str
        Experiment facet.
    timerange : str
        Timerange facet as start_year/end_year.
    trace_gas : str
        Trace gas.
    obs_name : str
        Name of the observational dataset.
    var : str
        Name of the plot.

    Returns
    -------
    record : dict
        Dictionary containing the provenance record.
    """
    captions = {
        "seas_maps_map": f"Seasonal maps of the {trace_gas.upper()} surface concentration "
        f"from {model_name} overlaid with station data points.",
        "seas_maps_scatter": f"Scatter plot that compares the {trace_gas.upper()} surface "
        f"concentration at the station locations of {obs_name} with the "
        f"model data from {model_name} extracted at these same locations "
        "for each season.",
        "timeserie_lat": f"Plots of {trace_gas.upper()} surface concentration time series "
        f"and seasonal cycle per latitude regions for {model_name} "
        f"and {obs_name}.",
        "sensitivity_ampl_growth": f"Plots of {trace_gas.upper()} surface concentration amplitude, "
        "growth, and corresponding sensitivity per latitude regions for "
        f"{model_name} and {obs_name}.",
        "taylor_diag": f"Taylor diagram of {trace_gas.upper()} surface concentration for "
        f"the observations {obs_name} and models {model_name}.",
    }
    plot_type = {
        "seas_maps_map": "map",
        "seas_maps_scatter": "scatter",
        "timeserie_lat": "times",
        "sensitivity_ampl_growth": "times",
        "taylor_diag": "taylor",
    }
    return {
        "caption": captions[var],
        "model": model_name,
        "project": project,
        "experiment": experiment,
        "timerange": timerange,
        "plot_type": plot_type[var],
        "authors": ["lenhardt_julien"],
        "ancestors": ancestors,
    }


def plot_trace_gas_mod_obs(
    fig,
    ax,
    md_data,
    obs_data,
    trace_gas_obs_cube,
    plot_dict,
):
    """Plot trace gas surface concentration.

    The function plots a contour for the model data overlaid
    with NOAA GML surface flask observations climatology.

    Parameters
    ----------
    fig: matplotlib.figure.Figure
        Figure object.
    ax: matplotlib.axes.Axis
        Axes of the figure.
    md_data : iris.cube.Cube
        Model trace gas surface concentration as a cube with latitude
        and longitude coordinates.
    obs_data : list
        Observations of trace gas surface concentration from each
        NOAA GML station.
    trace_gas_obs_cube : iris.cube.Cube
        Holds information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    plot_dict : dict
        Holds plotting settings.
    """
    # Plot model data
    cf_plot = iplt.contourf(
        md_data,
        plot_dict["Levels"],
        colors=plot_dict["Colours"],
        extend="both",
        axes=ax,
    )

    # Latitude and longitude of stations.
    noaa_gml_lats = trace_gas_obs_cube.coord("latitude").points
    noaa_gml_lons = (
        trace_gas_obs_cube.coord("longitude").points + 180
    ) % 360 - 180

    # Loop over stations
    for istn, stn_data in enumerate(obs_data):
        if ma.is_masked(stn_data):
            continue

        # Find position of the observed concentration on the colorscale.
        # np.searchsorted returns index at which inserting new value will
        # maintain a sorted array. We use the color to the left of index.
        cid = np.searchsorted(plot_dict["Levels"], stn_data)
        cid = max(0, cid - 1)  # filter out zero and max when seeking "left"
        cid = min(len(plot_dict["Colours"]) - 1, cid)
        pcol = plot_dict["Colours"][cid]

        # Overlay contourf with observations
        ax.plot(
            noaa_gml_lons[istn],
            noaa_gml_lats[istn],
            color=pcol,
            marker="o",
            markeredgecolor="k",
            markeredgewidth=2,
            markersize=9,
            transform=ccrs.Geodetic(),
        )

    # Decorate the plot
    ax.set_title(plot_dict["Title"], size=24)
    colbar = fig.colorbar(cf_plot, orientation="horizontal", ax=ax)
    colbar.set_label(plot_dict["cb_label"], fontsize=14)
    ax.coastlines(color="#525252")

    # Statistics on plot
    ax.text(
        0.5,
        -0.1,
        (
            f"Global mean {plot_dict['Mean']:.1f} - @Stations mean: mod="
            f"{plot_dict['Stn_mn_md']:.1f} / obs={plot_dict['Stn_mn_obs']:.1f}"
            f" (RMSE={plot_dict['RMS']:.1f})"
        ),
        ha="center",
        va="center",
        size=14,
        transform=ax.transAxes,
    )


def trace_gas_maps(
    model_dataset,
    attribute,
    colocated_datasets,
    clim_seas,
    timerange,
    trace_gas,
):
    """Plot model vs NOAA GML trace gas surface concentration.

    The function generates the plots and returns evaluation metrics.

    Parameters
    ----------
    model_dataset : str
        String containing the model dataset name.
    attribute: dict.
        Dictionary of model dataset entry in the recipe.
    colocated_datasets : dict
        Dictionary containing the preprocessed and colocated datasets.
    trace_gas_obs_cube : iris.cube.Cube
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    clim_seas : list
       Strings to denote climate seasons ["DJF", "MAM", "JJA", "SON"]
    timerange : str
        String containing the time range as "{start_year}/{end_year}".
    trace_gas : str
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.

    Returns
    -------
    figures : list
        Contains figure instances for the seasonal contour plots overlaid with
        observations of surface concentration from NOAA GML.
    fig_scatter : matplotlib.figure.Figure
        The scatter plot comparing model and obs surface concentrations.
    """
    # Set up seasonal contour plots
    figures = []

    # Model data
    trace_gas_at_noaa = colocated_datasets["seas_maps"][model_dataset][
        attribute["alias"]
    ]["colocated"]
    model_data = colocated_datasets["seas_maps"][model_dataset][
        attribute["alias"]
    ]["values"]
    # Obs data
    obs_cube = colocated_datasets["seas_maps"]["obs"]

    # Use global average as center value for the colormap
    grid_areas_weights = iris.analysis.cartography.area_weights(model_data)
    center_cmap_mean = np.round(
        model_data.collapsed(
            ["season_number", "latitude", "longitude"],
            iris.analysis.MEAN,
            weights=grid_areas_weights,
        ).data,
    )
    clevs = list(
        np.arange(
            center_cmap_mean - PLOT_PARAM[trace_gas]["cmap_width"],
            center_cmap_mean + PLOT_PARAM[trace_gas]["cmap_width"] + 1,
            PLOT_PARAM[trace_gas]["cmap_step"],
            dtype=int,
        ),
    )
    clabs = [
        str(lev) if (lev % PLOT_PARAM[trace_gas]["cmap_step"] == 0) else ""
        for lev in clevs
    ]
    cmapr = mpl.colormaps.get_cmap("viridis")
    cmap = colors.ListedColormap(cmapr(np.linspace(0, 1, len(clevs))))
    colours = cmap.colors
    cb_label = (
        f"{trace_gas.upper()} surface concentration "
        f"[{model_data.attributes['unit']}]"
    )

    # Set up the figure for scatter plotting
    fig_scatter = plt.figure(figsize=(10, 10))
    gs_scatter = gridspec.GridSpec(ncols=1, nrows=1)
    ax_scatter = fig_scatter.add_subplot(gs_scatter[0, 0])
    col_scatter = mpl.colormaps.get_cmap("Blues")(np.linspace(0.25, 1, 4))
    mark_scatter = ["o", "^", "s", "P"]
    leg_scatter = []
    min_scatter = 1e5
    max_scatter = 1e-5
    # Add 1:1 line
    line = mlines.Line2D([0, 1], [0, 1], color="#696969")
    transform = ax_scatter.transAxes
    line.set_transform(transform)
    ax_scatter.add_line(line)

    # If all seasons are plotted on the same plot:
    fig_cf, ax_cf = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(22, 16),
        dpi=300,
        subplot_kw={"projection": ccrs.Robinson()},
    )

    # Loop over seasons
    for s, season in enumerate(obs_cube.slices_over("clim_season")):
        # Match NOAA GML obs season with model season number
        model_sn = [c.lower() for c in clim_seas].index(
            season.coord("clim_season").points[0],
        )
        model_season = model_data[model_sn]

        # Generate statistics required - area-weighted mean
        grid_areas = iris.analysis.cartography.area_weights(model_season)
        global_mean = model_season.collapsed(
            ["latitude", "longitude"],
            iris.analysis.MEAN,
            weights=grid_areas,
        )

        # Extract model and obs data for season number (model_sn)
        seas_obs = season.data
        seas_md = np.array([x[model_sn] for x in trace_gas_at_noaa])

        # Match model data with valid obs data
        valid_indices = ~(seas_obs.mask | np.isnan(seas_md))
        valid_obs = seas_obs[valid_indices]
        valid_md = seas_md[valid_indices]

        # Model - obs statistics (diff, model mean and RMS, r2)
        diff = valid_md - valid_obs
        stn_mn_obs = np.mean(valid_obs)
        stn_mn_md = np.mean(valid_md)
        rms = np.sqrt(np.mean(diff**2))
        linreg = scipy.stats.linregress(valid_obs, valid_md)

        # Plot scatter of co-located model and obs data
        ax_scatter.scatter(
            valid_obs,
            valid_md,
            color=col_scatter[model_sn],
            marker=mark_scatter[model_sn],
        )
        min_scatter = np.min(
            [min_scatter, np.nanmin(valid_obs), np.nanmin(valid_md)],
        )
        max_scatter = np.max(
            [max_scatter, np.nanmax(valid_obs), np.nanmax(valid_md)],
        )

        # Legend
        label = f"{clim_seas[model_sn]} = {linreg.rvalue**2:.2f}"
        leg_scatter.append(
            mlines.Line2D(
                [0],
                [0],
                marker=mark_scatter[model_sn],
                color="w",
                label=label,
                markersize=15,
                markerfacecolor=col_scatter[model_sn],
            ),
        )

        # Plot contours overlaid with obs for this run and season
        n_stn = str(len(valid_obs))
        title = (
            f"\nSurface {trace_gas.upper()} concentration "
            + timerange
            + "\n"
            + model_dataset
            + ", "
            + clim_seas[model_sn]
            + r", $N_{stations}=$"
            + n_stn
        )

        # Plot dictionary
        plot_dict = {
            "Mean": global_mean.data,
            "Stn_mn_obs": stn_mn_obs,
            "Stn_mn_md": stn_mn_md,
            "RMS": rms,
            "Levels": clevs,
            "Colours": colours,
            "tick_labels": clabs,
            "cb_label": cb_label,
            "Title": title,
            "Season": clim_seas[model_sn],
            "trace_gas": trace_gas.upper(),
        }
        plot_trace_gas_mod_obs(
            fig_cf,
            ax_cf.flatten()[s],
            model_season,
            seas_obs,
            obs_cube,
            plot_dict,
        )

    figures = fig_cf

    # Decorate the scatter plot
    ax_scatter.set(
        xlim=(min_scatter - 2, max_scatter + 2),
        ylim=(min_scatter - 2, max_scatter + 2),
    )
    ax_scatter.set_xlabel(
        f"NOAA-GML-SURFACE-FLASK-{trace_gas.upper()} "
        f"[{model_data.attributes['unit']}]",
        fontsize=fontsizedict["axis"],
    )
    ax_scatter.set_ylabel(
        f"{model_dataset} {trace_gas.upper()} "
        f"[{model_data.attributes['unit']}]",
        fontsize=fontsizedict["axis"],
    )

    ax_scatter.tick_params(
        axis="both",
        which="major",
        labelsize=fontsizedict["ticklabel"],
        pad=10,
    )

    ax_scatter.set_title(
        f"Model vs observation @stations ($N={n_stn}$)\n"
        f"Surface {trace_gas.upper()} concentration - {timerange}\n",
        fontsize=fontsizedict["title"],
    )

    ax_scatter.legend(
        handles=leg_scatter,
        title="Seasonal R2",
        title_fontsize=fontsizedict["legend"],
        fontsize=fontsizedict["legend"],
    )

    return figures, fig_scatter


def trace_gas_timeserie_zonal(
    model_dataset,
    attribute,
    colocated_datasets,
    timerange,
    trace_gas,
):
    """Plot time series of zonal mean trace gas concentration.

    It uses different latitude slices for model and observational data:
    - 90S - 30S
    - 30S - 30N
    - 30N - 60N
    - 60N - 90N

    Parameters
    ----------
    model_dataset : str
        String containing the model dataset name.
    attribute: dict
        Dictionary of model dataset entry in the recipe.
    colocated_datasets : dict
        Dictionary containing the preprocessed and colocated datasets.
    obs_cube : iris.cube.Cube
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    timerange : str
        String containing the time range as "{start_year}/{end_year}".
    trace_gas : str
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.

    Returns
    -------
    figure : matplotlib.figure.Figure
        Figure for the zonal time series of the trace gas for model data
        and observations of surface concentration from NOAA GML.
        The left column of the figure contains the zonal time series.
        The center-left column contains the zonal scatter plots for model-obs.
        The center-right column contains the zonal seasonal anomaly.
        The right column contains the zonal max-min months.
    """
    # Create figure layout for different plots
    figure = plt.figure(figsize=(26, 26))
    widths = [1.0, 0.85, 1.0, 1.1]
    heights = [1.0, 1.0, 1.0, 1.0]
    gs = gridspec.GridSpec(
        4,
        4,
        left=0.05,
        right=0.95,
        width_ratios=widths,
        height_ratios=heights,
        hspace=0.35,
        wspace=0.25,
    )

    handles_mean = None
    labels_mean = None
    handles_minmax = None
    labels_minmax = None

    start_year, end_year = timerange.split("/")
    years = np.arange(int(start_year), int(end_year) + 1, 1)
    months = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    for l_i, (lat_range, _) in enumerate(LATITUDE_RANGES.items()):
        # Get list of valid obs and model colocated values
        valid_obs = [
            v
            for y in colocated_datasets["latitude_slices"][lat_range]["obs"][
                "yearly"
            ]["valid"].values()
            for v in y
        ]
        valid_md = [
            v
            for y in colocated_datasets["latitude_slices"][lat_range][
                model_dataset
            ][attribute["alias"]]["yearly"]["valid"].values()
            for v in y
        ]
        # Get number of stations
        n_stations = np.max(
            [
                len(y)
                for y in colocated_datasets["latitude_slices"][lat_range][
                    "obs"
                ]["yearly"]["valid"].values()
            ],
        )

        # Plot left column
        # zonal-mean time series and corresponding +/- std range
        model_ts_mean = np.array(
            [
                np.mean(
                    colocated_datasets["latitude_slices"][lat_range][
                        model_dataset
                    ][attribute["alias"]]["yearly"]["valid"][str(y)],
                )
                for y in years
            ],
        )
        model_ts_std = np.array(
            [
                np.std(
                    colocated_datasets["latitude_slices"][lat_range][
                        model_dataset
                    ][attribute["alias"]]["yearly"]["valid"][str(y)],
                )
                for y in years
            ],
        )
        obs_ts_mean = np.array(
            [
                np.mean(
                    colocated_datasets["latitude_slices"][lat_range]["obs"][
                        "yearly"
                    ]["valid"][str(y)],
                )
                for y in years
            ],
        )
        obs_ts_std = np.array(
            [
                np.std(
                    colocated_datasets["latitude_slices"][lat_range]["obs"][
                        "yearly"
                    ]["valid"][str(y)],
                )
                for y in years
            ],
        )
        _ = plt.subplot(gs[l_i, 0])
        plt.fill_between(
            years,
            model_ts_mean - model_ts_std,
            model_ts_mean + model_ts_std,
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.fill_between(
            years,
            obs_ts_mean - obs_ts_std,
            obs_ts_mean + obs_ts_std,
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        plt.plot(
            years,
            model_ts_mean,
            linestyle="solid",
            color=COLORS_MARKERS[0],
            marker=MARKERS[0],
            markersize=10,
            label=f"{model_dataset} @stations: Mean +/- 1 std"
            if l_i == 0
            else None,
        )
        plt.plot(
            years,
            obs_ts_mean,
            linestyle="solid",
            color=COLORS_MARKERS[1],
            marker=MARKERS[1],
            markersize=10,
            label="Observations: Mean +/- 1 std" if l_i == 0 else None,
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} annual mean "
            f"[{TRACE_GASES_UNITS[trace_gas]}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        text = f"No. of sites = {n_stations}"
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)
        # Get legend handles for the first latitude band
        if l_i == 0:
            axes = plt.gcf().axes
            handles_mean, labels_mean = axes[0].get_legend_handles_labels()

        # Plot centre-left column = scatter plot of model-obs
        linreg = scipy.stats.linregress(valid_obs, valid_md)
        _ = plt.subplot(gs[l_i, 1])
        min_axes = np.min(valid_obs + valid_md) - 2
        max_axes = np.max(valid_obs + valid_md) + 2
        plt.scatter(valid_obs, valid_md, color="black")
        plt.axline(
            [0, linreg.intercept],
            slope=linreg.slope,
            linestyle="dashed",
            color="gray",
        )
        plt.plot(
            np.arange(min_axes - 1, max_axes + 1, 1),
            np.arange(min_axes - 1, max_axes + 1, 1),
            color="black",
        )
        plt.xlim([min_axes, max_axes])
        plt.ylim([min_axes, max_axes])
        plt.xlabel(
            f"NOAA GML {trace_gas.upper()} [{TRACE_GASES_UNITS[trace_gas]}]",
            fontsize=20,
        )
        plt.ylabel(
            f"{model_dataset} {trace_gas.upper()} "
            f"[{TRACE_GASES_UNITS[trace_gas]}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", which="major", labelsize=16, pad=10)
        plt.tick_params(axis="both", labelsize=16)
        text = f"$r^{2}$ = {linreg.rvalue**2:.2f}"
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)

        # Plot centre-right column
        # multi-annual mean seasonal variation
        model_seas_anom_mean = np.array(
            [
                np.mean(
                    np.concatenate(
                        [
                            colocated_datasets["latitude_slices"][lat_range][
                                model_dataset
                            ][attribute["alias"]]["anomaly"][str(y)][m]
                            for y in years
                        ],
                    ),
                )
                for m in months
            ],
        )
        model_seas_anom_std = np.array(
            [
                np.std(
                    np.concatenate(
                        [
                            colocated_datasets["latitude_slices"][lat_range][
                                model_dataset
                            ][attribute["alias"]]["anomaly"][str(y)][m]
                            for y in years
                        ],
                    ),
                )
                for m in months
            ],
        )
        obs_seas_anom_mean = np.array(
            [
                np.mean(
                    np.concatenate(
                        [
                            colocated_datasets["latitude_slices"][lat_range][
                                "obs"
                            ]["anomaly"][str(y)][m]
                            for y in years
                        ],
                    ),
                )
                for m in months
            ],
        )
        obs_seas_anom_std = np.array(
            [
                np.std(
                    np.concatenate(
                        [
                            colocated_datasets["latitude_slices"][lat_range][
                                "obs"
                            ]["anomaly"][str(y)][m]
                            for y in years
                        ],
                    ),
                )
                for m in months
            ],
        )
        _ = plt.subplot(gs[l_i, 2])
        plt.fill_between(
            np.arange(0, 12, 1),
            model_seas_anom_mean - model_seas_anom_std,
            model_seas_anom_mean + model_seas_anom_std,
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.fill_between(
            np.arange(0, 12, 1),
            obs_seas_anom_mean - obs_seas_anom_std,
            obs_seas_anom_mean + obs_seas_anom_std,
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        plt.plot(
            np.arange(0, 12, 1),
            model_seas_anom_mean,
            linestyle="solid",
            color=COLORS_MARKERS[0],
            marker=MARKERS[0],
            markersize=10,
        )
        plt.plot(
            np.arange(0, 12, 1),
            obs_seas_anom_mean,
            linestyle="solid",
            color=COLORS_MARKERS[1],
            marker=MARKERS[1],
            markersize=10,
        )
        plt.xticks(np.arange(0, 12, 1), months)
        plt.xlabel("Month", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} seasonal anomaly "
            f"[{TRACE_GASES_UNITS[trace_gas]}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)

        # Plot right column
        # seasonal cycle timing
        seas_max_model = np.array(
            [
                np.argmax(
                    [
                        colocated_datasets["latitude_slices"][lat_range][
                            model_dataset
                        ][attribute["alias"]]["max"][str(y)][m]
                        for m in months
                    ],
                )
                for y in years
            ],
        )
        seas_min_model = np.array(
            [
                np.argmin(
                    [
                        colocated_datasets["latitude_slices"][lat_range][
                            model_dataset
                        ][attribute["alias"]]["min"][str(y)][m]
                        for m in months
                    ],
                )
                for y in years
            ],
        )
        seas_max_obs = np.array(
            [
                np.argmax(
                    [
                        colocated_datasets["latitude_slices"][lat_range][
                            "obs"
                        ]["max"][str(y)][m]
                        for m in months
                    ],
                )
                for y in years
            ],
        )
        seas_min_obs = np.array(
            [
                np.argmin(
                    [
                        colocated_datasets["latitude_slices"][lat_range][
                            "obs"
                        ]["min"][str(y)][m]
                        for m in months
                    ],
                )
                for y in years
            ],
        )
        month_offset = 4
        months_mod = months[month_offset:] + months[:month_offset]
        seas_max_model[seas_max_model < month_offset] += 12
        seas_min_model[seas_min_model < month_offset] += 12
        seas_max_obs[seas_max_obs < month_offset] += 12
        seas_min_obs[seas_min_obs < month_offset] += 12
        _ = plt.subplot(gs[l_i, 3])
        plt.plot(
            years,
            seas_max_model,
            linestyle="solid",
            color=COLORS_MARKERS[0],
            marker=MARKERS[0],
            markersize=10,
            alpha=0.7,
            label=" " if l_i == 0 else None,
        )
        plt.plot(
            years,
            seas_max_obs,
            linestyle="solid",
            color=COLORS_MARKERS[1],
            marker=MARKERS[1],
            markersize=10,
            alpha=0.7,
            label=" " if l_i == 0 else None,
        )
        plt.plot(
            years,
            seas_min_model,
            linestyle="dashed",
            color=COLORS_MARKERS[0],
            marker=MARKERS[0],
            markersize=10,
            fillstyle="none",
            alpha=0.7,
            label=" " if l_i == 0 else None,
        )
        plt.plot(
            years,
            seas_min_obs,
            linestyle="dashed",
            color=COLORS_MARKERS[1],
            marker=MARKERS[1],
            markersize=10,
            fillstyle="none",
            alpha=0.7,
            label=" " if l_i == 0 else None,
        )
        plt.xlabel("Year", fontsize=20)
        plt.yticks(np.arange(month_offset, month_offset + 12, 1), months_mod)
        plt.ylim([month_offset - 1, month_offset + 12])
        plt.ylabel("Month", fontsize=20)
        plt.tick_params(axis="both", labelsize=16)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)
        # Get legend handles for the first latitude band
        if l_i == 0:
            axes = plt.gcf().axes
            handles_minmax, labels_minmax = axes[3].get_legend_handles_labels()
            labels_minmax[0] = ""
            labels_minmax[1] = ""
            labels_minmax[2] = f"{model_dataset} @stations: Max/Min"
            labels_minmax[3] = "Observations: Max/Min"

    figure.legend(
        handles_mean,
        labels_mean,
        loc="upper center",
        bbox_to_anchor=[0.37, 0.08],
        ncol=2,
        fontsize=20,
        borderaxespad=0.01,
        frameon=True,
    )
    figure.legend(
        handles_minmax,
        labels_minmax,
        loc="lower right",
        bbox_to_anchor=[0.95, 0.05],
        columnspacing=0.01,
        handletextpad=0.01,
        ncol=2,
        fontsize=20,
        borderaxespad=0.01,
        frameon=True,
    )

    return figure


def trace_gas_seas_ampl_growth_rate(
    model_dataset,
    attribute,
    colocated_datasets,
    timerange,
    trace_gas,
):
    """Plot amplitude, growth rate and sensitivity between the two quantities.

    It uses different latitude slices for model and observational data:
    - 90S - 30S
    - 30S - 30N
    - 30N - 60N
    - 60N - 90N

    Parameters
    ----------
    model_dataset : str
        String containing the model dataset name.
    attribute: dict
        Dictionary of model dataset entry in the recipe.
    colocated_datasets : dict
        Dictionary containing the preprocessed and colocated datasets.
    obs_cube : iris.cube.Cube
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    timerange : str
        String containing the time range as "{start_year}/{end_year}".
    trace_gas : str
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.

    Returns
    -------
    figure : matplotlib.figure.Figure
        Figure for the zonal time series of the trace gas for model data
        and observations of surface concentration from NOAA GML.
        The left column of the figure contains the zonal amplitude.
        The center column contains the zonal growth rate.
        The right column contains the zonal sensitivity amplitude/growth.
    """
    start_year, end_year = timerange.split("/")
    years = np.arange(int(start_year), int(end_year) + 1, 1)

    # Plots per latitude range
    # Create figure layout for different plots
    figure = plt.figure(figsize=(26, 26))
    widths = [1.0, 1.0, 1.0]
    heights = [1.0, 1.0, 1.0, 1.0]
    gs = gridspec.GridSpec(
        4,
        3,
        left=0.05,
        right=0.95,
        width_ratios=widths,
        height_ratios=heights,
        hspace=0.35,
        wspace=0.25,
    )
    # Loop over latitude ranges
    for l_i, (lat_range, _) in enumerate(LATITUDE_RANGES.items()):
        # Get number of stations
        n_stations = np.max(
            [
                len(y)
                for y in colocated_datasets["latitude_slices"][lat_range][
                    "obs"
                ]["amplitude"]["valid"].values()
            ],
        )

        # Plots
        # Plot left column = amplitude time series
        model_a_mean = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range][
                    model_dataset
                ][attribute["alias"]]["amplitude"]["mean"][str(y)]
                for y in years
            ],
        )
        model_a_std = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range][
                    model_dataset
                ][attribute["alias"]]["amplitude"]["std"][str(y)]
                for y in years
            ],
        )
        obs_a_mean = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range]["obs"][
                    "amplitude"
                ]["mean"][str(y)]
                for y in years
            ],
        )
        obs_a_std = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range]["obs"][
                    "amplitude"
                ]["std"][str(y)]
                for y in years
            ],
        )
        _ = plt.subplot(gs[l_i, 0])
        plt.plot(
            years,
            model_a_mean,
            color=COLORS_MARKERS[0],
            linestyle="solid",
            marker=MARKERS[0],
            markersize=10,
            label=f"{model_dataset} @stations: Mean +/- 1 std"
            if l_i == 0
            else None,
        )
        plt.fill_between(
            years,
            model_a_mean - model_a_std,
            model_a_mean + model_a_std,
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.plot(
            years,
            obs_a_mean,
            color=COLORS_MARKERS[1],
            linestyle="solid",
            marker=MARKERS[1],
            markersize=10,
            label="Observations: Mean +/- 1 std" if l_i == 0 else None,
        )
        plt.fill_between(
            years,
            obs_a_mean - obs_a_std,
            obs_a_mean + obs_a_std,
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        plt.xlim([years[0] - 1, years[-1] + 1])
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            ),
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} seasonal amplitude "
            f"[{TRACE_GASES_UNITS[trace_gas]}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        text = f"No. of sites = {n_stations}"
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)

        # Plot center column = relative growth
        model_g_mean = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range][
                    model_dataset
                ][attribute["alias"]]["growth"]["mean"][str(y)]
                for y in years[1:]
            ],
        )
        model_g_std = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range][
                    model_dataset
                ][attribute["alias"]]["growth"]["std"][str(y)]
                for y in years[1:]
            ],
        )
        obs_g_mean = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range]["obs"][
                    "growth"
                ]["mean"][str(y)]
                for y in years[1:]
            ],
        )
        obs_g_std = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range]["obs"][
                    "growth"
                ]["std"][str(y)]
                for y in years[1:]
            ],
        )
        _ = plt.subplot(gs[l_i, 1])
        plt.plot(
            years[1:],
            model_g_mean,
            color=COLORS_MARKERS[0],
            linestyle="solid",
            marker=MARKERS[0],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            model_g_mean - model_g_std,
            model_g_mean + model_g_std,
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.plot(
            years[1:],
            obs_g_mean,
            color=COLORS_MARKERS[1],
            linestyle="solid",
            marker=MARKERS[1],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            obs_g_mean - obs_g_std,
            obs_g_mean + obs_g_std,
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        plt.xlim([years[0] - 1, years[-1] + 1])
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            ),
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            trace_gas.upper() + r" relative growth [$\%.yr^{-1}$]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)

        # Plot sensitivity between seasonal amplitude and growth
        model_s_mean = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range][
                    model_dataset
                ][attribute["alias"]]["sensitivity"]["mean"][str(y)]
                for y in years[1:]
            ],
        )
        model_s_std = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range][
                    model_dataset
                ][attribute["alias"]]["sensitivity"]["std"][str(y)]
                for y in years[1:]
            ],
        )
        obs_s_mean = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range]["obs"][
                    "sensitivity"
                ]["mean"][str(y)]
                for y in years[1:]
            ],
        )
        obs_s_std = np.array(
            [
                colocated_datasets["latitude_slices"][lat_range]["obs"][
                    "sensitivity"
                ]["std"][str(y)]
                for y in years[1:]
            ],
        )
        _ = plt.subplot(gs[l_i, 2])
        plt.plot(
            years[1:],
            model_s_mean,
            color=COLORS_MARKERS[0],
            linestyle="solid",
            marker=MARKERS[0],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            model_s_mean - model_s_std,
            model_s_mean + model_s_std,
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.plot(
            years[1:],
            obs_s_mean,
            color=COLORS_MARKERS[1],
            linestyle="solid",
            marker=MARKERS[1],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            obs_s_mean - obs_s_std,
            obs_s_mean + obs_s_std,
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        plt.xlim([years[0] - 1, years[-1] + 1])
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            ),
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} sensitivity\n(amplitude/growth)",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        ax = plt.gca()
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        plt.title(LATITUDE_TITLES[l_i], fontsize=28)

    figure.legend(
        loc="upper center",
        bbox_to_anchor=[0.5, 0.08],
        ncol=3,
        fontsize=20,
        borderaxespad=0.01,
    )

    return figure


def plot_taylor_diagram(
    n_stations,
    std_ref,
    std_models,
    corrs,
    labels,
    rmse_contours=None,
    obs_label=None,
):
    """Plot Taylor diagram for models and reference observation.

    Parameters
    ----------
    n_stations : int
        Number of stations in the observational dataset.
    std_ref : float
        Standard deviation of the observational dataset.
    std_models : numpy.array
        Array containing the standard deviations for each model.
    corrs : numpy.array
        Array containing the correlation w/ observations for each model.
    labels : numpy.array
        Array containing the labels for each model
    rmse_contours : numpy.array, list
        Levels at which to plot the RMSE contours.
    obs_label: str
        Label for the observational dataset.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object containing the Taylor diagram.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, polar=True)

    # Limit to top-right quadrant
    ax.set_thetamin(0)
    ax.set_thetamax(90)

    # Convert correlations to angles
    theta = np.arccos(corrs)

    # Plot each model as a point in the diagram
    for r, t, label in zip(std_models, theta, labels, strict=True):
        p = ax.plot(t, r, "o", label=label, markersize=10)
        ax.annotate(
            label,
            xy=(t, r),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=10,
            ha="left",
            va="bottom",
            fontweight="bold",
            c=p[0]._color,
        )

    # Plot the observation point
    ax.plot(
        0,
        std_ref,
        "ks",
        markersize=15,
        label=obs_label if obs_label else "Observations",
    )
    ax.annotate(
        "Obs",
        xy=(0, std_ref),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=10,
        ha="left",
        va="bottom",
        fontweight="bold",
        c="black",
    )

    # Set limits and labels for radial and angular axes
    rmax = max(std_models + [std_ref]) * 1.2
    ax.set_rlim(0, rmax)
    corr_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0]
    corr_angles = np.degrees(np.arccos(corr_ticks))
    ax.set_thetagrids(corr_angles, labels=[f"{c:.2f}" for c in corr_ticks])
    ax.text(
        0.0,
        0.5 * rmax,
        "Standard deviation",
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=12,
    )
    ax.text(
        np.pi / 2,
        1.05 * rmax,
        "Correlation coefficient",
        horizontalalignment="left",
        verticalalignment="bottom",
        fontsize=12,
    )
    fig.suptitle(
        f"Taylor Diagram of model vs. obs {obs_label} ($N=${n_stations})",
        fontsize=15,
    )

    # Add RMSE contours
    if rmse_contours:
        r_range = np.linspace(0, ax.get_rmax(), 100)
        theta_range = np.linspace(0, np.pi / 2, 100)
        r, theta = np.meshgrid(r_range, theta_range, indexing="ij")
        # Compute RMSE on the grid
        rmse_grid = np.sqrt(
            std_ref**2 + r**2 - 2 * std_ref * r * np.cos(theta),
        )
        # Plot the RMSE contours
        contour = ax.contour(
            theta,
            r,
            rmse_grid,
            rmse_contours,
            cmap="copper_r",
            linestyles="--",
        )
        fmt_dict = {level: f"RMSE={level:.2f}" for level in rmse_contours}
        labels = contour.clabel(inline=True, fontsize=8, fmt=fmt_dict)
        for label in labels:
            label.set_fontweight("bold")

    # Add legend
    ax.legend(title="Datasets", loc="upper right", bbox_to_anchor=(1.1, 1.1))

    return fig


def trace_gas_taylor_diag(datasets, obs, trace_gas, config):
    """Plot Taylor diagram for models and reference observation.

    Parameters
    ----------
    datasets : dict
        Dictionary containing the different model entries from the recipe.
    obs : iris.cube.Cube
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    trace_gas : str
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.
    config : dict
        The ESMValTool configuration.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object containing the Taylor diagram.
    """
    # Time coordinate
    time = obs.coord("time", dim_coords=True)
    # Station index coordinate
    ids = obs.coord("Station index (arbitrary)", dim_coords=True).points
    # Latitude coordinate
    lats = obs.coord("latitude", dim_coords=False).points

    # Colocate datasets and observations
    colocated_datasets = {
        "obs": {
            id_s: np.zeros(shape=len(time.points.tolist())) * np.nan
            for id_s in ids
        },
    }
    # Flag for saving obs data only once
    saved_obs = False
    # Looping over model datasets for the colocation process
    for model_dataset, group in datasets.items():
        # Set up dictionary for model data
        colocated_datasets[model_dataset] = {
            attr["alias"]: {} for attr in group
        }
        # Looping over variables in dataset group
        for attr in group:
            # Load cube
            mod_data = _quick_fix_cube(attr["filename"], trace_gas)
            # Set up dictionary for model data
            colocated_datasets[model_dataset][attr["alias"]] = {
                id_s: np.zeros(shape=len(time.points.tolist())) * np.nan
                for id_s in ids
            }
            # Looping over all time steps
            for t, _ts in enumerate(time.units.num2date(time.points)):
                # Extract monthly
                obs_ts = obs[t]
                v_obs, v_mod, v_id = _colocate_obs_model(
                    obs_ts,
                    mod_data["cube"][t],
                    w_id=True,
                )
                # Save outputs to dictionary
                for i, id_s in enumerate(v_id):
                    colocated_datasets[model_dataset][attr["alias"]][id_s][
                        t
                    ] = v_mod[i]
                if not saved_obs:
                    for i, id_s in enumerate(v_id):
                        colocated_datasets["obs"][id_s][t] = v_obs[i]
            saved_obs = True
    # Setup outputs as numpy arrays
    obs_array = np.stack([colocated_datasets["obs"][id_s] for id_s in ids])
    mod_array = []
    for model_dataset, group in datasets.items():
        for attr in group:
            mod_array.append(
                np.stack(
                    [
                        colocated_datasets[model_dataset][attr["alias"]][id_s]
                        for id_s in ids
                    ],
                ),
            )

    # Compute per station statistic for observations
    std_ref_per_station = np.nanstd(obs_array, axis=1)
    weights = _latitude_weights(lats)
    std_ref_global = np.nansum(weights * std_ref_per_station)
    n_stations = obs_array.shape[0]

    # Prepare lists to store statistics
    all_std_models = []
    all_corrs = []
    all_labels = []

    cnt = 0
    for _, group in datasets.items():
        for attr in group:
            std_model, corr = _aggregate_model_stats(
                obs_array,
                mod_array[cnt],
                lats,
            )
            all_std_models.append(std_model)
            all_corrs.append(corr)
            all_labels.append(attr["alias"])
            cnt += 1

    # Define RMSE contours
    rmse_contours = [2, 4, 6] + np.arange(10, 50, 10).tolist()

    fig = plot_taylor_diagram(
        n_stations,
        std_ref_global,
        all_std_models,
        all_corrs,
        all_labels,
        rmse_contours,
        config["observational_dataset"],
    )

    return fig


def preprocess_obs_dataset(obs_dataset, config):
    """Calculate a multiannual seasonal mean surface climatology of the trace gas.

    Observational surface timeseries data from NOAA GML are used to generate
    a multiannual seasonal mean climatology for each NOAA GML station. The
    user sets thresholds (or uses the default settings) to specify the
    amount of valid data required for the climatology. At this stage
    ESMValTool preprocessors are unsuitable for pre-processing the NOAA GML
    observations because of the bespoke nature and application of the
    filtering thresholds.

    Parameters
    ----------
    obs_dataset : dict
        ESMValTool dictionary. Holds meta data for the observational
        trace gas dataset.
    config : dict
        ESMValTool recipe configuration.

    Returns
    -------
    obs_cube : iris.cube.Cube
        Observational data cube w/ filtered stations according to the
        thresholds used for the multi-annual seasonal mean.
    multiannual_seaonal_mean : iris.cube.Cube.
        Preprocessed observational climatology of the trace gas.
    """
    obs_cube = iris.load_cube(obs_dataset[0]["filename"])

    trace_gas = config["trace_gas"]

    # If trace gas is N2O = units are [mol mol-1] according to the CMOR table.
    # Applying the scale factor is necessary for ppb.
    if trace_gas == "n2o":
        obs_cube = TRACE_GASES_FACTOR[trace_gas] * obs_cube

    # Change units to ppm/ppb
    obs_cube.units = Unit(TRACE_GASES_UNITS[trace_gas])

    # Add the clim_season and season_year coordinates.
    iris.coord_categorisation.add_year(obs_cube, "time", name="year")
    iris.coord_categorisation.add_month(obs_cube, "time", name="month")
    iris.coord_categorisation.add_season(obs_cube, "time", name="clim_season")
    iris.coord_categorisation.add_season_year(
        obs_cube,
        "time",
        name="season_year",
    )

    # Set up thresholds for generating the multi annual seasonal mean
    min_mon_per_seas = config["min_mon_per_seas"]
    min_seas_per_year = config["min_seas_per_year"]
    min_seas_per_clim = config["min_seas_per_clim"]

    # Copy obs cube and mask all months with NaNs
    masked_months_obs_cube = obs_cube.copy(
        data=ma.masked_where(np.isnan(obs_cube.data), obs_cube.data),
    )

    # Aggregate (mean) by season.
    # The number of unmasked months per season is counted,
    # and where there are fewer unmasked months than the
    # given threshold, the computed mean is masked.
    annual_seasonal_mean = masked_months_obs_cube.aggregated_by(
        ["clim_season", "season_year"],
        iris.analysis.MEAN,
    )
    annual_seasonal_count = masked_months_obs_cube.aggregated_by(
        ["clim_season", "season_year"],
        iris.analysis.COUNT,
        function=lambda values: ~ma.getmask(values),
    )
    annual_seasonal_mean.data = ma.masked_where(
        annual_seasonal_count.data < min_mon_per_seas,
        annual_seasonal_mean.data,
    )

    # Aggregate (mean) by multi-annual season.
    # The number of unmasked seasons per multi-annual season
    # is counted, and where there are fewer unmasked seasons
    # than the given threshold, the computed multi-annual
    # season is masked.
    multi_annual_seasonal_mean = annual_seasonal_mean.aggregated_by(
        "clim_season",
        iris.analysis.MEAN,
    )
    clim_season_agg_count = annual_seasonal_mean.aggregated_by(
        "clim_season",
        iris.analysis.COUNT,
        function=lambda values: ~ma.getmask(values),
    )
    multi_annual_seasonal_mean.data = ma.masked_where(
        clim_season_agg_count.data < min_seas_per_clim,
        multi_annual_seasonal_mean.data,
    )
    year_agg_count = multi_annual_seasonal_mean.aggregated_by(
        "year",
        iris.analysis.COUNT,
        function=lambda values: ~ma.getmask(values),
    )

    counter = range(
        len(
            multi_annual_seasonal_mean.coord("clim_season").points,
        ),
    )
    valid_station_indices_for_seasons = []
    for iseas in counter:
        multi_annual_seasonal_mean.data[iseas, :] = ma.masked_where(
            year_agg_count.data[0, :] < min_seas_per_year,
            multi_annual_seasonal_mean.data[iseas, :],
        )
        valid_station_indices_for_seasons.append(
            set(np.where(~multi_annual_seasonal_mean.data.mask[iseas, :])[0]),
        )
    valid_stations_indices = sorted(
        list(set.intersection(*valid_station_indices_for_seasons)),
    )
    valid_stations = obs_cube.coord(
        "Station index (arbitrary)",
    ).points[valid_stations_indices]

    # Only keep stations for which these conditions are met
    station_indices = obs_cube.coord("Station index (arbitrary)").points
    mask = np.isin(station_indices, valid_stations)
    obs_cube = obs_cube.extract(
        iris.Constraint(
            coord_values={
                "Station index (arbitrary)": lambda cell: mask[
                    int(cell.point)
                ],
            },
        ),
    )

    # Realize the observational data to avoid reading it from disk multiple times.
    obs_cube.data  # noqa B018
    multi_annual_seasonal_mean.data  # noqa B018

    return obs_cube, multi_annual_seasonal_mean


def preprocess_colocated_datasets(
    config,
    obs_cube,
    obs_seasonal,
    mod_datasets,
    seasons,
):
    """Preprocess obs and model datasets to colocate.

    Parameters
    ----------
    config : dict
        ESMValTool recipe configuration.
    obs_cube : iris.cube.Cube
        Cube containing the observation data.
    obs_seasonal : iris.cube.Cube
        Cube containing the seasonal mean of the observation data.
    mod_datasets : dict
        Dictionary containing the model dataset entries from the recipe.
    seasons : list
        List of the seasons for which to compute the seasonal means.

    Returns
    -------
    preproc_datasets: dict
        Dictionary containing the preprocessed observation and model datasets.
        Entries are grouped by requested diagnostic:
            - "seas_maps" for seasonal maps
            - "latitude_slices" for timeseries per latitude slice
    """
    preproc_datasets = {}

    trace_gas = config["trace_gas"]

    # General parameters
    months = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    # Load model data cubes for later iterations w/ quick fix
    cube_model_datasets = {}
    for model_dataset, group in mod_datasets.items():
        cube_model_datasets[model_dataset] = {
            attr["alias"]: _quick_fix_cube(attr["filename"], trace_gas)
            for attr in group
        }

    # Preprocessing for: seasonal maps and scatter plot
    #   plots parameter = "seas_maps"
    if "seas_maps" in config["plots"]:
        # Setup dictionary
        preproc_datasets["seas_maps"] = {}

        # Seasonal obs dataset
        preproc_datasets["seas_maps"]["obs"] = obs_seasonal
        noaa_gml_lats_seas = obs_seasonal.coord("latitude").points.tolist()
        noaa_gml_lons_seas = obs_seasonal.coord("longitude").points.tolist()

        # Looping over model datasets
        for model_dataset, group in mod_datasets.items():
            preproc_datasets["seas_maps"][model_dataset] = {}
            for attr in group:
                # Process multi-annual seasonal mean
                cube_seasonal = climate_statistics(
                    cube_model_datasets[model_dataset][attr["alias"]]["cube"],
                    operator="mean",
                    period="season",
                    seasons=seasons,
                )
                # Save outputs to dictionary
                preproc_datasets["seas_maps"][model_dataset][attr["alias"]] = {
                    "values": cube_seasonal,
                    "colocated": _extract_pt(
                        cube_seasonal,
                        noaa_gml_lats_seas,
                        noaa_gml_lons_seas,
                        nearest=True,
                    ),
                }

    # Preprocessing for:
    #   - colocated monthly time series ("taylor_diag")
    #   - colocated timeseries/scatter per latitude slices
    #     ("timeserie_lat" or "sensitivity_ampl_growth")
    plots_coloc = [
        p for p in config["plots"] if p not in ("seas_maps", "taylor_diag")
    ]
    plots_coloc_check = ["timeserie_lat", "sensitivity_ampl_growth"]
    if any(map(lambda v: v in plots_coloc_check, plots_coloc)):
        years_obs = obs_cube.coord("year", dim_coords=False).points.tolist()

        # Setup base entry for latitude slices
        preproc_datasets["latitude_slices"] = {}

        # Setup dictionary to save output for plots w/ the latitude time series
        # or the amplitude/sensitivity time series
        preproc_datasets["latitude_slices"] = {
            key: {model_dataset: {} for model_dataset in mod_datasets.keys()}
            for key, _ in LATITUDE_RANGES.items()
        }
        for key, _ in LATITUDE_RANGES.items():
            preproc_datasets["latitude_slices"][key]["obs"] = {
                "max": {str(y): {m: [] for m in months} for y in years_obs},
                "min": {str(y): {m: [] for m in months} for y in years_obs},
                "yearly": {
                    "mean": {str(y): 0.0 for y in years_obs},
                    "std": {str(y): 0.0 for y in years_obs},
                    "valid": {str(y): [] for y in years_obs},
                },
                "anomaly": {
                    str(y): {m: [] for m in months} for y in years_obs
                },
                "amplitude": {
                    "mean": {str(y): 0.0 for y in years_obs},
                    "std": {str(y): 0.0 for y in years_obs},
                    "valid": {str(y): [] for y in years_obs},
                },
                "growth": {
                    "mean": {str(y): 0.0 for y in years_obs[:1]},
                    "std": {str(y): 0.0 for y in years_obs[:1]},
                    "valid": {str(y): [] for y in years_obs[:1]},
                },
                "sensitivity": {
                    "mean": {str(y): 0.0 for y in years_obs[:1]},
                    "std": {str(y): 0.0 for y in years_obs[:1]},
                    "valid": {str(y): [] for y in years_obs[:1]},
                },
            }
        # Setup yearly mean, amplitude, and growth cubes for later
        # colocation w/ latitude slices
        model_slices = {key: {} for key, _ in LATITUDE_RANGES.items()}
        obs_slices = {key: {} for key, _ in LATITUDE_RANGES.items()}
        for key, (lat_min, lat_max) in LATITUDE_RANGES.items():
            constraint = iris.Constraint(
                latitude=lambda cell, lat_min=lat_min, lat_max=lat_max: lat_min
                <= cell
                < lat_max,
            )
            # Extract latitude slices for observations
            obs_slices[key]["values"] = obs_cube.extract(constraint)
            obs_slices[key]["yearly"] = obs_slices[key][
                "values"
            ].aggregated_by("year", iris.analysis.MEAN)
            obs_slices[key]["amplitude"] = obs_cube.extract(
                constraint,
            ).aggregated_by("year", iris.analysis.MAX) - obs_cube.extract(
                constraint,
            ).aggregated_by("year", iris.analysis.MIN)
            obs_slices[key]["growth"] = _setup_growth_cube(
                obs_slices[key]["yearly"],
                config,
                obs_slices[key]["values"],
                "obs",
            )
            # Extract latitude slices for model data
            model_slices[key] = {
                "values": {
                    model_dataset: {} for model_dataset in mod_datasets.keys()
                },
                "yearly": {
                    model_dataset: {} for model_dataset in mod_datasets.keys()
                },
                "amplitude": {
                    model_dataset: {} for model_dataset in mod_datasets.keys()
                },
                "growth": {
                    model_dataset: {} for model_dataset in mod_datasets.keys()
                },
            }
            for model_dataset, group in mod_datasets.items():
                for attr in group:
                    mod_tmp = cube_model_datasets[model_dataset][
                        attr["alias"]
                    ]["cube"].extract(constraint)
                    model_slices[key]["values"][model_dataset][
                        attr["alias"]
                    ] = mod_tmp
                    model_slices[key]["yearly"][model_dataset][
                        attr["alias"]
                    ] = mod_tmp.aggregated_by("year", iris.analysis.MEAN)
                    model_slices[key]["amplitude"][model_dataset][
                        attr["alias"]
                    ] = mod_tmp.aggregated_by(
                        "year",
                        iris.analysis.MAX,
                    ) - mod_tmp.aggregated_by("year", iris.analysis.MIN)
                    model_slices[key]["growth"][model_dataset][
                        attr["alias"]
                    ] = _setup_growth_cube(
                        model_slices[key]["yearly"][model_dataset][
                            attr["alias"]
                        ],
                        config,
                        mod_tmp,
                        "model",
                    )

        # Looping over model datasets for the colocation process
        for model_dataset, group in mod_datasets.items():
            for key, _ in LATITUDE_RANGES.items():
                preproc_datasets["latitude_slices"][key][model_dataset] = {}
            # Looping over variables in dataset group
            for attr in group:
                # Select time, cube, and time steps
                time = cube_model_datasets[model_dataset][attr["alias"]][
                    "time"
                ]
                cube = cube_model_datasets[model_dataset][attr["alias"]][
                    "cube"
                ]
                years = cube_model_datasets[model_dataset][attr["alias"]][
                    "years"
                ]
                months_ts = cube.coord("month", dim_coords=False).points
                years_ts = cube.coord("year", dim_coords=False).points
                for key, _ in LATITUDE_RANGES.items():
                    preproc_datasets["latitude_slices"][key][model_dataset][
                        attr["alias"]
                    ] = {
                        "max": {
                            str(y): {m: [] for m in months} for y in years
                        },
                        "min": {
                            str(y): {m: [] for m in months} for y in years
                        },
                        "yearly": {
                            "mean": {str(y): 0.0 for y in years},
                            "std": {str(y): 0.0 for y in years},
                            "valid": {str(y): [] for y in years},
                        },
                        "anomaly": {
                            str(y): {m: [] for m in months} for y in years
                        },
                        "amplitude": {
                            "mean": {str(y): 0.0 for y in years},
                            "std": {str(y): 0.0 for y in years},
                            "valid": {str(y): [] for y in years},
                        },
                        "growth": {
                            "mean": {str(y): 0.0 for y in years[:1]},
                            "std": {str(y): 0.0 for y in years[:1]},
                            "valid": {str(y): [] for y in years[:1]},
                        },
                        "sensitivity": {
                            "mean": {str(y): 0.0 for y in years[:1]},
                            "std": {str(y): 0.0 for y in years[:1]},
                            "valid": {str(y): [] for y in years[:1]},
                        },
                    }

                # Looping over all time steps
                for i, _ts in enumerate(time.units.num2date(time.points)):
                    if "timeserie_lat" in plots_coloc:
                        # Looping over latitude slices
                        for key, _ in LATITUDE_RANGES.items():
                            # Obs values
                            obs_ts = obs_slices[key]["values"][i]
                            obs_yearly = obs_slices[key]["yearly"].extract(
                                iris.Constraint(year=years_ts[i]),
                            )
                            # Model values
                            mod_ts = model_slices[key]["values"][
                                model_dataset
                            ][attr["alias"]][i]
                            mod_yearly = model_slices[key]["yearly"][
                                model_dataset
                            ][attr["alias"]].extract(
                                iris.Constraint(year=years_ts[i]),
                            )

                            # Extract monthly min, max
                            v_obs, v_mod, _ = _colocate_obs_model(
                                obs_ts,
                                mod_ts,
                            )
                            # Save outputs to dictionary
                            preproc_datasets["latitude_slices"][key][
                                model_dataset
                            ][attr["alias"]]["max"][str(years_ts[i])][
                                months_ts[i]
                            ] = np.max(v_mod) if len(v_mod) > 0 else np.nan
                            preproc_datasets["latitude_slices"][key][
                                model_dataset
                            ][attr["alias"]]["min"][str(years_ts[i])][
                                months_ts[i]
                            ] = np.min(v_mod) if len(v_mod) > 0 else np.nan
                            preproc_datasets["latitude_slices"][key]["obs"][
                                "max"
                            ][str(years_ts[i])][months_ts[i]] = (
                                np.max(v_obs) if len(v_obs) > 0 else np.nan
                            )
                            preproc_datasets["latitude_slices"][key]["obs"][
                                "min"
                            ][str(years_ts[i])][months_ts[i]] = (
                                np.min(v_obs) if len(v_obs) > 0 else np.nan
                            )

                            # Extract monthly anomaly
                            obs_anomaly = obs_ts - obs_yearly
                            mod_anomaly = mod_ts - mod_yearly
                            v_obs, v_mod, _ = _colocate_obs_model(
                                obs_anomaly,
                                mod_anomaly,
                            )
                            # Save outputs to dictionary
                            preproc_datasets["latitude_slices"][key][
                                model_dataset
                            ][attr["alias"]]["anomaly"][str(years_ts[i])][
                                months_ts[i]
                            ] = np.array(v_mod)
                            preproc_datasets["latitude_slices"][key]["obs"][
                                "anomaly"
                            ][str(years_ts[i])][months_ts[i]] = np.array(v_obs)

                # Looping over years
                for i, y in enumerate(years):
                    # Looping over latitude slices
                    for key, _ in LATITUDE_RANGES.items():
                        # Obs values
                        obs_y = obs_slices[key]["yearly"].extract(
                            iris.Constraint(year=y),
                        )
                        obs_a = obs_slices[key]["amplitude"].extract(
                            iris.Constraint(year=y),
                        )
                        # Model values
                        mod_y = model_slices[key]["yearly"][model_dataset][
                            attr["alias"]
                        ].extract(iris.Constraint(year=y))
                        mod_a = model_slices[key]["amplitude"][model_dataset][
                            attr["alias"]
                        ].extract(iris.Constraint(year=y))

                        # Extract yearly
                        v_obs_y, v_mod_y, _ = _colocate_obs_model(obs_y, mod_y)
                        # Save outputs to dictionary
                        preproc_datasets["latitude_slices"][key]["obs"][
                            "yearly"
                        ]["mean"][str(y)] = np.mean(v_obs_y)
                        preproc_datasets["latitude_slices"][key]["obs"][
                            "yearly"
                        ]["std"][str(y)] = np.std(v_obs_y)
                        preproc_datasets["latitude_slices"][key]["obs"][
                            "yearly"
                        ]["valid"][str(y)] = v_obs_y
                        preproc_datasets["latitude_slices"][key][
                            model_dataset
                        ][attr["alias"]]["yearly"]["mean"][str(y)] = np.mean(
                            v_mod_y,
                        )
                        preproc_datasets["latitude_slices"][key][
                            model_dataset
                        ][attr["alias"]]["yearly"]["std"][str(y)] = np.std(
                            v_mod_y,
                        )
                        preproc_datasets["latitude_slices"][key][
                            model_dataset
                        ][attr["alias"]]["yearly"]["valid"][str(y)] = v_mod_y

                        if "sensitivity_ampl_growth" in plots_coloc:
                            # Extract amplitude
                            v_obs_a, v_mod_a, _ = _colocate_obs_model(
                                obs_a,
                                mod_a,
                            )
                            # Save outputs to dictionary
                            preproc_datasets["latitude_slices"][key]["obs"][
                                "amplitude"
                            ]["mean"][str(y)] = np.mean(v_obs_a)
                            preproc_datasets["latitude_slices"][key]["obs"][
                                "amplitude"
                            ]["std"][str(y)] = np.std(v_obs_a)
                            preproc_datasets["latitude_slices"][key]["obs"][
                                "amplitude"
                            ]["valid"][str(y)] = v_obs_a
                            preproc_datasets["latitude_slices"][key][
                                model_dataset
                            ][attr["alias"]]["amplitude"]["mean"][
                                str(y)
                            ] = np.mean(v_mod_a)
                            preproc_datasets["latitude_slices"][key][
                                model_dataset
                            ][attr["alias"]]["amplitude"]["std"][
                                str(y)
                            ] = np.std(v_mod_a)
                            preproc_datasets["latitude_slices"][key][
                                model_dataset
                            ][attr["alias"]]["amplitude"]["valid"][
                                str(y)
                            ] = v_mod_a

                            # Extract growth & sensitivity = skip if first year as no data
                            if i > 0:
                                # Obs values
                                obs_g = obs_slices[key]["growth"].extract(
                                    iris.Constraint(year=y),
                                )
                                # Model values
                                mod_g = model_slices[key]["growth"][
                                    model_dataset
                                ][attr["alias"]].extract(
                                    iris.Constraint(year=y),
                                )

                                # Extract growth
                                v_obs_g, v_mod_g, _ = _colocate_obs_model(
                                    obs_g,
                                    mod_g,
                                )
                                # Compute relative growth
                                v_obs_g_relative = [
                                    100 * g / v_obs_y[j - 1]
                                    for j, g in enumerate(v_obs_g)
                                ]
                                v_mod_g_relative = [
                                    100 * g / v_mod_y[j - 1]
                                    for j, g in enumerate(v_mod_g)
                                ]
                                # Save outputs to dictionary
                                preproc_datasets["latitude_slices"][key][
                                    "obs"
                                ]["growth"]["mean"][str(y)] = np.mean(
                                    v_obs_g_relative,
                                )
                                preproc_datasets["latitude_slices"][key][
                                    "obs"
                                ]["growth"]["std"][str(y)] = np.std(
                                    v_obs_g_relative,
                                )
                                preproc_datasets["latitude_slices"][key][
                                    "obs"
                                ]["growth"]["valid"][str(y)] = v_obs_g_relative
                                preproc_datasets["latitude_slices"][key][
                                    model_dataset
                                ][attr["alias"]]["growth"]["mean"][
                                    str(y)
                                ] = np.mean(v_mod_g_relative)
                                preproc_datasets["latitude_slices"][key][
                                    model_dataset
                                ][attr["alias"]]["growth"]["std"][
                                    str(y)
                                ] = np.std(v_mod_g_relative)
                                preproc_datasets["latitude_slices"][key][
                                    model_dataset
                                ][attr["alias"]]["growth"]["valid"][
                                    str(y)
                                ] = v_mod_g_relative

                                # Compute sensitivity between amplitude and growth
                                v_obs_s = [
                                    v_obs_a[j] / v_obs_g[j]
                                    for j in range(1, len(v_obs_g))
                                ]
                                v_mod_s = [
                                    v_mod_a[j] / v_mod_g[j]
                                    for j in range(1, len(v_mod_g))
                                ]
                                # Save outputs to dictionary
                                preproc_datasets["latitude_slices"][key][
                                    "obs"
                                ]["sensitivity"]["mean"][str(y)] = np.mean(
                                    v_obs_s,
                                )
                                preproc_datasets["latitude_slices"][key][
                                    "obs"
                                ]["sensitivity"]["std"][str(y)] = np.std(
                                    v_obs_s,
                                )
                                preproc_datasets["latitude_slices"][key][
                                    "obs"
                                ]["sensitivity"]["valid"][str(y)] = v_obs_s
                                preproc_datasets["latitude_slices"][key][
                                    model_dataset
                                ][attr["alias"]]["sensitivity"]["mean"][
                                    str(y)
                                ] = np.mean(v_mod_s)
                                preproc_datasets["latitude_slices"][key][
                                    model_dataset
                                ][attr["alias"]]["sensitivity"]["std"][
                                    str(y)
                                ] = np.std(v_mod_s)
                                preproc_datasets["latitude_slices"][key][
                                    model_dataset
                                ][attr["alias"]]["sensitivity"]["valid"][
                                    str(y)
                                ] = v_mod_s

    return preproc_datasets


def main(config):
    """Produce surface trace gas climatology diagnostic.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    # Produce climatology for observational dataset
    seasons = ["DJF", "MAM", "JJA", "SON"]
    start_year, end_year = (
        datasets[config["observational_dataset"]][0]["start_year"],
        datasets[config["observational_dataset"]][0]["end_year"],
    )
    obs_dataset = datasets.pop(config["observational_dataset"])
    obs_cube, obs_cube_seasonal = preprocess_obs_dataset(obs_dataset, config)

    # Default value for the parameter plots:
    #   - plots = "seas_maps", "timeserie_lat",
    #       "sensitivity_ampl_growth", "taylor_diag"
    if "plots" not in config:
        config["plots"] = [
            "seas_maps",
            "timeserie_lat",
            "sensitivity_ampl_growth",
            "taylor_diag",
        ]

    # Colocating model and obs datasets
    logger.info("Pre-processing datasets: colocation w/ surface obs")
    colocated_datasets = preprocess_colocated_datasets(
        config,
        obs_cube,
        obs_cube_seasonal,
        datasets,
        seasons,
    )

    for model_dataset, group in datasets.items():
        # "model_dataset" is the name of the model dataset.
        # "group" is a list of dictionaries containing metadata.
        logger.info("Plotting for dataset %s", model_dataset)

        for attributes in group:
            # Set up for analysis and plotting
            plot_file_prefix = (
                model_dataset
                + "_"
                + attributes["mip"]
                + "_"
                + attributes["exp"]
                + "_"
                + attributes["short_name"]
                + "_"
                + str(attributes["start_year"])
                + "_"
                + str(attributes["end_year"])
                + "_"
            )

            if "seas_maps" in config["plots"]:
                logger.info("Producing seasonal maps...")
                # Analysis and plotting for model-obs seasonal map comparison
                figures_seas_map, fig_scatter_seas = trace_gas_maps(
                    model_dataset,
                    attributes,
                    colocated_datasets,
                    seasons,
                    attributes["timerange"],
                    config["trace_gas"],
                )
                # Saving plots
                output_file = plot_file_prefix + "scatter"
                output_path = get_plot_filename(output_file, config)
                fig_scatter_seas.savefig(output_path)
                provenance = get_provenance_record(
                    ancestors=[
                        attributes["filename"],
                        obs_dataset[0]["filename"],
                    ],
                    model_name=attributes["dataset"],
                    project=attributes["project"],
                    experiment=attributes["exp"],
                    timerange=attributes["timerange"],
                    trace_gas=config["trace_gas"],
                    obs_name=obs_dataset[0]["dataset"],
                    var="seas_maps_scatter",
                )
                with ProvenanceLogger(config) as provenance_logger:
                    provenance_logger.log(output_path, provenance)

                output_file = plot_file_prefix + "seas_map"
                output_path = get_plot_filename(output_file, config)
                figures_seas_map.savefig(output_path, bbox_inches="tight")
                provenance = get_provenance_record(
                    ancestors=[
                        attributes["filename"],
                        obs_dataset[0]["filename"],
                    ],
                    model_name=attributes["dataset"],
                    project=attributes["project"],
                    experiment=attributes["exp"],
                    timerange=attributes["timerange"],
                    trace_gas=config["trace_gas"],
                    obs_name=obs_dataset[0]["dataset"],
                    var="seas_maps_map",
                )
                with ProvenanceLogger(config) as provenance_logger:
                    provenance_logger.log(output_path, provenance)

            if "timeserie_lat" in config["plots"]:
                logger.info(
                    "Producing time series of trend and seasonal cycle for latitude slices...",
                )
                # Analysis and plotting for zonal-mean time series for model-obs
                figure_timeserie_zonal = trace_gas_timeserie_zonal(
                    model_dataset,
                    attributes,
                    colocated_datasets,
                    attributes["timerange"],
                    config["trace_gas"],
                )
                # Save time series plots per latitude range
                output_file = plot_file_prefix + "timeseries_latitude"
                output_path = get_plot_filename(output_file, config)
                figure_timeserie_zonal.savefig(
                    output_path,
                    bbox_inches="tight",
                )
                provenance = get_provenance_record(
                    ancestors=[
                        attributes["filename"],
                        obs_dataset[0]["filename"],
                    ],
                    model_name=attributes["dataset"],
                    project=attributes["project"],
                    experiment=attributes["exp"],
                    timerange=attributes["timerange"],
                    trace_gas=config["trace_gas"],
                    obs_name=obs_dataset[0]["dataset"],
                    var="timeserie_lat",
                )
                with ProvenanceLogger(config) as provenance_logger:
                    provenance_logger.log(output_path, provenance)

            if "sensitivity_ampl_growth" in config["plots"]:
                logger.info(
                    "Producing time series of amplitude, growth, and sensitivity for latitude slices...",
                )
                # Analysis and plotting for sensitivity between
                # seasonal cycle amplitude and growth rate
                figure_sens_ampl_trend = trace_gas_seas_ampl_growth_rate(
                    model_dataset,
                    attributes,
                    colocated_datasets,
                    attributes["timerange"],
                    config["trace_gas"],
                )
                # Save sensitivity of amplitude and growth trend plot
                output_file = plot_file_prefix + "sensitivity_ampl_growth"
                output_path = get_plot_filename(output_file, config)
                figure_sens_ampl_trend.savefig(
                    output_path,
                    bbox_inches="tight",
                )
                provenance = get_provenance_record(
                    ancestors=[
                        attributes["filename"],
                        obs_dataset[0]["filename"],
                    ],
                    model_name=attributes["dataset"],
                    project=attributes["project"],
                    experiment=attributes["exp"],
                    timerange=attributes["timerange"],
                    trace_gas=config["trace_gas"],
                    obs_name=obs_dataset[0]["dataset"],
                    var="sensitivity_ampl_growth",
                )
                with ProvenanceLogger(config) as provenance_logger:
                    provenance_logger.log(output_path, provenance)

    if "taylor_diag" in config["plots"]:
        logger.info("Producing Taylor diagram...")
        # Analysis and plotting for Taylor diagram
        figure_taylor_diag = trace_gas_taylor_diag(
            datasets,
            obs_cube,
            config["trace_gas"],
            config,
        )
        # Save Taylor diagram plot
        model_names = [
            att["dataset"] for gr in datasets.values() for att in gr
        ]
        output_file = (
            "trace_gas_"
            + config["trace_gas"]
            + "_"
            + "_".join(model_names)
            + "_"
            + str(start_year)
            + "_"
            + str(end_year)
            + "_taylor_diag"
        )
        output_path = get_plot_filename(output_file, config)
        figure_taylor_diag.savefig(output_path, bbox_inches="tight", dpi=300)
        provenance = get_provenance_record(
            ancestors=[
                att["filename"] for gr in datasets.values() for att in gr
            ]
            + [obs_dataset[0]["filename"]],
            model_name=model_names,
            project=list(
                {att["project"] for gr in datasets.values() for att in gr},
            ),
            experiment=list(
                {att["exp"] for gr in datasets.values() for att in gr},
            ),
            timerange=obs_dataset[0]["timerange"],
            trace_gas=config["trace_gas"],
            obs_name=obs_dataset[0]["dataset"],
            var="taylor_diag",
        )
        with ProvenanceLogger(config) as provenance_logger:
            provenance_logger.log(output_path, provenance)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
