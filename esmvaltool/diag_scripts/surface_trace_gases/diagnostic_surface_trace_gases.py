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
import iris.plot as iplt
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy
from esmvalcore.preprocessor._multimodel import multi_model_statistics
from esmvalcore.preprocessor._time import climate_statistics
from matplotlib import colors, gridspec
from numpy import ma

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from esmvaltool.diag_scripts.surface_trace_gases.utils_surface_trace_gases import (
    add_bounds,
    extract_pt,
)

logger = logging.getLogger(Path(__file__).stem)
fontsizedict = {"title": 25, "axis": 20, "legend": 18, "ticklabel": 18}
TRACE_GASES_FACTOR = {
    "ch4": 1e9,
    "co2": 1e6,
    "n2o": 1e9,
}
TRACE_GASES_UNITS = {
    "ch4": "ppb",
    "co2": "ppm",
    "n2o": "ppb",
}

# Markers and corresponding colors to be used for plots
COLORS_MARKERS = ["cornflowerblue", "royalblue", "lightsteelblue"]
MARKERS = ["o", "^", "s"]


def weighted_std_dev(cube, dim, weights=None):
    """Compute a weighted standard deviation.

    Code is taken from https://github.com/SciTools/iris/issues/3208.

    Parameters
    ----------
    cube : iris.cube.Cube
        The cube containing the input data.
    dim : list
        List containing one or more dimensions of the cube over which to
        compute the standard deviation.
    weights : optional
        None or output of iris.analysis.cartography.area_weights.

    Returns
    -------
    Collapsed weighted standard deviation cube over the dimension(s) dim.
    """
    if weights is None:
        return (
            cube.collapsed(dim, iris.analysis.RMS) ** 2
            - cube.collapsed(dim, iris.analysis.MEAN) ** 2
        ) ** 0.5
    return (
        cube.collapsed(dim, iris.analysis.RMS, weights=weights) ** 2
        - cube.collapsed(dim, iris.analysis.MEAN, weights=weights) ** 2
    ) ** 0.5


def plot_trace_gas_mod_obs(
    fig, ax, md_data, obs_data, trace_gas_obs_cube, plot_dict
):
    """Plot trace gas surface concentration.

    The function plots a contour for the model data overlaid
    with NOAA GML surface flask observations climatology.

    Parameters
    ----------
    fig: matplotlib figure
        Figure object.
    ax: matplotlib axes
        Axes of the figure.
    md_data : Iris cube
        Model trace gas surface concentration as a cube with latitude
        and longitude coordinates.
    obs_data : List.
        Observations of trace gas surface concentration from each
        NOAA GML station.
    trace_gas_obs_cube : Iris cube.
        Holds information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    plot_dict : Dictionary.
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
    colbar.set_ticks(plot_dict["Levels"])
    colbar.set_ticklabels(plot_dict["tick_labels"], fontsize=14)
    colbar.set_label(plot_dict["cb_label"], fontsize=14)
    ax.coastlines(color="#525252")

    # Statistics on plot
    ax.text(
        0.5,
        -0.1,
        (
            f"Global mean {plot_dict['Mean']:.1f} - @Stations mean: mod="
            f"{plot_dict['Stn_mn_md']:.1f} & obs={plot_dict['Stn_mn_obs']:.1f}"
            f" (RMSE={plot_dict['RMS']:.1f})"
        ),
        ha="center",
        va="center",
        size=14,
        transform=ax.transAxes,
    )


def trace_gas_maps(
    model_data, trace_gas_obs_cube, clim_seas, timerange, trace_gas
):
    """Plot model vs NOAA GML trace gas surface concentration.

    The function generates the plots and returns evaluation metrics.

    Parameters
    ----------
    model_data : Iris Cube.
        Contains model output of surface concentration with coordinates;
        time, latitude and longitude.
    trace_gas_obs_cube : Iris Cube.
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    clim_seas : List.
       Strings to denote climate seasons ["DJF", "MAM", "JJA", "SON"]
    timerange : String.
        String containing the time range as "{year_start}/{year_end}".
    trace_gas : String.
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.

    Returns
    -------
    figures : List.
        Contains figure instances for the seasonal contour plots overlaid with
        observations of surface concentration from NOAA GML.
    fig_scatter : Figure object.
        The scatter plot comparing model and obs surface concentrations.
    """
    # Get model run id
    if "parent_source_id" in model_data.attributes:
        model_id = model_data.attributes["parent_source_id"]
    else:
        model_id = "Multi-Model-Mean"

    # Co-locate model grid points with measurement sites
    noaa_gml_lats = trace_gas_obs_cube.coord("latitude").points.tolist()
    noaa_gml_lons = trace_gas_obs_cube.coord("longitude").points.tolist()
    trace_gas_at_noaa = extract_pt(
        model_data, noaa_gml_lats, noaa_gml_lons, nearest=True
    )

    # Set up seasonal contour plots
    figures = []

    # Use global average as center value for the colormap
    grid_areas_weights = iris.analysis.cartography.area_weights(model_data)
    center_cmap_mean = np.round(
        model_data.collapsed(
            ["season_number", "latitude", "longitude"],
            iris.analysis.MEAN,
            weights=grid_areas_weights,
        ).data,
    )
    step = 5
    clevs = list(
        np.arange(
            center_cmap_mean - center_cmap_mean % 5 - 15,
            center_cmap_mean - center_cmap_mean % 5 + 16,
            step,
            dtype=int,
        )
    )
    clabs = [str(lev) if (lev % 5 == 0) else "" for lev in clevs]
    cmapr = mpl.colormaps.get_cmap("Blues")
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
    col_scatter = cmapr(np.linspace(0.25, 1, 4))
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
    for s, season in enumerate(trace_gas_obs_cube.slices_over("clim_season")):
        # Match NOAA GML obs season with model season number
        model_sn = [c.lower() for c in clim_seas].index(
            season.coord("clim_season").points[0]
        )
        model_season = model_data[model_sn]

        logger.info(
            "Analysing %s for %s: %s",
            trace_gas.upper(),
            model_id,
            clim_seas[model_sn],
        )

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
            )
        )

        # Plot contours overlaid with obs for this run and season
        n_stn = str(len(valid_obs))
        title = (
            f"\nSurface {trace_gas.upper()} concentration "
            + timerange
            + "\n"
            + model_id
            + ", "
            + clim_seas[model_sn]
            + ", N stations="
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
            trace_gas_obs_cube,
            plot_dict,
        )

    figures = fig_cf

    # Decorate the scatter plot
    ax_scatter.set(
        xlim=(min_scatter - 2, max_scatter + 2),
        xticks=np.arange(
            min_scatter - min_scatter % 10,
            max_scatter + min_scatter % 10,
            step=10,
        ),
        ylim=(min_scatter - 2, max_scatter + 2),
        yticks=np.arange(
            min_scatter - min_scatter % 10,
            max_scatter + min_scatter % 10,
            step=10,
        ),
    )
    ax_scatter.set_xlabel(
        f"NOAA GML Flask {trace_gas.upper()} "
        f"[{model_data.attributes['unit']}]",
        fontsize=fontsizedict["axis"],
    )
    ax_scatter.set_ylabel(
        f"{model_id} {trace_gas.upper()} [{model_data.attributes['unit']}]",
        fontsize=fontsizedict["axis"],
    )

    ax_scatter.tick_params(
        axis="both", which="major", labelsize=fontsizedict["ticklabel"], pad=10
    )

    ax_scatter.set_title(
        "Model vs observations @stations\n"
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
    model_data,
    obs_cube,
    trace_gas,
    include_global=False,
):
    """Plot time series of zonal mean trace gas concentration.

    It uses different latitude slices for model and observational data:
    - 90S - 30S
    - 30S - 30N
    - 30N - 60N
    - 60N - 90N

    Parameters
    ----------
    model_data : Iris Cube.
        Contains model output of surface concentration with coordinates;
         time, latitude and longitude.
    obs_cube : Iris Cube.
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    trace_gas : String.
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.
    include_global : Bool
        Boolean flag to indicate if the global values should be overlaid.

    Returns
    -------
    figure : Figure object.
        Figure for the zonal time series of the trace gas for model data
        and observations of surface concentration from NOAA GML.
        The left column of the figure contains the zonal time series.
        The center-left column contains the zonal scatter plots for model-obs.
        The center-right column contains the zonal seasonal anomaly.
        The right column contains the zonal max-min months.
    """
    # Get model run id
    if "parent_source_id" in model_data.attributes:
        model_id = model_data.attributes["parent_source_id"]
    else:
        model_id = "Multi-Model-Mean"

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
    latitude_titles = [
        r"Latitudes 60$^\circ$N - 90$^\circ$N",
        r"Latitudes 30$^\circ$N - 60$^\circ$N",
        r"Latitudes 30$^\circ$S - 30$^\circ$N",
        r"Latitudes 90$^\circ$S - 30$^\circ$S",
    ]
    handles_mean = None
    labels_mean = None
    handles_minmax = None
    labels_minmax = None

    # Define latitude ranges to slice the model and obs cubes
    latitude_ranges = {
        "60N - 90N": (60, 90),
        "30N - 60N": (30, 60),
        "30S - 30N": (-30, 30),
        "90S - 30S": (-90, -30),
    }
    model_slices = {}
    obs_slices = {}
    for key, (lat_min, lat_max) in latitude_ranges.items():
        constraint = iris.Constraint(
            latitude=lambda cell, lat_min=lat_min, lat_max=lat_max: lat_min
            <= cell
            < lat_max,
        )
        model_slices[key] = model_data.extract(constraint)
        obs_slices[key] = obs_cube.extract(constraint)

    for l_i, lat_range in enumerate(latitude_ranges.keys()):
        # Select relevant latitude range
        obs = obs_slices[lat_range]
        mod = model_slices[lat_range]

        # Add year/month in the time coordinate of the model data
        time = mod.coord("time", dim_coords=True)
        try:
            iris.coord_categorisation.add_year(mod, "time")
            iris.coord_categorisation.add_month(mod, "time")
        except ValueError:
            err = f"Year/month coordinates already present in {model_id}."
            logger.debug(err)
        years = list(set(mod.coord("year", dim_coords=False).points))
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

        # Co-locate model grid points with measurement sites
        noaa_gml_lats = obs.coord("latitude").points.tolist()
        noaa_gml_lons = obs.coord("longitude").points.tolist()
        # Yearly time series
        trace_gas_at_model_y = {}
        trace_gas_at_station_y = {}
        model_yearly = mod.aggregated_by("year", iris.analysis.MEAN)
        obs_yearly = obs.aggregated_by("year", iris.analysis.MEAN)
        valid_obs = []
        valid_md = []
        for y in years:
            obs_y = obs_yearly.extract(iris.Constraint(year=y))
            extracted_gas = extract_pt(
                model_yearly.extract(iris.Constraint(year=y)),
                noaa_gml_lats,
                noaa_gml_lons,
                nearest=True,
            )
            # Select only stations with valid data
            valid_indices = ~(obs_y.data.mask | np.isnan(extracted_gas))
            v_obs = obs_y.data[valid_indices].tolist()
            v_mod = [
                tg.item()
                for i, tg in enumerate(extracted_gas)
                if (tg.item() is not None) and valid_indices[i]
            ]
            # Save values in bulk lists for scatter plot
            valid_obs += v_obs
            valid_md += v_mod
            # Save values in dictionary for time series plot
            trace_gas_at_model_y[str(y)] = np.array(v_mod)
            trace_gas_at_station_y[str(y)] = np.array(v_obs)
        # Time series of monthly data and seasonal anomalies
        trace_gas_at_model = {
            "max": {str(y): {m: [] for m in months} for y in years},
            "min": {str(y): {m: [] for m in months} for y in years},
        }
        trace_gas_at_station = {
            "max": {str(y): {m: [] for m in months} for y in years},
            "min": {str(y): {m: [] for m in months} for y in years},
        }
        trace_gas_at_model_anom = {
            str(y): {m: [] for m in months} for y in years
        }
        trace_gas_at_station_anom = {
            str(y): {m: [] for m in months} for y in years
        }
        months_ts = mod.coord("month", dim_coords=False).points
        years_ts = mod.coord("year", dim_coords=False).points
        for i, ts in enumerate(time.units.num2date(time.points)):
            # Monthly max/min time series
            obs_ts = obs.extract(iris.Constraint(time=ts))
            extracted_gas = extract_pt(
                mod.extract(iris.Constraint(time=ts)),
                noaa_gml_lats,
                noaa_gml_lons,
                nearest=True,
            )
            valid_indices = ~(obs_ts.data.mask | np.isnan(extracted_gas))
            v_obs = obs_ts.data[valid_indices].tolist()
            v_mod = [
                tg.item()
                for i, tg in enumerate(extracted_gas)
                if (tg.item() is not None) and valid_indices[i]
            ]
            trace_gas_at_model["max"][str(years_ts[i])][months_ts[i]] = np.max(
                v_mod
            )
            trace_gas_at_station["max"][str(years_ts[i])][months_ts[i]] = (
                np.max(v_obs)
            )
            trace_gas_at_model["min"][str(years_ts[i])][months_ts[i]] = np.min(
                v_mod
            )
            trace_gas_at_station["min"][str(years_ts[i])][months_ts[i]] = (
                np.min(v_obs)
            )
            # Monthly anomaly time series
            obs_anomaly = obs.extract(
                iris.Constraint(time=ts)
            ) - obs_yearly.extract(iris.Constraint(year=years_ts[i]))
            mod_anomaly = mod.extract(
                iris.Constraint(time=ts)
            ) - model_yearly.extract(iris.Constraint(year=years_ts[i]))
            extracted_anom = extract_pt(
                mod_anomaly,
                noaa_gml_lats,
                noaa_gml_lons,
                nearest=True,
            )
            valid_indices = ~(obs_anomaly.data.mask | np.isnan(extracted_anom))
            v_obs = obs_anomaly.data[valid_indices].tolist()
            v_mod = [
                tg.item()
                for i, tg in enumerate(extracted_anom)
                if (tg.item() is not None) and valid_indices[i]
            ]
            trace_gas_at_model_anom[str(years_ts[i])][months_ts[i]] = np.array(
                v_mod
            )
            trace_gas_at_station_anom[str(years_ts[i])][months_ts[i]] = (
                np.array(v_obs)
            )
        # Compute aggregated values for the entire latitude range
        if include_global:
            mod_latitude = {
                "mean": {
                    "yearly": np.zeros(shape=len(years)),
                    "monthly": np.zeros(shape=len(years) - 1),
                },
                "std": {
                    "yearly": np.zeros(shape=len(years)),
                    "monthly": np.zeros(shape=len(years) - 1),
                },
            }
            weights_yearly = iris.analysis.cartography.area_weights(
                model_yearly,
            )
            mod_latitude["mean"]["yearly"] = model_yearly.collapsed(
                ["latitude", "longitude"],
                iris.analysis.MEAN,
                weights=weights_yearly,
            ).data
            mod_latitude["std"]["yearly"] = weighted_std_dev(
                cube=model_yearly,
                dim=["latitude", "longitude"],
                weights=weights_yearly,
            ).data
            anomaly = mod.copy(data=np.full_like(mod.data, np.nan))
            for i, ts in enumerate(time.units.num2date(time.points)):
                anomaly.data[i, ...] = (
                    mod.extract(iris.Constraint(time=ts)).data
                    - model_yearly.extract(
                        iris.Constraint(year=years_ts[i])
                    ).data
                )
            anomaly = anomaly.aggregated_by("month", iris.analysis.MEAN)
            weights_anom = iris.analysis.cartography.area_weights(anomaly)
            mod_latitude["mean"]["monthly"] = (
                anomaly.collapsed(
                    ["latitude", "longitude"],
                    iris.analysis.MEAN,
                    weights=weights_anom,
                )
                .aggregated_by("month", iris.analysis.MEAN)
                .data
            )
            mod_latitude["std"]["monthly"] = (
                anomaly.collapsed(
                    ["latitude", "longitude"],
                    iris.analysis.MEAN,
                    weights=weights_anom,
                )
                .aggregated_by("month", iris.analysis.STD_DEV)
                .data
            )

        # Plot left column
        # zonal-mean time series and corresponding +/- std range
        model_ts_mean = np.array(
            [np.mean(trace_gas_at_model_y[str(y)]) for y in years]
        )
        model_ts_std = np.array(
            [np.std(trace_gas_at_model_y[str(y)]) for y in years]
        )
        obs_ts_mean = np.array(
            [np.mean(trace_gas_at_station_y[str(y)]) for y in years]
        )
        obs_ts_std = np.array(
            [np.std(trace_gas_at_station_y[str(y)]) for y in years]
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
        if include_global:
            plt.fill_between(
                years,
                mod_latitude["mean"]["yearly"] - mod_latitude["std"]["yearly"],
                mod_latitude["mean"]["yearly"] + mod_latitude["std"]["yearly"],
                color=COLORS_MARKERS[2],
                alpha=0.1,
            )
        plt.plot(
            years,
            model_ts_mean,
            linestyle="solid",
            color=COLORS_MARKERS[0],
            marker=MARKERS[0],
            markersize=10,
            label=f"{model_id} @stations: Mean +/- 1 std"
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
        if include_global:
            plt.plot(
                years,
                mod_latitude["mean"]["yearly"],
                linestyle="solid",
                color=COLORS_MARKERS[2],
                marker=MARKERS[2],
                markersize=10,
                label="{model_id}: Mean +/- 1 std" if l_i == 0 else None,
            )
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            )
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} annual mean "
            f"[{model_data.attributes['unit']}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        text = (
            f"No. of sites = {obs.coord('Station index (arbitrary)').shape[0]}"
        )
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)
        # Get legend handles for the first latitude band
        if l_i == 0:
            axes = plt.gcf().axes
            handles_mean, labels_mean = axes[0].get_legend_handles_labels()

        # Plot central-left column = scatter plot of model-obs
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
            f"NOAA GML {trace_gas.upper()} [{model_data.attributes['unit']}]",
            fontsize=20,
        )
        plt.ylabel(
            f"{model_id} {trace_gas.upper()} "
            f"[{model_data.attributes['unit']}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", which="major", labelsize=16, pad=10)
        plt.tick_params(axis="both", labelsize=16)
        text = f"$r^{2}$ = {linreg.rvalue**2:.2f}"
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)

        # Plot center-right column
        # multi-annual mean seasonal variation
        model_seas_anom_mean = np.array(
            [
                np.mean(
                    np.concatenate(
                        [trace_gas_at_model_anom[str(y)][m] for y in years],
                    ),
                )
                for m in months
            ]
        )
        model_seas_anom_std = np.array(
            [
                np.std(
                    np.concatenate(
                        [trace_gas_at_model_anom[str(y)][m] for y in years],
                    ),
                )
                for m in months
            ]
        )
        obs_seas_anom_mean = np.array(
            [
                np.mean(
                    np.concatenate(
                        [trace_gas_at_station_anom[str(y)][m] for y in years],
                    ),
                )
                for m in months
            ]
        )
        obs_seas_anom_std = np.array(
            [
                np.std(
                    np.concatenate(
                        [trace_gas_at_station_anom[str(y)][m] for y in years],
                    ),
                )
                for m in months
            ]
        )
        model_ts_mean = model_seas_anom_mean
        model_ts_std = model_seas_anom_std
        obs_ts_mean = obs_seas_anom_mean
        obs_ts_std = obs_seas_anom_std
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
        if include_global:
            plt.fill_between(
                np.arange(0, 12, 1),
                mod_latitude["mean"]["monthly"]
                - mod_latitude["std"]["monthly"],
                mod_latitude["mean"]["monthly"]
                + mod_latitude["std"]["monthly"],
                color=COLORS_MARKERS[2],
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
        if include_global:
            plt.plot(
                np.arange(0, 12, 1),
                mod_latitude["mean"]["monthly"],
                linestyle="solid",
                color=COLORS_MARKERS[2],
                marker=MARKERS[2],
                markersize=10,
            )
        plt.xticks(np.arange(0, 12, 1), months)
        plt.xlabel("Month", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} seasonal anomaly "
            f"[{model_data.attributes['unit']}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)

        # Plot right column
        # seasonal cycle timing
        seas_max_model = np.array(
            [
                np.argmax(
                    [trace_gas_at_model["max"][str(y)][m] for m in months]
                )
                for y in years
            ]
        )
        seas_min_model = np.array(
            [
                np.argmin(
                    [trace_gas_at_model["min"][str(y)][m] for m in months]
                )
                for y in years
            ]
        )
        seas_max_obs = np.array(
            [
                np.argmax(
                    [trace_gas_at_station["max"][str(y)][m] for m in months]
                )
                for y in years
            ]
        )
        seas_min_obs = np.array(
            [
                np.argmin(
                    [trace_gas_at_station["min"][str(y)][m] for m in months]
                )
                for y in years
            ]
        )
        month_offset = 4
        months = months[month_offset:] + months[:month_offset]
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
            label="" if l_i == 0 else None,
        )
        plt.plot(
            years,
            seas_max_obs,
            linestyle="solid",
            color=COLORS_MARKERS[1],
            marker=MARKERS[1],
            markersize=10,
            alpha=0.7,
            label="" if l_i == 0 else None,
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
            label="" if l_i == 0 else None,
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
            label="" if l_i == 0 else None,
        )
        if include_global:
            plt.scatter(
                years,
                1,
                linestyle="solid",
                color=COLORS_MARKERS[2],
                marker=MARKERS[2],
                markersize=10,
                label="{model_id}: Max" if l_i == 0 else None,
            )
            plt.scatter(
                years,
                1,
                linestyle="dashed",
                color=COLORS_MARKERS[2],
                marker=MARKERS[2],
                markersize=10,
                fillstyle="none",
                label="{model_id}: Min" if l_i == 0 else None,
            )
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            )
        )
        plt.xlabel("Year", fontsize=20)
        plt.yticks(np.arange(month_offset, month_offset + 12, 1), months)
        plt.ylim([month_offset - 1, month_offset + 12])
        plt.ylabel("Month", fontsize=20)
        plt.tick_params(axis="both", labelsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)
        # Get legend handles for the first latitude band
        if l_i == 0:
            axes = plt.gcf().axes
            handles_minmax, labels_minmax = axes[3].get_legend_handles_labels()
            labels_minmax[2] = f"{model_id} @stations: Max/Min"
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
    model_data,
    obs_cube,
    trace_gas,
    include_global=False,
):
    """Plot amplitude, growth rate and sensitivity between the two quantities.

    It uses different latitude slices for model and observational data:
    - 90S - 30S
    - 30S - 30N
    - 30N - 60N
    - 60N - 90N

    Parameters
    ----------
    model_data : Iris Cube.
        Contains model output of surface concentration with coordinates;
         time, latitude and longitude.
    obs_cube : Iris Cube.
        Contains information about NOAA GML measurement stations including
        station names, station latitude and station longitude.
    trace_gas : String.
        String containing the trace gas concentration to be plotted.
        Should be among CH4, CO2, or N2O.
    include_global : Bool
        Boolean flag to indicate if the global values should be overlaid.

    Returns
    -------
    figure : Figure object.
        Figure for the zonal time series of the trace gas for model data
        and observations of surface concentration from NOAA GML.
        The left column of the figure contains the zonal amplitude.
        The center column contains the zonal growth rate.
        The right column contains the zonal sensitivity amplitude/growth.
    """
    # Get model run id
    if "parent_source_id" in model_data.attributes:
        model_id = model_data.attributes["parent_source_id"]
    else:
        model_id = "Multi-Model-Mean"

    # Time coordinates
    try:
        iris.coord_categorisation.add_year(model_data, "time")
    except ValueError:
        err = f"Year coordinate already present in {model_id}."
        logger.debug(err)
    years = list(set(model_data.coord("year", dim_coords=False).points))

    # Latitude ranges
    latitude_titles = [
        r"Latitudes 60$^\circ$N - 90$^\circ$N",
        r"Latitudes 30$^\circ$N - 60$^\circ$N",
        r"Latitudes 30$^\circ$S - 30$^\circ$N",
        r"Latitudes 90$^\circ$S - 30$^\circ$S",
    ]
    latitude_ranges = {
        "60N - 90N": (60, 90),
        "30N - 60N": (30, 60),
        "30S - 30N": (-30, 30),
        "90S - 30S": (-90, -30),
    }

    # Model and observations preprocessing
    # Prepare iris cubes for model data
    model_yearly = model_data.aggregated_by("year", iris.analysis.MEAN)
    model_amplitude = model_data.aggregated_by(
        "year", iris.analysis.MAX
    ) - model_data.aggregated_by("year", iris.analysis.MIN)
    model_growth = iris.cube.Cube(
        np.diff(model_yearly.data, axis=0),
        long_name=f"{trace_gas}_growth",
        units=model_data.attributes["unit"],
        dim_coords_and_dims=[
            (
                iris.coords.DimCoord.from_coord(
                    model_yearly.coord("year")[1:]
                ),
                0,
            ),
            (model_data.coord("latitude"), 1),
            (model_data.coord("longitude"), 2),
        ],
    )
    model_growth.attributes = model_data.attributes
    # And observations
    obs_yearly = obs_cube.aggregated_by("year", iris.analysis.MEAN)
    obs_amplitude = obs_cube.aggregated_by(
        "year", iris.analysis.MAX
    ) - obs_cube.aggregated_by("year", iris.analysis.MIN)
    obs_growth = iris.cube.Cube(
        np.diff(obs_yearly.data, axis=0),
        long_name=f"{trace_gas}_growth",
        units=model_data.attributes["unit"],
        dim_coords_and_dims=[
            (iris.coords.DimCoord.from_coord(obs_yearly.coord("year")[1:]), 0),
            (obs_cube.coord("Station index (arbitrary)"), 1),
        ],
    )
    obs_growth.attributes = obs_cube.attributes
    for aux_coord_name in [
        "altitude",
        "latitude",
        "longitude",
        "platform_name",
    ]:
        aux_coord = obs_cube.coord(aux_coord_name)
        if obs_cube.coord_dims(aux_coord):
            obs_growth.add_aux_coord(
                aux_coord.copy(), obs_cube.coord_dims(aux_coord)
            )
    # Preprocess latitude slices
    model_slices = {}
    obs_slices = {}
    for key, (lat_min, lat_max) in latitude_ranges.items():
        constraint = iris.Constraint(
            latitude=lambda cell, lat_min=lat_min, lat_max=lat_max: lat_min
            <= cell
            < lat_max,
        )
        model_slices[key] = {}
        obs_slices[key] = {}
        # Yearly
        model_slices[key]["yearly"] = model_yearly.extract(constraint)
        obs_slices[key]["yearly"] = obs_yearly.extract(constraint)
        # Amplitude
        model_slices[key]["amplitude"] = model_amplitude.extract(constraint)
        obs_slices[key]["amplitude"] = obs_amplitude.extract(constraint)
        # Growth
        model_slices[key]["growth"] = model_growth.extract(constraint)
        obs_slices[key]["growth"] = obs_growth.extract(constraint)

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
    for l_i, lat_range in enumerate(latitude_ranges.keys()):
        # Select relevant latitude range
        obs_yearly_ls = obs_slices[lat_range]["yearly"]
        model_yearly_ls = model_slices[lat_range]["yearly"]
        obs_amplitude_ls = obs_slices[lat_range]["amplitude"]
        model_amplitude_ls = model_slices[lat_range]["amplitude"]
        obs_growth_ls = obs_slices[lat_range]["growth"]
        model_growth_ls = model_slices[lat_range]["growth"]
        # Compute aggregated values for the entire latitude range
        if include_global:
            mod_latitude = {
                "mean": {
                    "amplitude": np.zeros(shape=len(years)),
                    "growth": np.zeros(shape=len(years) - 1),
                    "sensitivity": np.zeros(shape=len(years) - 1),
                },
                "std": {
                    "amplitude": np.zeros(shape=len(years)),
                    "growth": np.zeros(shape=len(years) - 1),
                    "sensitivity": np.zeros(shape=len(years) - 1),
                },
            }
            weights_amp = iris.analysis.cartography.area_weights(
                model_amplitude_ls
            )
            mod_latitude["mean"]["amplitude"] = model_amplitude_ls.collapsed(
                ["latitude", "longitude"],
                iris.analysis.MEAN,
                weights=weights_amp,
            ).data
            mod_latitude["std"]["amplitude"] = weighted_std_dev(
                cube=model_amplitude_ls,
                dim=["latitude", "longitude"],
                weights=weights_amp,
            ).data
            tmp_relative_growth = iris.cube.Cube(
                100 * model_growth_ls.data / model_yearly_ls.data[1:, ...],
                long_name=f"{trace_gas}_relative_growth",
                units="%/yr",
                dim_coords_and_dims=[
                    (
                        iris.coords.DimCoord.from_coord(
                            model_growth_ls.coord("year")
                        ),
                        0,
                    ),
                    (model_growth_ls.coord("latitude"), 1),
                    (model_growth_ls.coord("longitude"), 2),
                ],
            )
            weights_growth = iris.analysis.cartography.area_weights(
                tmp_relative_growth
            )
            mod_latitude["mean"]["growth"] = tmp_relative_growth.collapsed(
                ["latitude", "longitude"],
                iris.analysis.MEAN,
                weights=weights_growth,
            ).data
            mod_latitude["std"]["growth"] = weighted_std_dev(
                cube=tmp_relative_growth,
                dim=["latitude", "longitude"],
                weights=weights_growth,
            ).data
            tmp_sensitivity = iris.cube.Cube(
                model_amplitude_ls.data[1:, ...] / model_growth_ls.data,
                long_name=f"{trace_gas}_sensitivity_ampl_grow",
                units=f"{model_data.attributes['unit']}.yr",
                dim_coords_and_dims=[
                    (
                        iris.coords.DimCoord.from_coord(
                            model_growth_ls.coord("year")
                        ),
                        0,
                    ),
                    (model_growth_ls.coord("latitude"), 1),
                    (model_growth_ls.coord("longitude"), 2),
                ],
            )
            weights_sensitivity = iris.analysis.cartography.area_weights(
                tmp_sensitivity
            )
            mod_latitude["mean"]["growth"] = tmp_sensitivity.collapsed(
                ["latitude", "longitude"],
                iris.analysis.MEAN,
                weights=weights_sensitivity,
            ).data
            mod_latitude["std"]["growth"] = weighted_std_dev(
                cube=tmp_sensitivity,
                dim=["latitude", "longitude"],
                weights=weights_sensitivity,
            ).data
        # Co-locate model grid points with measurement sites
        noaa_gml_lats = obs_yearly_ls.coord("latitude").points.tolist()
        noaa_gml_lons = obs_yearly_ls.coord("longitude").points.tolist()
        valid_obs = {
            "mean": {
                "amplitude": np.zeros(shape=len(years)),
                "growth": np.zeros(shape=len(years) - 1),
                "sensitivity": np.zeros(shape=len(years) - 1),
            },
            "std": {
                "amplitude": np.zeros(shape=len(years)),
                "growth": np.zeros(shape=len(years) - 1),
                "sensitivity": np.zeros(shape=len(years) - 1),
            },
        }
        valid_md = {
            "mean": {
                "amplitude": np.zeros(shape=len(years)),
                "growth": np.zeros(shape=len(years) - 1),
                "sensitivity": np.zeros(shape=len(years) - 1),
            },
            "std": {
                "amplitude": np.zeros(shape=len(years)),
                "growth": np.zeros(shape=len(years) - 1),
                "sensitivity": np.zeros(shape=len(years) - 1),
            },
        }
        for i, y in enumerate(years):
            obs_y = obs_yearly_ls.extract(iris.Constraint(year=y))
            obs_a = obs_amplitude_ls.extract(iris.Constraint(year=y))
            obs_g = obs_growth_ls.extract(iris.Constraint(year=y))
            # Extract yearly
            extracted_yearly = extract_pt(
                model_yearly_ls.extract(iris.Constraint(year=y)),
                noaa_gml_lats,
                noaa_gml_lons,
                nearest=True,
            )
            valid_indices_y = ~(obs_y.data.mask | np.isnan(extracted_yearly))
            v_obs_y = obs_y.data[valid_indices_y].tolist()
            v_mod_y = [
                tg.item()
                for j, tg in enumerate(extracted_yearly)
                if (tg.item() is not None) and valid_indices_y[j]
            ]
            # Extract amplitude
            extracted_amplitude = extract_pt(
                model_amplitude_ls.extract(iris.Constraint(year=y)),
                noaa_gml_lats,
                noaa_gml_lons,
                nearest=True,
            )
            valid_indices_a = ~(
                obs_a.data.mask | np.isnan(extracted_amplitude)
            )
            v_obs_a = obs_a.data[valid_indices_a]
            v_mod_a = [
                tg.item()
                for j, tg in enumerate(extracted_amplitude)
                if (tg.item() is not None) and valid_indices_a[j]
            ]
            valid_obs["mean"]["amplitude"][i] = np.mean(v_obs_a)
            valid_obs["std"]["amplitude"][i] = np.std(v_obs_a)
            valid_md["mean"]["amplitude"][i] = np.mean(v_mod_a)
            valid_md["std"]["amplitude"][i] = np.std(v_mod_a)
            # Extract growth & sensitivity = skip if first year as no data
            if i > 0:
                # Growth
                extracted_growth = extract_pt(
                    model_growth_ls.extract(iris.Constraint(year=y)),
                    noaa_gml_lats,
                    noaa_gml_lons,
                    nearest=True,
                )
                valid_indices_g = ~(
                    obs_g.data.mask | np.isnan(extracted_growth)
                )
                v_obs_g = obs_g.data[valid_indices_g]
                v_mod_g = [
                    tg.item()
                    for j, tg in enumerate(extracted_growth)
                    if (tg.item() is not None) and valid_indices_g[j]
                ]
                # Relative growth
                v_obs_g_relative = [
                    100 * g / v_obs_y[j - 1] for j, g in enumerate(v_obs_g)
                ]
                v_mod_g_relative = [
                    100 * g / v_mod_y[j - 1] for j, g in enumerate(v_mod_g)
                ]
                valid_obs["mean"]["growth"][i - 1] = np.mean(v_obs_g_relative)
                valid_obs["std"]["growth"][i - 1] = np.std(v_obs_g_relative)
                valid_md["mean"]["growth"][i - 1] = np.mean(v_mod_g_relative)
                valid_md["std"]["growth"][i - 1] = np.std(v_mod_g_relative)
                # Compute sensitivity between amplitude and growth
                v_obs_s = [
                    v_obs_a[j] / v_obs_g[j] for j in range(1, len(v_obs_g))
                ]
                v_mod_s = [
                    v_mod_a[j] / v_mod_g[j] for j in range(1, len(v_mod_g))
                ]
                valid_obs["mean"]["sensitivity"][i - 1] = np.mean(v_obs_s)
                valid_obs["std"]["sensitivity"][i - 1] = np.std(v_obs_s)
                valid_md["mean"]["sensitivity"][i - 1] = np.mean(v_mod_s)
                valid_md["std"]["sensitivity"][i - 1] = np.std(v_mod_s)
        # Plots
        # Plot left column = amplitude time series
        _ = plt.subplot(gs[l_i, 0])
        plt.plot(
            years,
            valid_md["mean"]["amplitude"],
            color=COLORS_MARKERS[0],
            linestyle="solid",
            marker=MARKERS[0],
            markersize=10,
            label=f"{model_id} @stations: Mean +/- 1 std"
            if l_i == 0
            else None,
        )
        plt.fill_between(
            years,
            valid_md["mean"]["amplitude"] - valid_md["std"]["amplitude"],
            valid_md["mean"]["amplitude"] + valid_md["std"]["amplitude"],
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.plot(
            years,
            valid_obs["mean"]["amplitude"],
            color=COLORS_MARKERS[1],
            linestyle="solid",
            marker=MARKERS[1],
            markersize=10,
            label="Observations: Mean +/- 1 std" if l_i == 0 else None,
        )
        plt.fill_between(
            years,
            valid_obs["mean"]["amplitude"] - valid_obs["std"]["amplitude"],
            valid_obs["mean"]["amplitude"] + valid_obs["std"]["amplitude"],
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        if include_global:
            plt.plot(
                years,
                mod_latitude["mean"]["amplitude"],
                color=COLORS_MARKERS[2],
                linestyle="solid",
                marker=MARKERS[2],
                markersize=10,
                label=f"{model_id}: Mean +/- 1 std" if l_i == 0 else None,
            )
            plt.fill_between(
                years,
                mod_latitude["mean"]["amplitude"]
                - mod_latitude["std"]["amplitude"],
                mod_latitude["mean"]["amplitude"]
                + mod_latitude["std"]["amplitude"],
                color=COLORS_MARKERS[2],
                alpha=0.1,
            )
        plt.xlim([years[0] - 1, years[-1] + 1])
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            )
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} seasonal amplitude "
            f"[{model_data.attributes['unit']}]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        text = "No. of sites = {}".format(
            obs_yearly_ls.coord("Station index (arbitrary)").shape[0],
        )
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)
        # Plot center column = relative growth
        _ = plt.subplot(gs[l_i, 1])
        plt.plot(
            years[1:],
            valid_md["mean"]["growth"],
            color=COLORS_MARKERS[0],
            linestyle="solid",
            marker=MARKERS[0],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            valid_md["mean"]["growth"] - valid_md["std"]["growth"],
            valid_md["mean"]["growth"] + valid_md["std"]["growth"],
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.plot(
            years[1:],
            valid_obs["mean"]["growth"],
            color=COLORS_MARKERS[1],
            linestyle="solid",
            marker=MARKERS[1],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            valid_obs["mean"]["growth"] - valid_obs["std"]["growth"],
            valid_obs["mean"]["growth"] + valid_obs["std"]["growth"],
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        if include_global:
            plt.plot(
                years[1:],
                mod_latitude["mean"]["growth"],
                color=COLORS_MARKERS[2],
                linestyle="solid",
                marker=MARKERS[2],
                markersize=10,
            )
            plt.fill_between(
                years[1:],
                mod_latitude["mean"]["growth"] - mod_latitude["std"]["growth"],
                mod_latitude["mean"]["growth"] + mod_latitude["std"]["growth"],
                color=COLORS_MARKERS[2],
                alpha=0.1,
            )
        plt.xlim([years[0] - 1, years[-1] + 1])
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            )
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            trace_gas.upper() + r" relative growth [$\%.yr^{-1}$]",
            fontsize=20,
        )
        plt.tick_params(axis="both", labelsize=16)
        text = "No. of sites = {}".format(
            obs_yearly_ls.coord("Station index (arbitrary)").shape[0],
        )
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)
        # Plot sensitivity between seasonal amplitude and growth
        _ = plt.subplot(gs[l_i, 2])
        plt.plot(
            years[1:],
            valid_md["mean"]["sensitivity"],
            color=COLORS_MARKERS[0],
            linestyle="solid",
            marker=MARKERS[0],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            valid_md["mean"]["sensitivity"] - valid_md["std"]["sensitivity"],
            valid_md["mean"]["sensitivity"] + valid_md["std"]["sensitivity"],
            color=COLORS_MARKERS[0],
            alpha=0.1,
        )
        plt.plot(
            years[1:],
            valid_obs["mean"]["sensitivity"],
            color=COLORS_MARKERS[1],
            linestyle="solid",
            marker=MARKERS[1],
            markersize=10,
        )
        plt.fill_between(
            years[1:],
            valid_obs["mean"]["sensitivity"] - valid_obs["std"]["sensitivity"],
            valid_obs["mean"]["sensitivity"] + valid_obs["std"]["sensitivity"],
            color=COLORS_MARKERS[1],
            alpha=0.1,
        )
        if include_global:
            plt.plot(
                years[1:],
                mod_latitude["mean"]["sensitivity"],
                color=COLORS_MARKERS[2],
                linestyle="solid",
                marker=MARKERS[2],
                markersize=10,
            )
            plt.fill_between(
                years[1:],
                mod_latitude["mean"]["sensitivity"]
                - mod_latitude["std"]["sensitivity"],
                mod_latitude["mean"]["sensitivity"]
                + mod_latitude["std"]["sensitivity"],
                color=COLORS_MARKERS[2],
                alpha=0.1,
            )
        plt.xlim([years[0] - 1, years[-1] + 1])
        plt.xticks(
            np.arange(
                years[0],
                years[-1] + 1,
                1 if years[-1] - years[0] < 10 else 4,
            )
        )
        plt.xlabel("Year", fontsize=20)
        plt.ylabel(
            f"{trace_gas.upper()} sensitivity\n(amplitude/growth)", fontsize=20
        )
        plt.tick_params(axis="both", labelsize=16)
        text = "No. of sites = {}".format(
            obs_yearly_ls.coord("Station index (arbitrary)").shape[0],
        )
        plt.annotate(text, (0.05, 0.9), xycoords="axes fraction", fontsize=16)
        plt.title(latitude_titles[l_i], fontsize=28)

    figure.legend(
        loc="upper center",
        bbox_to_anchor=[0.5, 0.08],
        ncol=3,
        fontsize=20,
        borderaxespad=0.01,
    )

    return figure


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
    obs_dataset : Dictionary
        ESMValTool dictionary. Holds meta data for the observational
        trace gas dataset.
    config : Dictionary
        ESMValTool recipe configuration.

    Returns
    -------
     multiannual_seaonal_mean : iris.cube.Cube. Preprocessed observational
        climatology of the trace gas.
    """
    obs_cube = iris.load_cube(obs_dataset[0]["filename"])

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
        data=ma.masked_where(np.isnan(obs_cube.data), obs_cube.data)
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
        )
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
        list(set.intersection(*valid_station_indices_for_seasons))
    )
    valid_stations = obs_cube.coord(
        "Station index (arbitrary)",
    ).points[valid_stations_indices]

    # Only keep stations for which these conditions are met
    station_indices = obs_cube.coord("Station index (arbitrary)").points
    mask = np.isin(station_indices, valid_stations)
    obs_cube_valid = obs_cube.extract(
        iris.Constraint(
            coord_values={
                "Station index (arbitrary)": lambda cell: mask[
                    int(cell.point)
                ],
            },
        ),
    )

    return obs_cube_valid, multi_annual_seasonal_mean


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
    obs_dataset = datasets.pop(config["observational_dataset"])
    obs_cube, obs_cube_seasonal = preprocess_obs_dataset(obs_dataset, config)
    logger.info(obs_cube)

    datasets_preproc = {}
    datasets_seasonal = {}

    # Default values for the parameters plots and only_multimodel:
    #   - plots = "seas_maps", "timeserie_lat", "sensitivity_ampl_trend"
    #   - only_multimodel = False
    if "plots" not in config:
        config["plots"] = [
            "seas_maps",
            "timeserie_lat",
            "sensitivity_ampl_trend",
        ]
    if "only_multimodel" not in config:
        config["only_multimodel"] = False

    # Only consider the multi-model mean
    if config["only_multimodel"]:
        model_dataset = "MultiModelMean"
        group = datasets[model_dataset]
        # "model_dataset" is the name of the model dataset.
        # "group" is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        logger.info(group)

        for attributes in group:
            logger.info(attributes["filename"])
            input_file = attributes["filename"]
            cube = iris.load_cube(input_file)

            # Add bounds for lat and lon if not present
            cube = add_bounds(cube)

            # Put observations and model data on same scale for ppm/ppb
            cube = TRACE_GASES_FACTOR[config["trace_gas"]] * cube
            # Change units accordingly
            cube.attributes["unit"] = TRACE_GASES_UNITS[config["trace_gas"]]

            # Process multi-annual seasonal mean
            cube_seasonal = climate_statistics(
                cube,
                operator="mean",
                period="season",
                seasons=["DJF", "MAM", "JJA", "SON"],
            )

            # Store datasets for multi-model average
            datasets_preproc[model_dataset] = {
                "attributes": attributes,
                "cube": cube,
            }
            datasets_seasonal[model_dataset] = {
                "attributes": attributes,
                "cube": cube_seasonal,
            }

            # Set up for analysis and plotting
            seasons = ["DJF", "MAM", "JJA", "SON"]
            plot_file_prefix = (
                model_dataset
                + "_"
                + attributes["activity"]
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

            # Get activity, mip, exp, variable , time range
            attrs = datasets_preproc[model_dataset]["attributes"]
            plot_file_prefix = (
                "Multi_model_mean_"
                + attrs["activity"]
                + "_"
                + attrs["mip"]
                + "_"
                + attrs["exp"]
                + "_"
                + attrs["short_name"]
                + "_"
                + str(attrs["start_year"])
                + "_"
                + str(attrs["end_year"])
                + "_"
            )
            # Compute multi-model averages
            multi_model = multi_model_statistics(
                products=[
                    datasets_preproc[m]["cube"] for m in datasets_preproc
                ],
                span="overlap",
                statistics=["mean", "std_dev"],
            )
            multi_model_seas = multi_model_statistics(
                products=[
                    datasets_seasonal[m]["cube"] for m in datasets_preproc
                ],
                span="overlap",
                statistics=["mean", "std_dev"],
            )

            if "seas_maps" in config["plots"]:
                # Analysis and plotting for model-obs seasonal map comparison
                figures_seas_map, fig_scatter_seas = trace_gas_maps(
                    multi_model_seas["mean"],
                    obs_cube_seasonal,
                    seasons,
                    attrs["timerange"],
                    config["trace_gas"],
                )
                # Save scatter plot and contour plots for model-obs seasonal comparison
                output_file = plot_file_prefix + "scatter"
                output_path = get_plot_filename(output_file, config)
                fig_scatter_seas.savefig(output_path, bbox_inches="tight")
                output_file = plot_file_prefix + "seas_map"
                output_path = get_plot_filename(output_file, config)
                figures_seas_map.savefig(output_path, bbox_inches="tight")

            if "timeserie_lat" in config["plots"]:
                # Analysis and plotting for zonal-mean time series for model-obs
                figure_timeserie_zonal = trace_gas_timeserie_zonal(
                    multi_model["mean"],
                    obs_cube,
                    config["trace_gas"],
                )
                # Save time series plots per latitude range
                output_file = plot_file_prefix + "timeseries_latitude"
                output_path = get_plot_filename(output_file, config)
                figure_timeserie_zonal.savefig(
                    output_path, bbox_inches="tight"
                )

            if "sensitivity_ampl_trend" in config["plots"]:
                # Analysis and plotting for sensitivity between
                # seasonal cycle amplitude and growth rate
                figure_sens_ampl_trend = trace_gas_seas_ampl_growth_rate(
                    multi_model["mean"],
                    obs_cube,
                    config["trace_gas"],
                )
                # Save sensitivity seasonal amplitude and trend plot
                output_file = plot_file_prefix + "seas_amplitude_trend"
                output_path = get_plot_filename(output_file, config)
                figure_sens_ampl_trend.savefig(
                    output_path, bbox_inches="tight"
                )

    # or run for each seperate model + multi-model mean
    else:
        for model_dataset, group in datasets.items():
            # "model_dataset" is the name of the model dataset.
            # "group" is a list of dictionaries containing metadata.
            logger.info("Processing data for %s", model_dataset)
            logger.info(group)

            for attributes in group:
                logger.info(attributes["filename"])
                input_file = attributes["filename"]
                cube = iris.load_cube(input_file)

                # Add bounds for lat and lon if not present
                cube = add_bounds(cube)

                # Put observations and model data on same scale for ppm/ppb
                cube = TRACE_GASES_FACTOR[config["trace_gas"]] * cube
                # Change units accordingly
                cube.attributes["unit"] = TRACE_GASES_UNITS[
                    config["trace_gas"]
                ]

                # Process multi-annual seasonal mean
                cube_seasonal = climate_statistics(
                    cube,
                    operator="mean",
                    period="season",
                    seasons=["DJF", "MAM", "JJA", "SON"],
                )

                # Store datasets for multi-model average
                datasets_preproc[model_dataset] = {
                    "attributes": attributes,
                    "cube": cube,
                }
                datasets_seasonal[model_dataset] = {
                    "attributes": attributes,
                    "cube": cube_seasonal,
                }

                # Set up for analysis and plotting
                seasons = ["DJF", "MAM", "JJA", "SON"]
                plot_file_prefix = (
                    model_dataset
                    + "_"
                    + attributes["activity"]
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
                    # Analysis and plotting for model-obs seasonal map comparison
                    figures_seas_map, fig_scatter_seas = trace_gas_maps(
                        cube_seasonal,
                        obs_cube_seasonal,
                        seasons,
                        attributes["timerange"],
                        config["trace_gas"],
                    )

                if "timeserie_lat" in config["plots"]:
                    # Analysis and plotting for zonal-mean time series for model-obs
                    figure_timeserie_zonal = trace_gas_timeserie_zonal(
                        cube,
                        obs_cube,
                        config["trace_gas"],
                    )

                if "sensitivity_ampl_trend" in config["plots"]:
                    # Analysis and plotting for sensitivity between
                    # seasonal cycle amplitude and growth rate
                    figure_sens_ampl_trend = trace_gas_seas_ampl_growth_rate(
                        cube,
                        obs_cube,
                        config["trace_gas"],
                    )

            # Saving plots
            if "seas_maps" in config["plots"]:
                # Save scatter plot and contour plots for model-obs seasonal comparison
                output_file = plot_file_prefix + "scatter"
                output_path = get_plot_filename(output_file, config)
                fig_scatter_seas.savefig(output_path)
                output_file = plot_file_prefix + "seas_map"
                output_path = get_plot_filename(output_file, config)
                figures_seas_map.savefig(output_path, bbox_inches="tight")
            if "timeserie_lat" in config["plots"]:
                # Save time series plots per latitude range
                output_file = plot_file_prefix + "timeseries_latitude"
                output_path = get_plot_filename(output_file, config)
                figure_timeserie_zonal.savefig(
                    output_path, bbox_inches="tight"
                )
            if "sensitivity_ampl_trend" in config["plots"]:
                # Save sensitivity seasonal amplitude and trend plot
                output_file = plot_file_prefix + "seas_amplitude_trend"
                output_path = get_plot_filename(output_file, config)
                figure_sens_ampl_trend.savefig(
                    output_path, bbox_inches="tight"
                )


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
