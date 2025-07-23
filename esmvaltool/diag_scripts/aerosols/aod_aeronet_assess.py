"""Implement the AOD climatology metric from ground-based
AeroNet observations.
"""

import logging
import os

import iris
import iris.plot as iplt
import matplotlib.cm as mpl_cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy
from matplotlib import colors, gridspec
from numpy import ma

from esmvaltool.diag_scripts.aerosols.aero_utils import add_bounds, extract_pt
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename

logger = logging.getLogger(os.path.basename(__file__))
fontsizedict = {"title": 25, "axis": 20, "legend": 18, "ticklabel": 18}


def get_provenance_record(filenames):
    """Return a provenance record describing the metric.

    Parameters
    ----------
    filenames : List of strings
        The filenames containing the data used to create the metric.

    Returns
    -------
    dictionary
        The provenance record describing the metric.
    """
    record = {
        "ancestors": filenames,
    }

    return record


def plot_aod_mod_obs(md_data, obs_data, aeronet_obs_cube, plot_dict):
    """Plot AOD contour overlaid with Aeronet climatology.

    Parameters
    ----------
    md_data : Iris cube
        Model AOD as a cube with latitude and longitude coordinates.
    obs_data : List.
        Observations of AOD from each AeroNET station.
    aeronet_obs_cube : Iris cube.
        Holds information about Aeronet measurement stations including
        station names, station latitude and station longitude.
    plot_dict : Dictionary.
        Holds plotting settings.
    """
    # Plot model data
    cf_plot = iplt.contourf(
        md_data, plot_dict["Levels"], colors=plot_dict["Colours"], extend="max"
    )

    # Latitude and longitude of stations.
    anet_aod_lats = aeronet_obs_cube.coord("latitude").points
    anet_aod_lons = (
        aeronet_obs_cube.coord("longitude").points + 180
    ) % 360 - 180

    # Loop over stations
    for istn, stn_data in enumerate(obs_data):
        if ma.is_masked(stn_data):
            continue

        # Find position of the observed AOD on the colorscale.
        # np.searchsorted returns index at which inserting new value will
        # maintain a sorted array. We use the color to the left of index.
        cid = np.searchsorted(plot_dict["Levels"], stn_data)
        cid = max(0, cid - 1)  # filter out zero and max when seeking 'left'
        cid = min(len(plot_dict["Colours"]) - 1, cid)
        pcol = plot_dict["Colours"][cid]

        # Overlay contourf with observations
        plt.plot(
            anet_aod_lons[istn],
            anet_aod_lats[istn],
            color=pcol,
            marker="o",
            markeredgecolor="k",
            markeredgewidth=2,
            markersize=9,
        )

    # Decorate the plot
    plt.title(plot_dict["Title"], size=24)
    colbar = plt.colorbar(cf_plot, orientation="horizontal")
    colbar.set_ticks(plot_dict["Levels"])
    colbar.set_ticklabels(plot_dict["tick_labels"])
    plt.gca().coastlines(color="#525252")

    # Statistics on plot
    plt.figtext(
        0.12,
        0.27,
        (
            f"""Global mean AOD={plot_dict["Mean_aod"]:.3f}; RMSE="""
            f"""{plot_dict["RMS_aod"]:.3f}; Stn mean: md="""
            f"""{plot_dict["Stn_mn_md"]:.3f}; obs="""
            f"""{plot_dict["Stn_mn_obs"]:.3f}"""
        ),
        size=16,
    )


def aod_analyse(model_data, aeronet_obs_cube, clim_seas, wavel):
    """Evaluates AOD vs Aeronet, generates plots and returns evaluation
    metrics.

    Parameters
    ----------
    model_data : Iris Cube.
        Contains model output of AOD with coordinates; time, latitude and
        longitude.
    aeronet_obs_cube : Iris Cube.
        Contains information about Aeronet measurement stations including
        station names, station latitude and station longitude.
    clim_seas : List.
       Strings to denote climate seasons ["DJF", "MAM", "JJA", "SON"]
    wavel : String.
        AOD wavelength, default = 440nm - translates to pseudo-level.

    Returns
    -------
    figures : List.
        Contains figure instances for the seasonal contour plots overlaid with
        observations of AOD from AeroNET.
    fig_scatter : Figure object.
        The scatter plot comparing modelled and observed AOD at 440nm.
    """
    # Convert wave length nm -> um
    wv_mi = str(float(wavel) / 1000.0)

    # Get model run id
    if "parent_source_id" in model_data.attributes:
        model_id = model_data.attributes["parent_source_id"]
    else:
        model_id = "Multi-Model-Mean"

    # Add bounds for lat and lon if not present
    model_data = add_bounds(model_data)

    # Co-locate model grid points with measurement sites --func from aero_utils
    anet_aod_lats = aeronet_obs_cube.coord("latitude").points.tolist()
    anet_aod_lons = aeronet_obs_cube.coord("longitude").points.tolist()
    aod_at_anet = extract_pt(model_data, anet_aod_lats, anet_aod_lons)

    # Set up seasonal contour plots
    figures = []

    clevs = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 2.0]
    clabs = [
        "0.0",
        "",
        "0.1",
        "",
        "0.2",
        "",
        "0.3",
        "",
        "0.4",
        "",
        "0.5",
        "2.0",
    ]
    cmapr = mpl_cm.get_cmap("brewer_Spectral_11")
    cmap = colors.ListedColormap(cmapr(range(cmapr.N))[-1::-1])
    colours = cmap.colors

    # Set up the figure for scatter plotting
    fig_scatter = plt.figure(figsize=(10, 10))
    gs_scatter = gridspec.GridSpec(ncols=1, nrows=1)
    ax_scatter = fig_scatter.add_subplot(gs_scatter[0, 0])
    col_scatter = ["#081d58", "#41ab5d", "#fe9929", "#7f0000"]
    leg_scatter = []

    # Loop over seasons
    for season in aeronet_obs_cube.slices_over("clim_season"):
        # Match Aeronet obs season with model season number
        model_sn = [c.lower() for c in clim_seas].index(
            season.coord("clim_season").points[0]
        )
        model_season = model_data[model_sn]

        logger.info(
            "Analysing AOD for %s: %s", {model_id}, {clim_seas[model_sn]}
        )

        # Generate statistics required - area-weighted mean
        grid_areas = iris.analysis.cartography.area_weights(model_season)
        global_mean = model_season.collapsed(
            ["latitude", "longitude"],
            iris.analysis.MEAN,
            weights=grid_areas,
        )

        # Extract model and obs data for season number (model_sn)
        seas_anet_obs = season.data
        seas_anet_md = np.array([x[model_sn] for x in aod_at_anet])

        # Match model data with valid obs data
        valid_indices = ma.where(seas_anet_obs)
        valid_obs = seas_anet_obs[valid_indices]
        valid_md = seas_anet_md[valid_indices]

        # Model - obs statistics (diff, model mean and RMS, r2)
        diff = valid_md - valid_obs
        stn_mn_obs = np.mean(valid_obs)
        stn_mn_md = np.mean(valid_md)
        rms_aod = np.sqrt(np.mean(diff**2))
        linreg = scipy.stats.linregress(valid_obs, valid_md)

        # Plot scatter of co-located model and obs data
        ax_scatter.scatter(valid_obs, valid_md, color=col_scatter[model_sn])

        # Legend
        label = f"{clim_seas[model_sn]} = {linreg.rvalue**2:.2f}"
        leg_scatter.append(
            mlines.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=label,
                markersize=15,
                markerfacecolor=col_scatter[model_sn],
            )
        )

        # Plot contours overlaid with obs for this run and season
        fig_cf = plt.figure(figsize=(11, 8), dpi=300)

        n_stn = str(len(valid_obs))
        title = (
            "\nTotal Aerosol Optical Depth at "
            + wv_mi
            + " microns"
            + "\n"
            + model_id
            + ", "
            + clim_seas[model_sn]
            + ", N stations="
            + n_stn
        )

        # Plot dictionary
        plot_dict = {
            "Mean_aod": global_mean.data,
            "Stn_mn_obs": stn_mn_obs,
            "Stn_mn_md": stn_mn_md,
            "RMS_aod": rms_aod,
            "Levels": clevs,
            "Colours": colours,
            "tick_labels": clabs,
            "Title": title,
            "Season": clim_seas[model_sn],
        }
        plot_aod_mod_obs(
            model_season, seas_anet_obs, aeronet_obs_cube, plot_dict
        )

        figures.append(fig_cf)

    # Decorate the scatter plot
    line = mlines.Line2D([0, 1], [0, 1], color="#696969")
    transform = ax_scatter.transAxes
    line.set_transform(transform)
    ax_scatter.add_line(line)

    ax_scatter.set(
        xlim=(0, 1),
        xticks=np.linspace(0.0, 1.0, num=6),
        ylim=(0, 1),
        yticks=np.linspace(0.0, 1.0, num=6),
    )
    ax_scatter.set_xlabel("AeroNET AOD", fontsize=fontsizedict["axis"])
    ax_scatter.set_ylabel(model_id + " AOD", fontsize=fontsizedict["axis"])

    ax_scatter.tick_params(
        axis="both", which="major", labelsize=fontsizedict["ticklabel"]
    )

    ax_scatter.set_title(
        "Model vs obs: Total Aerosol Optical Depth \n at "
        + wv_mi
        + " microns",
        fontsize=fontsizedict["title"],
    )

    ax_scatter.legend(
        handles=leg_scatter,
        loc="lower right",
        title="Seasonal R2",
        title_fontsize=fontsizedict["legend"],
        fontsize=fontsizedict["legend"],
    )

    return figures, fig_scatter


def preprocess_aod_obs_dataset(obs_dataset):
    """Calculate a multiannual seasonal mean AOD climatology.

    Observational AOD timeseries data from AeroNET are used to generate a
    multiannual seasonal mean climatology for each AeroNET station. The
    user sets thresholds (or uses the default settings) to specify the
    amount of valid data required for the climatology. At this stage
    ESMValTool preprocessors are unsuitable for pre-processing the AeroNET
    AOD observations because of the bespoke nature and application of the
    filtering thresholds.

    Parameters
    ----------
    obs_dataset : ESMValTool dictionary. Holds meta data for the observational
        AOD dataset.

    Returns
    -------
     multiannual_seaonal_mean : Iris cube. Preprocessed observational
         AOD climatology.
    """
    obs_cube = iris.load_cube(obs_dataset[0]["filename"])

    # Set up thresholds for generating the multi annual seasonal mean
    min_days_per_mon = 1
    min_mon_per_seas = 3
    min_seas_per_year = 4
    min_seas_per_clim = 5

    # Add the clim_season and season_year coordinates.
    iris.coord_categorisation.add_year(obs_cube, "time", name="year")

    iris.coord_categorisation.add_season(obs_cube, "time", name="clim_season")

    iris.coord_categorisation.add_season_year(
        obs_cube, "time", name="season_year"
    )

    # Copy obs cube and mask all months with fewer
    # "Number of days" than given threshold.
    num_days_var = obs_cube.ancillary_variable("Number of days")
    masked_months_obs_cube = obs_cube.copy(
        data=ma.masked_where(
            num_days_var.data < min_days_per_mon, obs_cube.data
        )
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
        len(multi_annual_seasonal_mean.coord("clim_season").points)
    )
    for iseas in counter:
        multi_annual_seasonal_mean.data[iseas, :] = ma.masked_where(
            year_agg_count.data[0, :] < min_seas_per_year,
            multi_annual_seasonal_mean.data[iseas, :],
        )

    return multi_annual_seasonal_mean


def main(config):
    """Produce the AOD climatology metric from ground-based AeroNet
    observations.

    Parameters
    ----------
    wavel : String.
        User defined. Default is "440".
    config : dict
        The ESMValTool configuration.
    """
    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    # Default wavelength
    wavel = "440"

    # Produce climatology for observational dataset
    obs_dataset = datasets.pop(config["observational_dataset"])
    obs_cube = preprocess_aod_obs_dataset(obs_dataset)

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        logger.info(group)

        for attributes in group:
            logger.info(attributes["filename"])

            input_file = attributes["filename"]
            provenance_record = get_provenance_record(input_file)
            logger.info(provenance_record)
            cube = iris.load_cube(input_file)

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

            # Analysis and plotting for model-obs comparison
            figures, fig_scatter = aod_analyse(
                cube, obs_cube, seasons, wavel=wavel
            )

        # Save the scatter plot
        output_file = plot_file_prefix + "scatter"
        output_path = get_plot_filename(output_file, config)
        fig_scatter.savefig(output_path)

        # Save the contour plots
        for ifig, seas_fig in enumerate(figures):
            output_file = plot_file_prefix + seasons[ifig]
            output_path = get_plot_filename(output_file, config)
            seas_fig.savefig(output_path)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
