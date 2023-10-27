"""Implement the AOD climatology metric from ground-based AeroNet
observations."""

import logging

import iris
import iris.plot as iplt
import matplotlib.cm as mpl_cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy
from aero_utils import AeroAnsError, add_bounds, extract_pt
from matplotlib import colors, gridspec
from numpy import ma

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename

fontsizedict = {"title": 25, "axis": 20, "legend": 18, "ticklabel": 18}


def get_provenance_record(filenames):
    """Return a provenance record describing the metric.

    Parameters
    ----------
    filenames : list of strings
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


def plot_aod_mod_obs(md_data, obs_data, aeronet_ref_cube, plot_dict):
    """Plot AOD contour overlaid with Aeronet values for given period.

    Parameters
    ----------
    md_data : Iris cube
        Model AOD as a cube with latitude and longitude coordinates.
    obs_data : List.
        Observations of AOD from each AeroNET station.
    aeronet_ref_cube : Iris cube.
        Contains information about Aeronet measurement stations including
        station names, station latitude and station longitude.
    plot_dict : Dictionary.
        Contains plotting settings.
    """

    # Plot model data
    cf_plot = iplt.contourf(md_data,
                            plot_dict["Levels"],
                            colors=plot_dict["Colours"],
                            extend="max")

    # Latitude and longitude of stations.
    anet_aod_lats = aeronet_ref_cube.coord("latitude").points
    anet_aod_lons = (
        (aeronet_ref_cube.coord("longitude").points + 180) % 360 - 180
    )

    # Loop over stations
    for istn, stn in enumerate(aeronet_ref_cube.coord("platform_name").points):
        if obs_data.mask[istn]:
            continue

        # Find position of the observed AOD on the colorscale.
        # np.searchsorted returns index at which inserting new value will
        # maintain a sorted array. We use the color to the left of index.
        cid = np.searchsorted(plot_dict["Levels"], obs_data[istn])
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
        0.25,
        0.30,
        "MEAN: {0:6.3f}. RMSE v Aeronet: {1:6.3f} ".format(
            plot_dict["Mean_aod"], plot_dict["RMS_aod"]),
        size=20,
    )


def read_aeronet_clim(aeronet_dir, years):
    """Read in pre-processed AERONET AOD values.

    Method
    ------
        The processed AERONET data is held as NetCDF files. An
        individual file contains a time series of pre-processed monthly
        mean AOD data from AeroNET observations at a single station.
        The meta data holds information about data's provenance, the Station
        (e.g. name, altitude) and about AeroNet. The filenames have the
        structure: OBS_AERONET_ground_*_AERmon_od440aer_199301-202201.nc
        where * = the station id, e.g. Mauna_Loa.

    Parameters
    ----------
    aeronet_dir : String. The file path to the AeroNET observation data files.
    years : List of strings for the start and end year of the time series of
        interest

    Returns
    -------
    anet_aod : Dictionary.
        Contains station names, station latitude, station longitude and
        multiannual seasonal mean AOD climatologies for each AeroNET station.
    """

    # Open observational AOD data files
    try:
        cubes = iris.load(aeronet_dir + "*")

    except AeroAnsError:
        print('READ_AERONET:Error loading Aeronet data: "{}"'.format(
            aeronet_dir))

    # Lists to hold station data
    stnlist = []
    lats = []
    lons = []
    aod_seas = []

    # Set up time constraint
    start_yr = int(years[0]) - 1
    end_yr = int(years[1])
    pdt1 = iris.time.PartialDateTime(year=start_yr, month=11, day=1)
    pdt2 = iris.time.PartialDateTime(year=end_yr, month=12, day=1)
    time_constr = iris.Constraint(time=lambda cell: pdt1 < cell.point < pdt2)

    # Loop over station files
    for cube in cubes:

        cube = iris.util.squeeze(cube)

        # Mask missing data
        masked_data = ma.masked_equal(cube.data, -999.0, copy=True)
        points_to_mask = ma.getmask(masked_data)
        iris.util.mask_cube(cube, points_to_mask)

        # Extract station name and lat lon
        stnlist.append(cube.attributes["station_id"])
        lats.append(cube.coord("latitude").points[0])
        lons.append(cube.coord("longitude").points[0])

        # Add coordinates to obs data
        iris.coord_categorisation.add_season(cube, "time", name="clim_season")
        iris.coord_categorisation.add_season_year(cube,
                                                  "time",
                                                  name="season_year")

        # Extract time period of interest
        cube = cube.extract(time_constr)

        # Calculate multi-annual seasonal means
        season = cube.aggregated_by(["clim_season"], iris.analysis.MEAN)
        aod_seas.append(season.data)

    # Analyse and plot the AOD metrics
    anet_aod = {
        "Stations": stnlist,
        "Lats": lats,
        "Lons": lons,
        "MA_season": aod_seas,
    }

    return anet_aod


def aod_analyse(md_data, aeronet_ref_cube, clim_seas, years, wavel):
    """Evaluates AOD vs Aeronet, generates plots and returns evaluation
    metrics.

    Parameters
    ----------
    md_data : Iris Cube.
        Contains model output of AOD with coordinates; time, latitude and
        longitude.
    aeronet_dir : String.
        The file path for the AeroNET observations stored as NetCDF files.
    seas : List.
       Strings to denote climate seasons ["DJF", "MAM", "JJA", "SON"]
    wavel : String.
        AOD wavelength, default = 440nm - translates to pseudo-lev.

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

    # Dictionary for output metric/s
    metric = {}

    # Get model run id
    md_id = md_data.attributes["parent_source_id"]

    # Add bounds for lat and lon if not present
    md_data = add_bounds(md_data)

    # Co-locate model grid points with measurement sites --func from aero_utils
    anet_aod_lats = aeronet_ref_cube.coord("latitude").points.tolist()
    anet_aod_lons = aeronet_ref_cube.coord("longitude").points.tolist()
    aod_at_anet = extract_pt(md_data, anet_aod_lats, anet_aod_lons)

    # Set up seasonal contour plots
    figures = []

    clevs = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 2.0]
    clabs = [
        "0.0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5", "2.0"
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
    for iseas, seas in enumerate(md_data.slices_over("season_number")):
        print("Analysing AOD for ", md_id, ":", clim_seas[iseas])

        # Generate statistics required - area-weighted mean
        grid_areas = iris.analysis.cartography.area_weights(seas)
        global_mean = seas.collapsed(["latitude", "longitude"],
                                     iris.analysis.MEAN,
                                     weights=grid_areas)

        # Extract model and obs data for iseas
        seas_anet_obs = aeronet_ref_cube.data[iseas]
        seas_anet_md = np.array([x[iseas] for x in aod_at_anet])

        # Match model data with valid obs data
        valid_indices = ma.where(seas_anet_obs)
        valid_obs = seas_anet_obs[valid_indices]
        valid_md = seas_anet_md[valid_indices]

        # Model - obs statistics (diff, model mean and RMS, r2)
        diff = valid_md - valid_obs
        rms_aod = np.sqrt(np.mean(diff**2))
        linreg = scipy.stats.linregress(valid_obs, valid_md)

        # Populate metric
        metric["Aerosol optical depth at " + wv_mi + " um " +
               clim_seas[iseas]] = rms_aod

        # Plot scatter of co-located model and obs data
        ax_scatter.scatter(valid_obs, valid_md, color=col_scatter[iseas])

        # Legend
        label = clim_seas[iseas] + " = " + str("%.2f" % linreg.rvalue**2)
        leg_scatter.append(
            mlines.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=label,
                markersize=15,
                markerfacecolor=col_scatter[iseas],
            ))

        # Plot contours overlaid with obs for this run and season
        fig_cf = plt.figure(figsize=(11, 8), dpi=300)

        n_stn = str(len(valid_obs))
        title = (md_id + ", N stations = " + n_stn + ", " + clim_seas[iseas] +
                 "\nTotal Aerosol Optical Depth at " + wv_mi + " microns")

        # Plot dictionary
        plot_dict = {
            "Mean_aod": global_mean.data,
            "RMS_aod": rms_aod,
            "Levels": clevs,
            "Colours": colours,
            "tick_labels": clabs,
            "Title": title,
            "Season": clim_seas[iseas],
        }
        plot_aod_mod_obs(seas, seas_anet_obs, aeronet_ref_cube, plot_dict)

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
    ax_scatter.set_ylabel(md_id + " AOD", fontsize=fontsizedict["axis"])
    ax_scatter.tick_params(axis="both",
                           which="major",
                           labelsize=fontsizedict["ticklabel"])

    ax_scatter.set_title(
        "Model vs obs: Total Aerosol Optical Depth \n at " + wv_mi +
        " microns",
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


def main(config):
    """Produce the AOD climatology metric from ground-based AeroNet
    observations.

    Parameters
    ----------
    wavel : String.
        User defined. Default is "440".
    aeronet_dir : String.
        The directory containing the Aeronet observational climatology
        This is currently set by the user because the Aeronet climatologies
        are not yet CMORized or ready for use with the ESMValTool
        pre-processors.
    config : dict
        The ESMValTool configuration.
    """

    logger = logging.getLogger(__name__)

    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    # Default wavelength
    wavel = "440"

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        logger.info(group)

        for attributes in group:
            aeronet_ref = attributes["reference_dataset"]
            if aeronet_ref == attributes["dataset"]:
                continue

            aeronet_ref_cube = iris.load_cube(
                datasets[aeronet_ref][0]["filename"]
            )

            logger.info(attributes["filename"])

            input_file = attributes["filename"]
            provenance_record = get_provenance_record(input_file)
            logger.info(provenance_record)
            cube = iris.load_cube(input_file)

            # Set up for analysis and plotting
            seasons = ["DJF", "MAM", "JJA", "SON"]
            years = [attributes["start_year"], attributes["end_year"]]

            plot_file_prefix = (model_dataset + "_" + attributes["activity"] +
                                "_" + attributes["mip"] + "_" +
                                attributes["exp"] + "_" +
                                attributes["short_name"] + "_" +
                                attributes["grid"] + "_" +
                                str(attributes["start_year"]) + "_" +
                                str(attributes["end_year"]) + "_")

            # Analysis and plotting for model-obs comparison
            figures, fig_scatter = aod_analyse(cube,
                                               aeronet_ref_cube,
                                               seasons,
                                               years,
                                               wavel=wavel)

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
