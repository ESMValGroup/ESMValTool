"""Implement the AOD climatology metric from ground-based AeroNet
observations."""

import logging

import iris
import iris.plot as iplt
import matplotlib.cm as mpl_cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from aero_utils import AeroAnsError, add_bounds, extract_pt

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename


def plot_aod_mod_obs(md_data, obs_data, anet_aod, plot_dict):
    """Plot AOD contour overlaid with Aeronet values for given period.

    Parameters
    ----------
    md_data : Model AOD as a cube with latitude and longitude coordinates
    obs_data : List of AOD observations at each station
    anet_aod : Dictionary for AOD observations data.
        The dictionary includes station names, station latitude, station
        longitude and model-obs comparison statistics.
    plot_dict : Dictionary of plotting settings
    """

    # Plot model data
    cf = iplt.contourf(md_data,
                       plot_dict["Levels"],
                       colors=plot_dict["Colours"])

    # Loop over stations
    for istn, stn in enumerate(anet_aod["Stations"]):
        # Find position of the observed AOD on the colorscale.
        # np.searchsorted returns index at which inserting new value will
        # maintain a sorted array. We use the color to the left of index
        cid = np.searchsorted(plot_dict["Levels"], obs_data[istn])
        cid = max(0, cid - 1)  # filter out zero and max when seeking 'left'
        cid = min(len(plot_dict["Colours"]) - 1, cid)
        pcol = plot_dict["Colours"][cid]

        # Overlay contourf with observations
        plt.plot(
            anet_aod["Lons"][istn],
            anet_aod["Lats"][istn],
            color=pcol,
            marker="o",
            markeredgecolor="k",
            markeredgewidth=2,
            markersize=9,
        )

    # Decorate the plot
    plt.title(plot_dict["Title"], size=24)
    plt.colorbar(cf, orientation="horizontal")
    plt.gca().coastlines(color="#525252")

    # Statistics on plot
    plt.figtext(0.25,
                0.30,
                "MEAN: {0:6.3f}. RMSE v Aeronet: {1:6.3f} ".format(
                    plot_dict["Mean_aod"], plot_dict["RMS_aod"]),
                size=20)


def read_aeronet_clim(aeronet_dir, years):
    """Read in pre-processed AERONET AOD values."""
    """
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
    cubes : Cube list for observational AOD data.
        Individual cubes (with time, lat, lon coordinates) contain a time
        series of pre-processed monthly mean AOD data from AeroNET
        observations at a single station. The cube's attributes holds metadata
        for the Station and for AeroNet.
    aeronet_dir : The file path to the AeroNET observation data files.

    Returns
    -------
    anet_aod : Dictionary for AOD observations data.
        This dictionary includes station names, station latitude, station
        longitude and model-obs comparison statistics and AOD values for a
        specified period i.e. all 12 months, 4 seasons or 1 annual

    """

    # Open observational AOD data files
    try:
        cubes = iris.load(aeronet_dir + "*")

    except Exception:
        raise AeroAnsError(
            'READ_AERONET:Error loading Aeronet data: "{}"'.format(
                aeronet_dir))

    # Lists to hold station data
    stnlist = []
    lats = []
    lons = []
    aod_month = []
    aod_seas = []
    aod_ann = []

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

        start_yr = years[0] - 1
        end_yr = years[1]
        pdt1 = iris.time.PartialDateTime(year=start_yr, month=11, day=1)
        pdt2 = iris.time.PartialDateTime(year=end_yr, month=12, day=1)
        time_constr = iris.Constraint(
            time=lambda cell: pdt1 < cell.point < pdt2)
        cube = cube.extract(time_constr)

        # Calculate multi-annual annual mean
        ann = cube.collapsed("time", iris.analysis.MEAN)
        aod_ann.append(ann.data)

        # Calculate multi-annual seasonal means
        season = cube.aggregated_by(["clim_season"], iris.analysis.MEAN)
        aod_seas.append(season.data)

        # Calculate multi-annual monthly means
        iris.coord_categorisation.add_month(cube, "time", name="month")
        month = cube.aggregated_by(["month"], iris.analysis.MEAN)
        aod_month.append(month.data)

    # Analyse and plot the AOD metrics
    anet_aod = {
        "Stations": stnlist,
        "Lats": lats,
        "Lons": lons,
        "MA_month": aod_month,
        "MA_season": aod_seas,
        "MA_annual": aod_ann,
    }

    return anet_aod


def aod_analyse(MD, aeronet_dir, seasons, years, wavel):
    """Evaluates AOD vs Aeronet, generates plots and returns metric."""
    """
    Parameters
    ----------
    MD : Iris Cube. Iris cube containing a modelled AOD climatology with
         coordinates (time, lat, lon).
    aeronet_dir : String. The file path for the AeroNET observations stored
         as NetCDF files.
    wavel : AOD wavelength, default = 440nm -translates to pseudo-lev.

    Returns
    -------
    metric : Dictionary for model-obs AOD comparison statistics.
        These statistics include root mean square error (RMS).

   """

    # Convert wave length nm -> um
    wv_mi = str(float(wavel) / 1000.0)

    # Dictionary for output metric/s
    metric = dict()

    # Get model run id
    md_id = MD.attributes["parent_source_id"]

    # Add bounds for lat and lon if not present
    MD = add_bounds(MD)

    # Plotting settings
    # Contour levels and colours
    clevs = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 2.0]
    cmapr = mpl_cm.get_cmap("brewer_Spectral_11")
    cmap = colors.ListedColormap(cmapr(range(cmapr.N))[-1::-1])
    colours = cmap.colors

    # Read AOD measurements
    anet_aod = read_aeronet_clim(aeronet_dir, years)

    # Seasonal plots and statistics
    # Co-locate model grid points with measurement sites --func from aero_utils
    aod_at_anet = extract_pt(MD, anet_aod["Lats"], anet_aod["Lons"])

    # Set up the figure
    figures = []

    # Loop over seasons
    for iseas, seas in enumerate(MD.slices_over("season_number")):
        print("Analysing AOD for ", md_id, ":", seasons[iseas])

        # Generate statistics required - area-weighted mean
        grid_areas = iris.analysis.cartography.area_weights(seas)
        global_mean = seas.collapsed(["latitude", "longitude"],
                                     iris.analysis.MEAN,
                                     weights=grid_areas)

        # Locate valid obs data
        masked_obs = [ma.getmask(x[iseas]) for x in anet_aod["MA_season"]]

        # Extract model and obs data for iseas
        seas_anet_obs = [x[iseas] for x in anet_aod["MA_season"]]
        seas_anet_md = [x[iseas] for x in aod_at_anet]

        # Match model data with valid obs data
        valid_obs = [b for a, b in zip(masked_obs, seas_anet_obs) if not a]
        valid_md = [b for a, b in zip(masked_obs, seas_anet_md) if not a]

        # Diff, model mean and RMS
        diff = [a - b for a, b in zip(valid_md, valid_obs)]
        rms_aod = np.sqrt(np.sum([x * x for x in diff]) / len(diff))

        # Populate metric
        metric["Aerosol optical depth at " + wv_mi + " um " +
               seasons[iseas]] = rms_aod

        # Plot this run and season
        figure = plt.figure(figsize=(11, 8), dpi=300)

        title = (md_id + ":\n Total Aerosol Optical Depth at " + wv_mi +
                 " microns " + seasons[iseas])

        # Plot dictionary
        plot_dict = {
            "Mean_aod": global_mean.data,
            "RMS_aod": rms_aod,
            "Levels": clevs,
            "Colours": colours,
            "Title": title,
            "Season": seasons[iseas],
        }
        plot_aod_mod_obs(seas, seas_anet_obs, anet_aod, plot_dict)

        figures.append(figure)

    return metric, figures


def main(config):
    """Produce the AOD climatology metric from ground-based AeroNet
    observations.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    logger = logging.getLogger(__name__)

    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    # Settings
    wavel = "440"
    home_dir = "/home/h04/chardacr/"
    aeronet_dir = home_dir + "aod_climatology/OBS/single_site_nc_limited/"

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)

        for attributes in group:

            logger.info("Loading filename: ", attributes["filename"])

            input_file = attributes["filename"]
            cube = iris.load_cube(input_file)

            logger.info("Cube:", cube)

            # Read in the pre-processed AERONET AOD values
            seasons = ["DJF", "MAM", "JJA", "SON"]
            years = [attributes["start_year"], attributes["end_year"]]
            metric, figures = aod_analyse(cube,
                                          aeronet_dir,
                                          seasons,
                                          years,
                                          wavel=wavel)

            for ifig, fig in enumerate(figures):
                output_file = (model_dataset + "_" + attributes["activity"] +
                               "_" + attributes["mip"] + "_" +
                               attributes["exp"] + "_" +
                               attributes["short_name"] + "_" +
                               attributes["grid"] + "_" +
                               str(attributes["start_year"]) + "_" +
                               str(attributes["end_year"]) + "_" +
                               seasons[ifig])

                output_path = get_plot_filename(output_file, config)
                fig.savefig(output_path, close=True)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
