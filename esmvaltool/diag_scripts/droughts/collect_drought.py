"""Compares SPI/SPEI data from models with observations/reanalysis.

Description
-----------
This diagnostic applies drought characteristics based on Martin (2018)
to data produced by spei.R. These characteristics and differences to
a reference dataset or between different time periods are plotted for each
dataset and multi-model mean. The calculated frequency, duration and severity
of drought events are saved to netcdf files for further use.
It expects multiple datasets for a particular index as input. The reference
dataset can be specified with ``reference_dataset`` and is not part of the
multi-model mean.

.. note:: Previous Version:
   With ESMValTool v2.12 and previous, multiple collect_drought_*.py
   diagnostics and a collect_drought_func.py diagnostic existed in the
   `droughtindex` folders. Those have been archived and replaced by this
   diagnostic in v2.13.

Configuration options
---------------------
indexname: str
    The indexname is used to generate filenames, plot titles and captions.
    Should be ``SPI`` or ``SPEI``.
reference_dataset: str
    Dataset name to use for comparison (excluded from MMM). With
    ``compare_intervals=True`` this option has no effect.
threshold: float, optional (default: -2.0)
    Threshold for an event to be considered as drought.
compare_intervals: bool, false
    If true, start and end of the time periods are compared instead of
    models and reference. The lengths of start and end period is given by
    ``comparison_period``.
comparison_period: int
    Number of years from start and end of the full period to be compared.
    Should be < (end_year - start_year)/2.
    If ``compare_intervals=False`` this option has no effect.
plot_models: bool, false
    Save plots for each individual model, in addition to multi-model mean.
start_year: int
    This option is used to select the time slices for comparison if
    ``compare_intervals=True``.
end_year: int
    This option is used to select the time slices for comparison if
    ``compare_intervals=True``.
"""

import datetime as dt
import logging
from pathlib import Path
from pprint import pformat

import cartopy.crs as cart
import iris
import matplotlib.pyplot as plt
import numpy as np
from iris.analysis import Aggregator

from esmvaltool.diag_scripts.droughts.utils import (
    count_spells,
    create_cube_from_data,
)
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

log = logging.getLogger(Path(__file__).name)


def _get_provenance_record(
    ancestor_files,
    caption,
    domains,
    refs,
    plot_type="geo",
) -> dict:
    """Get Provenance record."""
    return {
        "caption": caption,
        "statistics": ["mean"],
        "domains": domains,
        "plot_type": plot_type,
        "themes": ["phys"],
        "authors": [
            "weigel_katja",
            "adeniyi_kemisola",
        ],
        "references": refs,
        "ancestors": ancestor_files,
    }


def _get_drought_data(cfg, cube):
    """Prepare data and calculate characteristics."""
    # make a new cube to increase the size of the data array
    # Make an aggregator from the user function.
    spell_no = Aggregator(
        "spell_count",
        count_spells,
        units_func=lambda units: 1,
    )
    new_cube = _make_new_cube(cube)

    # calculate the number of drought events and their average duration
    drought_show = new_cube.collapsed(
        "time",
        spell_no,
        threshold=cfg["threshold"],
    )
    drought_show.rename("Drought characteristics")
    # length of time series
    time_length = len(new_cube.coord("time").points) / 12.0
    # Convert number of droughtevents to frequency (per year)
    drought_show.data[:, :, 0] = drought_show.data[:, :, 0] / time_length
    return drought_show


def _provenance_map_spei(cfg, name_dict, spei, dataset_name):
    """Set provenance for plot_map_spei."""
    caption = (
        "Global map of "
        + name_dict["drought_char"]
        + " ["
        + name_dict["unit"]
        + "] "
        + "based on "
        + cfg["indexname"]
        + "."
    )

    if cfg["indexname"].lower == "spei":
        set_refs = ["martin18grl", "vicente10jclim"]
    elif cfg["indexname"].lower == "spi":
        set_refs = ["martin18grl", "mckee93proc"]
    else:
        set_refs = ["martin18grl"]

    provenance_record = _get_provenance_record(
        [name_dict["input_filenames"]],
        caption,
        ["global"],
        set_refs,
    )

    diagnostic_file = get_diagnostic_filename(
        cfg["indexname"]
        + "_map"
        + name_dict["add_to_filename"]
        + "_"
        + dataset_name,
        cfg,
    )
    plot_file = get_plot_filename(
        cfg["indexname"]
        + "_map"
        + name_dict["add_to_filename"]
        + "_"
        + dataset_name,
        cfg,
    )
    log.info("Saving analysis results to %s", diagnostic_file)
    cubesave = create_cube_from_data(spei, name_dict)
    iris.save(cubesave, target=diagnostic_file)
    log.info(
        "Recording provenance of %s:\n%s",
        diagnostic_file,
        pformat(provenance_record),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, provenance_record)
        provenance_logger.log(diagnostic_file, provenance_record)


def _provenance_map_spei_multi(cfg, data_dict, spei, input_filenames):
    """Set provenance for plot_map_spei_multi."""
    caption = (
        f"Global map of the multi-model mean of "
        f"{data_dict['drought_char']} [{data_dict['unit']}] based on "
        f"{cfg['indexname']}."
    )
    if cfg["indexname"].lower == "spei":
        set_refs = ["martin18grl", "vicente10jclim"]
    elif cfg["indexname"].lower == "spi":
        set_refs = ["martin18grl", "mckee93proc"]
    else:
        set_refs = ["martin18grl"]

    provenance_record = _get_provenance_record(
        input_filenames,
        caption,
        ["global"],
        set_refs,
    )

    diagnostic_file = get_diagnostic_filename(
        cfg["indexname"]
        + "_map"
        + data_dict["filename"]
        + "_"
        + data_dict["datasetname"],
        cfg,
    )
    plot_file = get_plot_filename(
        cfg["indexname"]
        + "_map"
        + data_dict["filename"]
        + "_"
        + data_dict["datasetname"],
        cfg,
    )
    log.info("Saving analysis results to %s", diagnostic_file)
    iris.save(create_cube_from_data(spei, data_dict), target=diagnostic_file)
    log.info(
        "Recording provenance of %s:\n%s",
        diagnostic_file,
        pformat(provenance_record),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, provenance_record)
        provenance_logger.log(diagnostic_file, provenance_record)


def _plot_multi_model_maps(
    cfg,
    all_drought_mean,
    lats_lons,
    input_filenames,
    tstype,
):
    """Prepare plots for multi-model mean."""
    data_dict = {
        "latitude": lats_lons[0],
        "longitude": lats_lons[1],
        "model_kind": tstype,
    }
    if tstype == "Difference":
        # RCP85 Percentage difference
        data_dict.update(
            {
                "data": all_drought_mean[:, :, 0],
                "var": "diffnumber",
                "datasetname": "Percentage",
                "drought_char": "Number of drought events",
                "unit": "%",
                "filename": "Percentage_difference_of_No_of_Events",
                "drought_numbers_level": np.arange(-100, 110, 10),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="rainbow",
        )

        data_dict.update(
            {
                "data": all_drought_mean[:, :, 1],
                "var": "diffduration",
                "drought_char": "Duration of drought events",
                "filename": "Percentage_difference_of_Dur_of_Events",
                "drought_numbers_level": np.arange(-100, 110, 10),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="rainbow",
        )

        data_dict.update(
            {
                "data": all_drought_mean[:, :, 2],
                "var": "diffseverity",
                "drought_char": "Severity Index of drought events",
                "filename": "Percentage_difference_of_Sev_of_Events",
                "drought_numbers_level": np.arange(-50, 60, 10),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="rainbow",
        )

        data_dict.update(
            {
                "data": all_drought_mean[:, :, 3],
                "var": "diff" + (cfg["indexname"]).lower(),
                "drought_char": "Average "
                + cfg["indexname"]
                + " of drought events",
                "filename": "Percentage_difference_of_Avr_of_Events",
                "drought_numbers_level": np.arange(-50, 60, 10),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="rainbow",
        )
    else:
        data_dict.update(
            {
                "data": all_drought_mean[:, :, 0],
                "var": "frequency",
                "unit": "year-1",
                "drought_char": "Number of drought events per year",
                "filename": tstype + "_No_of_Events_per_year",
                "drought_numbers_level": np.arange(0, 0.4, 0.05),
            }
        )
        if tstype == "Observations":
            data_dict["datasetname"] = "Mean"
        else:
            data_dict["datasetname"] = "MultiModelMean"
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="gnuplot",
        )

        data_dict.update(
            {
                "data": all_drought_mean[:, :, 1],
                "var": "duration",
                "unit": "month",
                "drought_char": "Duration of drought events [month]",
                "filename": tstype + "_Dur_of_Events",
                "drought_numbers_level": np.arange(0, 6, 1),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="gnuplot",
        )

        data_dict.update(
            {
                "data": all_drought_mean[:, :, 2],
                "var": "severity",
                "unit": "1",
                "drought_char": "Severity Index of drought events",
                "filename": tstype + "_Sev_index_of_Events",
                "drought_numbers_level": np.arange(0, 9, 1),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="gnuplot",
        )
        namehlp = "Average " + cfg["indexname"] + " of drought events"
        namehlp2 = tstype + "_Average_" + cfg["indexname"] + "_of_Events"
        data_dict.update(
            {
                "data": all_drought_mean[:, :, 3],
                "var": (cfg["indexname"]).lower(),
                "unit": "1",
                "drought_char": namehlp,
                "filename": namehlp2,
                "drought_numbers_level": np.arange(-2.8, -1.8, 0.2),
            }
        )
        plot_map_spei_multi(
            cfg,
            data_dict,
            input_filenames,
            colormap="gnuplot",
        )


def _plot_single_maps(cfg, cube2, drought_show, tstype, input_filenames):
    """Plot map of drought characteristics for individual models and times."""
    cube2.data = drought_show.data[:, :, 0]
    name_dict = {
        "add_to_filename": tstype + "_No_of_Events_per_year",
        "name": tstype + " Number of drought events per year",
        "var": "frequency",
        "unit": "year-1",
        "drought_char": "Number of drought events per year",
        "input_filenames": input_filenames,
    }
    plot_map_spei(cfg, cube2, np.arange(0, 0.4, 0.05), name_dict)
    # plot the average duration of drought events
    cube2.data = drought_show.data[:, :, 1]
    name_dict.update(
        {
            "add_to_filename": tstype + "_Dur_of_Events",
            "name": tstype + " Duration of drought events(month)",
            "var": "duration",
            "unit": "month",
            "drought_char": "Number of drought events per year",
            "input_filenames": input_filenames,
        }
    )
    plot_map_spei(cfg, cube2, np.arange(0, 6, 1), name_dict)
    # plot the average severity index of drought events
    cube2.data = drought_show.data[:, :, 2]
    name_dict.update(
        {
            "add_to_filename": tstype + "_Sev_index_of_Events",
            "name": tstype + " Severity Index of drought events",
            "var": "severity",
            "unit": "1",
            "drought_char": "Number of drought events per year",
            "input_filenames": input_filenames,
        }
    )
    plot_map_spei(cfg, cube2, np.arange(0, 9, 1), name_dict)
    # plot the average spei of drought events
    cube2.data = drought_show.data[:, :, 3]
    namehlp = tstype + "_Avr_" + cfg["indexname"] + "_of_Events"
    namehlp2 = tstype + "_Average_" + cfg["indexname"] + "_of_Events"
    name_dict.update(
        {
            "add_to_filename": namehlp,
            "name": namehlp2,
            "var": "severity",
            "unit": "1",
            "drought_char": "Number of drought events per year",
            "input_filenames": input_filenames,
        }
    )
    plot_map_spei(cfg, cube2, np.arange(-2.8, -1.8, 0.2), name_dict)


def plot_map_spei_multi(
    cfg,
    data_dict,
    input_filenames,
    colormap="jet",
) -> None:
    """Plot contour maps for multi model mean."""
    spei = np.ma.array(data_dict["data"], mask=np.isnan(data_dict["data"]))
    # Get latitudes and longitudes from cube
    lons = data_dict["longitude"]
    if max(lons) > 180.0:
        lons = np.where(lons > 180, lons - 360, lons)
        # sort the array
        index = np.argsort(lons)
        lons = lons[index]
        spei = spei[np.ix_(range(data_dict["latitude"].size), index)]
    # Plot data
    # Create figure and axes instances
    subplot_kw = {"projection": cart.PlateCarree(central_longitude=0.0)}
    fig, axx = plt.subplots(figsize=(6.5, 4), subplot_kw=subplot_kw)
    axx.set_extent(
        [-180.0, 180.0, -90.0, 90.0],
        cart.PlateCarree(central_longitude=0.0),
    )
    # Draw filled contours
    cnplot = plt.contourf(
        lons,
        data_dict["latitude"],
        spei,
        data_dict["drought_numbers_level"],
        transform=cart.PlateCarree(central_longitude=0.0),
        cmap=colormap,
        extend="both",
        corner_mask=False,
    )
    # Style plot
    axx.coastlines()
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation="horizontal")
    if data_dict["model_kind"] == "Difference":
        cbar.set_label(
            data_dict["model_kind"] + " " + data_dict["drought_char"] + " [%]",
        )
    else:
        cbar.set_label(
            data_dict["model_kind"] + " " + data_dict["drought_char"],
        )
    axx.set_xlabel("Longitude")
    axx.set_ylabel("Latitude")
    axx.set_title(
        f"{data_dict['datasetname']} {data_dict['model_kind']} "
        f"{data_dict['drought_char']}",
    )
    # set ticks
    axx.set_xticks(np.linspace(-180, 180, 7))
    axx.set_xticklabels(
        [
            "180°W",
            "120°W",
            "60°W",
            "0°",
            "60°E",
            "120°E",
            "180°E",
        ]
    )
    axx.set_yticks(np.linspace(-90, 90, 7))
    axx.set_yticklabels(["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])

    fig.tight_layout()
    fig.savefig(
        get_plot_filename(
            cfg["indexname"]
            + "_map"
            + data_dict["filename"]
            + "_"
            + data_dict["datasetname"],
            cfg,
        ),
        dpi=300,
    )
    plt.close()
    _provenance_map_spei_multi(cfg, data_dict, spei, input_filenames)


def plot_map_spei(cfg, cube, levels, name_dict) -> None:
    """Plot contour map."""
    mask = np.isnan(cube.data)
    spei = np.ma.array(cube.data, mask=mask)
    np.ma.masked_less_equal(spei, 0)
    # Get latitudes and longitudes from cube
    name_dict.update({"latitude": cube.coord("latitude").points})
    lons = cube.coord("longitude").points
    lons = np.where(lons > 180, lons - 360, lons)
    # sort the array
    index = np.argsort(lons)
    lons = lons[index]
    name_dict.update({"longitude": lons})
    spei = spei[np.ix_(range(len(cube.coord("latitude").points)), index)]
    # Get data set name from cube
    try:
        dataset_name = cube.metadata.attributes["model_id"]
    except KeyError:
        try:
            dataset_name = cube.metadata.attributes["source_id"]
        except KeyError:
            dataset_name = "Observations"
    # Plot data
    # Create figure and axes instances
    subplot_kw = {"projection": cart.PlateCarree(central_longitude=0.0)}
    fig, axx = plt.subplots(figsize=(8, 4), subplot_kw=subplot_kw)
    axx.set_extent(
        [-180.0, 180.0, -90.0, 90.0],
        cart.PlateCarree(central_longitude=0.0),
    )
    # Draw filled contours
    cnplot = plt.contourf(
        lons,
        cube.coord("latitude").points,
        spei,
        levels,
        transform=cart.PlateCarree(central_longitude=0.0),
        cmap="gnuplot",
        extend="both",
        corner_mask=False,
    )
    axx.coastlines()
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation="horizontal")
    cbar.set_label(name_dict["name"])
    axx.set_xlabel("Longitude")
    axx.set_ylabel("Latitude")
    axx.set_title(dataset_name + " " + name_dict["name"])

    # Set up x and y ticks
    axx.set_xticks(np.linspace(-180, 180, 7))
    axx.set_xticklabels(
        ["180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"],
    )
    axx.set_yticks(np.linspace(-90, 90, 7))
    axx.set_yticklabels(["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])
    fig.tight_layout()
    basename = (
        f"{cfg['indexname']}_map_{name_dict['add_to_filename']}_{dataset_name}"
    )
    fig.savefig(get_plot_filename(basename, cfg), dpi=300)
    plt.close()
    _provenance_map_spei(cfg, name_dict, spei, dataset_name)


def _plot_models_vs_obs(cfg, cube, mmm, obs, fnames):
    """Compare drought metrics of multi-model mean to observations."""
    latslons = [cube.coord(i).points for i in ["latitude", "longitude"]]
    perc_diff = (obs - mmm) / (obs + mmm) * 200
    _plot_multi_model_maps(cfg, mmm, latslons, fnames, "Historic")
    _plot_multi_model_maps(cfg, obs, latslons, fnames, "Observations")
    _plot_multi_model_maps(cfg, perc_diff, latslons, fnames, "Difference")


def _plot_future_vs_past(cfg, cube, slices, fnames):
    """Compare drought metrics of future and historic time slices."""
    latslons = [cube.coord(i).points for i in ["latitude", "longitude"]]
    slices["Difference"] = (
        (slices["Future"] - slices["Historic"])
        / (slices["Future"] + slices["Historic"])
        * 200
    )
    for tstype in ["Historic", "Future", "Difference"]:
        _plot_multi_model_maps(cfg, slices[tstype], latslons, fnames, tstype)


def _make_new_cube(cube):
    """Make a new cube with an extra dimension for result of spell count."""
    new_shape = (*cube.shape, 4)
    new_data = iris.util.broadcast_to_shape(cube.data, new_shape, [0, 1, 2])
    new_cube = iris.cube.Cube(new_data)
    new_cube.add_dim_coord(
        iris.coords.DimCoord(cube.coord("time").points, long_name="time"),
        0,
    )
    new_cube.add_dim_coord(
        iris.coords.DimCoord(
            cube.coord("latitude").points,
            long_name="latitude",
        ),
        1,
    )
    new_cube.add_dim_coord(
        iris.coords.DimCoord(
            cube.coord("longitude").points,
            long_name="longitude",
        ),
        2,
    )
    new_cube.add_dim_coord(
        iris.coords.DimCoord([0, 1, 2, 3], long_name="z"),
        3,
    )
    return new_cube


def _set_tscube(cfg, cube, time, tstype):
    """Time slice from a cube with start/end given by cfg."""
    if tstype == "Future":
        start_year = cfg["end_year"] - cfg["comparison_period"] + 1
        start = dt.datetime(start_year, 1, 15, 0, 0, 0)
        end = dt.datetime(cfg["end_year"], 12, 16, 0, 0, 0)
    elif tstype == "Historic":
        start = dt.datetime(cfg["start_year"], 1, 15, 0, 0, 0)
        end_year = cfg["start_year"] + cfg["comparison_period"] - 1
        end = dt.datetime(end_year, 12, 16, 0, 0, 0)
    else:
        raise ValueError("Unknown time slice type: " + tstype)
    stime = time.nearest_neighbour_index(time.units.date2num(start))
    etime = time.nearest_neighbour_index(time.units.date2num(end))
    return cube[stime:etime, :, :]


def main(cfg) -> None:
    """Run the diagnostic."""
    # Read input data
    drought_data = []
    drought_slices = {"Historic": [], "Future": []}
    fnames = []  # why do we need them?
    ref_data = None
    for meta in cfg["input_data"].values():
        fname = meta["filename"]
        cube = iris.load_cube(fname)
        fnames.append(fname)
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()
        cube_mean = cube.collapsed("time", iris.analysis.MEAN)
        if cfg.get("compare_intervals", False):
            # calculate and plot metrics per time slice
            for tstype in ["Historic", "Future"]:
                ts_cube = _set_tscube(cfg, cube, cube.coord("time"), tstype)
                drought_show = _get_drought_data(cfg, ts_cube)
                drought_slices[tstype].append(drought_show.data)
                if cfg.get("plot_models", False):
                    _plot_single_maps(
                        cfg,
                        cube_mean,
                        drought_show,
                        tstype,
                        fname,
                    )
        else:
            # calculate and plot metrics per dataset
            drought_show = _get_drought_data(cfg, cube)
            if meta["dataset"] == cfg["reference_dataset"]:
                ref_data = drought_show.data
            else:
                drought_data.append(drought_show.data)
            if cfg.get("plot_models", False):
                _plot_single_maps(
                    cfg,
                    cube_mean,
                    drought_show,
                    "Historic",
                    fname,
                )

    if cfg.get("compare_intervals", False):
        # calculate multi model mean for time slices
        slices = {k: np.array(v) for k, v in drought_slices.items()}
        mean_slices = {k: np.nanmean(v, axis=0) for k, v in slices.items()}
        _plot_future_vs_past(cfg, cube, mean_slices, fnames)
    else:
        # calculate multi model mean and compare with reference dataset
        drought_data = np.array(drought_data)
        mmm = np.nanmean(np.array(drought_data), axis=0)
        _plot_models_vs_obs(cfg, cube, mmm, ref_data, fnames)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
