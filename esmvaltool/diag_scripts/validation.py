"""Validation Diagnostic.

This diagnostic uses two datasets (control and experiment), applies
operations on their data, and plots one against the other. It can
optionally use a number of OBS, OBS4MIPS datasets.
"""

import logging
import os

import iris
import iris.analysis.maths as imath
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import extract_region, extract_season

from esmvaltool.diag_scripts.shared import (
    apply_supermeans,
    get_control_exper_obs,
    group_metadata,
    run_diagnostic,
)
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


def _get_provenance_record(cfg, plot_file, caption, loc):
    """Create a provenance record describing the diagnostic data and plot."""
    ancestor_files = []
    for dataset in cfg["input_data"].values():
        if (
            dataset["alias"] in plot_file
            and dataset["short_name"] in plot_file
        ):
            ancestor_files.append(dataset["filename"])
    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_types": ["map", "metrics"],
        "authors": [
            "predoi_valeriu",
        ],
        "references": [],
        "plot_file": plot_file,
        "ancestors": ancestor_files,
    }

    p_cfg = {}
    p_cfg["run_dir"] = loc
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, record)


def plot_contour(cube, cfg, plt_title, file_name):
    """Plot a contour with iris.quickplot (qplot)."""
    if len(cube.shape) == 2:
        qplt.contourf(cube, cmap="RdYlBu_r")
    else:
        qplt.contourf(cube[0], cmap="RdYlBu_r")
    plt.title(plt_title)
    plt.gca().coastlines()
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()
    _get_provenance_record(
        cfg, file_name, plt_title, loc=os.path.basename(file_name)
    )


def save_plotted_cubes(cube, cfg, plot_name):
    """Save cubes that were plotted to disk."""
    if "save_cubes" in cfg:
        if cfg["save_cubes"]:
            save_name = plot_name.replace("png", "nc")
            save_path = os.path.join(cfg["work_dir"], save_name)
            iris.save(cube, save_path)
            parsed_file_name = plot_name.replace(".png", "").replace("_", " ")
            prov_name = "Output file: " + parsed_file_name
            # files are overwritten once; provenance fails second time
            try:
                _get_provenance_record(
                    cfg, save_path, prov_name, loc=os.path.dirname(save_path)
                )
            except KeyError:
                pass


def plot_latlon_cubes(
    cube_1, cube_2, cfg, data_names, obs_name=None, season=None
):
    """Plot lat-lon vars for control, experiment, and obs.

    Also plot Difference plots (control-exper, control-obs)
    cube_1: first cube (dataset: dat1)
    cube_2: second cube (dataset: dat2)
    cfg: configuration dictionary
    data_names: var + '_' + dat1 + '_vs_' + dat2
    """
    if not season:
        plot_name = "_".join([cfg["analysis_type"], data_names]) + ".png"
        plot_title = "alltime " + cfg["analysis_type"] + ": " + data_names
        plot_file_path = os.path.join(
            cfg["plot_dir"], "alltime", "Difference_" + plot_name
        )
    else:
        plot_name = (
            "_".join([cfg["analysis_type"], data_names, season]) + ".png"
        )
        plot_title = season + " " + cfg["analysis_type"] + ": " + data_names
        plot_file_path = os.path.join(
            cfg["plot_dir"], season, "Difference_" + plot_name
        )
    cubes = [cube_1, cube_2]

    # plot difference: cube_1 - cube_2; use numpy.ma.abs()
    diffed_cube = imath.subtract(cube_1, cube_2)

    plot_contour(diffed_cube, cfg, "Difference " + plot_title, plot_file_path)
    save_plotted_cubes(diffed_cube, cfg, "Difference_" + plot_name)

    # plot each cube
    var = data_names.split("_")[0]
    if not obs_name:
        cube_names = [
            data_names.replace(f"{var}_", "").split("_vs_")[i]
            for i in range(2)
        ]
        for cube, cube_name in zip(cubes, cube_names, strict=True):
            if not season:
                plot_file_path = os.path.join(
                    cfg["plot_dir"],
                    "alltime",
                    "_".join([cube_name, var]) + ".png",
                )
                plot_contour(
                    cube,
                    cfg,
                    " ".join([cube_name, cfg["analysis_type"], var]),
                    plot_file_path,
                )
            else:
                plot_file_path = os.path.join(
                    cfg["plot_dir"],
                    season,
                    "_".join([cube_name, var, season]) + ".png",
                )
                plot_contour(
                    cube,
                    cfg,
                    " ".join([season, cube_name, cfg["analysis_type"], var]),
                    plot_file_path,
                )
            save_plotted_cubes(cube, cfg, os.path.basename(plot_file_path))
    else:
        # obs is always cube_2
        if not season:
            plot_file_path = os.path.join(
                cfg["plot_dir"], "alltime", "_".join([obs_name, var]) + ".png"
            )
            plot_contour(
                cube_2,
                cfg,
                " ".join([obs_name, cfg["analysis_type"], var]),
                plot_file_path,
            )
        else:
            plot_file_path = os.path.join(
                cfg["plot_dir"],
                season,
                "_".join([obs_name, var, season]) + ".png",
            )
            plot_contour(
                cube_2,
                cfg,
                " ".join([season, obs_name, cfg["analysis_type"], var]),
                plot_file_path,
            )
        save_plotted_cubes(cube_2, cfg, os.path.basename(plot_file_path))


def plot_zonal_cubes(cube_1, cube_2, cfg, plot_data):
    """Plot cubes data vs latitude or longitude when zonal meaning."""
    # xcoordinate: latotude or longitude (str)
    data_names, xcoordinate, period = plot_data
    var = data_names.split("_")[0]
    cube_names = data_names.replace(var + "_", "").split("_vs_")
    lat_points = cube_1.coord(xcoordinate).points
    plt.plot(lat_points, cube_1.data, label=cube_names[0])
    plt.plot(lat_points, cube_2.data, label=cube_names[1])
    plt.title(
        f"Annual Climatology of {var}"
        if period == "alltime"
        else f"{period} of {var}"
    )
    if xcoordinate == "latitude":
        axis = plt.gca()
        axis.set_xticks(
            [-60, -30, 0, 30, 60],
            labels=[
                "60\N{DEGREE SIGN} S",
                "30\N{DEGREE SIGN} S",
                "0\N{DEGREE SIGN}",
                "30\N{DEGREE SIGN} N",
                "60\N{DEGREE SIGN} N",
            ],
        )
    elif xcoordinate == "longitude":
        axis = plt.gca()
        axis.set_xticks(
            [0, 60, 120, 180, 240, 300, 360],
            labels=[
                "0\N{DEGREE SIGN} E",
                "60\N{DEGREE SIGN} E",
                "120\N{DEGREE SIGN} E",
                "180\N{DEGREE SIGN} E",
                "240\N{DEGREE SIGN} E",
                "300\N{DEGREE SIGN} E",
                "0\N{DEGREE SIGN} E",
            ],
        )
    plt.xlabel(xcoordinate + " (deg)")
    plt.ylabel(f"{var} [{str(cube_1.units)}]")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    png_name = f"{xcoordinate}_{period}_{data_names}.png"
    if xcoordinate == "latitude":
        png_name = "Zonal_Mean_" + png_name
    elif xcoordinate == "longitude":
        png_name = "Merid_Mean_" + png_name
    plot_file_path = os.path.join(cfg["plot_dir"], period, png_name)
    plt.savefig(plot_file_path)
    save_plotted_cubes(
        cube_1,
        cfg,
        "_".join([cube_names[0], os.path.basename(plot_file_path)]),
    )
    save_plotted_cubes(
        cube_2,
        cfg,
        "_".join([cube_names[1], os.path.basename(plot_file_path)]),
    )
    plt.close()
    caption = period + " Zonal/Meridional Mean for " + var + " " + data_names
    _get_provenance_record(
        cfg, plot_file_path, caption, loc=os.path.join(cfg["plot_dir"], period)
    )


def apply_seasons(data_set_dict):
    """Extract seaons and apply a time mean per season."""
    data_file = data_set_dict["filename"]
    logger.info("Loading %s for seasonal extraction", data_file)
    data_cube = iris.load_cube(data_file)
    seasons = ["DJF", "MAM", "JJA", "SON"]
    iris.coord_categorisation.add_season(
        data_cube, "time", "clim_season", seasons
    )
    iris.coord_categorisation.add_season_year(
        data_cube, "time", "season_year", seasons
    )
    season_cubes = [extract_season(data_cube, season) for season in seasons]
    season_meaned_cubes = [
        season_cube.collapsed("time", iris.analysis.MEAN)
        for season_cube in season_cubes
    ]

    return season_meaned_cubes


def coordinate_collapse(data_set, cfg):
    """Perform coordinate-specific collapse and (if) area slicing and mask."""
    # see what analysis needs performing
    analysis_type = cfg["analysis_type"]

    # if subset on LAT-LON
    if "lat_lon_slice" in cfg:
        start_longitude = cfg["lat_lon_slice"]["start_longitude"]
        end_longitude = cfg["lat_lon_slice"]["end_longitude"]
        start_latitude = cfg["lat_lon_slice"]["start_latitude"]
        end_latitude = cfg["lat_lon_slice"]["end_latitude"]
        data_set = extract_region(
            data_set,
            start_longitude,
            end_longitude,
            start_latitude,
            end_latitude,
        )

    # if apply mask
    if "2d_mask" in cfg:
        mask_file = os.path.join(cfg["2d_mask"])
        mask_cube = iris.load_cube(mask_file)
        if "mask_threshold" in cfg:
            thr = cfg["mask_threshold"]
            data_set.data = np.ma.masked_array(
                data_set.data, mask=mask_cube.data > thr
            )
        else:
            logger.warning("Could not find masking threshold")
            logger.warning("Please specify it if needed")
            logger.warning("Masking on 0-values = True (masked value)")
            data_set.data = np.ma.masked_array(
                data_set.data, mask=mask_cube.data == 0
            )

    # if zonal mean on LON
    if analysis_type == "zonal_mean":
        data_set = data_set.collapsed("longitude", iris.analysis.MEAN)

    # if zonal mean on LAT
    if analysis_type == "meridional_mean":
        data_set = data_set.collapsed("latitude", iris.analysis.MEAN)

    # if vertical mean
    elif analysis_type == "vertical_mean":
        data_set = data_set.collapsed("pressure", iris.analysis.MEAN)

    return data_set


def do_preamble(cfg):
    """Execute some preamble functionality."""
    # prepare output dirs
    time_chunks = ["alltime", "DJF", "MAM", "JJA", "SON"]
    time_plot_dirs = [
        os.path.join(cfg["plot_dir"], t_dir) for t_dir in time_chunks
    ]
    for time_plot_dir in time_plot_dirs:
        if not os.path.exists(time_plot_dir):
            os.makedirs(time_plot_dir)

    # get data
    input_data = cfg["input_data"].values()
    grouped_input_data = group_metadata(
        input_data, "short_name", sort="dataset"
    )

    return input_data, grouped_input_data


def plot_ctrl_exper(ctrl, exper, cfg, plot_key):
    """Call plotting functions and make plots depending on case."""
    if cfg["analysis_type"] == "lat_lon":
        plot_latlon_cubes(ctrl, exper, cfg, plot_key)
    elif cfg["analysis_type"] == "zonal_mean":
        plot_info = [plot_key, "latitude", "alltime"]
        plot_zonal_cubes(ctrl, exper, cfg, plot_info)
    elif cfg["analysis_type"] == "meridional_mean":
        plot_info = [plot_key, "longitude", "alltime"]
        plot_zonal_cubes(ctrl, exper, cfg, plot_info)


def plot_ctrl_exper_seasons(ctrl_seasons, exper_seasons, cfg, plot_key):
    """Call plotting functions and make plots with seasons."""
    seasons = ["DJF", "MAM", "JJA", "SON"]
    if cfg["analysis_type"] == "lat_lon":
        for control_season, experiment_season, season in zip(
            ctrl_seasons, exper_seasons, seasons, strict=True
        ):
            plot_latlon_cubes(
                control_season,
                experiment_season,
                cfg,
                plot_key,
                obs_name=None,
                season=season,
            )
    elif cfg["analysis_type"] == "zonal_mean":
        for control_season, experiment_season, season in zip(
            ctrl_seasons, exper_seasons, seasons, strict=True
        ):
            plot_info = [plot_key, "latitude", season]
            plot_zonal_cubes(control_season, experiment_season, cfg, plot_info)
    elif cfg["analysis_type"] == "meridional_mean":
        for control_season, experiment_season, season in zip(
            ctrl_seasons, exper_seasons, seasons, strict=True
        ):
            plot_info = [plot_key, "longitude", season]
            plot_zonal_cubes(control_season, experiment_season, cfg, plot_info)


def main(cfg):
    """Execute validation analysis and plotting."""
    logger.setLevel(cfg["log_level"].upper())
    input_data, grouped_input_data = do_preamble(cfg)

    # select variables and their corresponding obs files
    for short_name in grouped_input_data:
        logger.info("Processing variable %s", short_name)

        # get the control, experiment and obs dicts
        cmip_type = cfg["cmip_era"] if "cmip_era" in cfg else None
        ctrl, exper, obs = get_control_exper_obs(
            short_name, input_data, cfg, cmip_type=cmip_type
        )
        # set a plot key holding info on var and data set names
        plot_key = f"{short_name}_{ctrl['alias']}_vs_{exper['alias']}"
        control_dataset_name = ctrl["alias"]

        # get seasons if needed then apply analysis
        if cfg["seasonal_analysis"]:
            ctrl_seasons = apply_seasons(ctrl)
            exper_seasons = apply_seasons(exper)
            ctrl_seasons = [
                coordinate_collapse(cts, cfg) for cts in ctrl_seasons
            ]
            exper_seasons = [
                coordinate_collapse(exps, cfg) for exps in exper_seasons
            ]
            plot_ctrl_exper_seasons(ctrl_seasons, exper_seasons, cfg, plot_key)
            if obs:
                for iobs in obs:
                    obs_seasons = apply_seasons(iobs)
                    obs_seasons = [
                        coordinate_collapse(obss, cfg) for obss in obs_seasons
                    ]
                    plot_key_obs = (
                        f"{short_name}_{ctrl['alias']}"
                        + f"_vs_{iobs['alias']}"
                    )
                    plot_ctrl_exper_seasons(
                        ctrl_seasons, obs_seasons, cfg, plot_key_obs
                    )

        # apply the supermeans (MEAN on time), collapse a coord and plot
        ctrl, exper, obs_list = apply_supermeans(ctrl, exper, obs)
        ctrl = coordinate_collapse(ctrl, cfg)
        exper = coordinate_collapse(exper, cfg)
        plot_ctrl_exper(ctrl, exper, cfg, plot_key)

        # apply desired analysis on obs's
        if obs_list:
            for obs_i, obsfile in zip(obs_list, obs, strict=True):
                obs_analyzed = coordinate_collapse(obs_i, cfg)
                obs_name = obsfile["alias"]
                plot_key = f"{short_name}_{control_dataset_name}_vs_{obs_name}"
                if cfg["analysis_type"] == "lat_lon":
                    plot_latlon_cubes(
                        ctrl, obs_analyzed, cfg, plot_key, obs_name=obs_name
                    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
