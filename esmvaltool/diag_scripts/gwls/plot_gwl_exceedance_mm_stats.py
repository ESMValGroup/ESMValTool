"""Plot multimodel variable stats at Global Warming Level exceedance years."""

import logging
import os
from pathlib import Path
from pprint import pformat

import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from esmvalcore.preprocessor import extract_time
from matplotlib import colors

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    io,
    run_diagnostic,
    select_metadata,
)

SEP = "_"

logger = logging.getLogger(Path(__file__).stem)


def log_provenance(provenance, filename, cfg):
    """Create a provenance record for the output file."""
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)


def calculate_gwl_mm_cube(
    project_data,
    gwl_subset_df,
    window_size,
):
    """Calculate multimodel stats around GWL exceedance years."""
    cubes = iris.cube.CubeList()

    datasets = gwl_subset_df["Model"].unique()
    for dataset in datasets:
        # get every exp for which the data is available for
        all_exps_in_dataset = select_metadata(project_data, dataset=dataset)
        for el in all_exps_in_dataset:
            cube = iris.load_cube(el["filename"])
            logger.info("Cube var name is %s", cube.var_name)
            if cube.var_name == "tas":
                cube.convert_units("Celsius")
            elif cube.var_name == "pr":
                cube = cube * 86400
                cube.units = "mm/day"

            ih.prepare_cube_for_merging(cube, el["filename"])
            # extract window period
            year_of_exceedance = gwl_subset_df[
                (gwl_subset_df["Model"] == dataset)
                & (gwl_subset_df["Exp"] == el["exp"])
            ]["Exceedance_Year"].values[0]
            if np.isnan(year_of_exceedance):
                continue
            start_year = int(year_of_exceedance - ((window_size - 1) / 2))
            end_year = int(year_of_exceedance + ((window_size - 1) / 2))
            logger.info(
                "Model: %s, Exp : %s start year: %s, endyear: %s",
                dataset,
                el["exp"],
                str(start_year),
                str(end_year),
            )
            # Start year Jan 1st (1, 1) and end_year Dec 31st (12, 31)
            cube = extract_time(cube, start_year, 1, 1, end_year, 12, 31)
            cube = cube.collapsed("time", iris.analysis.MEAN)
            cubes.append(cube)
    return cubes


def calculate_mm_stats(cubes, mean_file, stdev_file):
    """Calculate mean and standard deviation from merged mm cubes."""
    # find index of time coord so this can be unified across data sets
    # before merging as the time data points are different across models/scenarios
    index = 0
    for coord in cubes[0].aux_coords:
        if coord.standard_name == "time":
            break
        index = index + 1

    for c in cubes[1:]:
        c.remove_coord("time")
        c.add_aux_coord(cubes[0].aux_coords[index])

    if not cubes:
        logger.info("No model instances exceed this GWL.")
    elif len(cubes) == 1:
        iris.save(cubes[0], mean_file)
        logger.info("No standard deviation calculated for a single instance.")
    else:
        try:
            mm_cube = cubes.merge_cube()
        except TypeError:
            logger.debug(pformat(cubes))
            raise

        mm_mean_cube = mm_cube.collapsed(["cube_label"], iris.analysis.MEAN)
        mm_stdev_cube = mm_cube.collapsed(
            ["cube_label"], iris.analysis.STD_DEV
        )
        iris.save(mm_mean_cube, mean_file)
        iris.save(mm_stdev_cube, stdev_file)


def plot_mean_stats(cfg, project, mean_cube_file, gwl):
    """Plot the multimodel mean."""
    cmap_mean_str = cfg["quickplot"]["cmap_mean"]
    mean_level_params = cfg["quickplot"]["mean_level_params"]
    filename = SEP.join([project, "mm_mean", str(gwl)]) + ".png"
    mean_plot_file = os.path.join(cfg["plot_dir"], filename)
    cube = iris.load_cube(mean_cube_file)
    plt.axes(projection=ccrs.Robinson())
    levels = np.arange(
        float(mean_level_params[0]),
        float(mean_level_params[1]),
        float(mean_level_params[2]),
    )
    if cfg["quickplot"]["title_var"] == "Temperature":
        qplt.contourf(
            cube,
            norm=colors.CenteredNorm(),
            cmap=cmap_mean_str,
            levels=levels,
            extend="max",
        )

    if cfg["quickplot"]["title_var"] == "Precipitation":
        qplt.contourf(
            cube,
            cmap=cmap_mean_str,
            levels=levels,
            extend="max",
        )
    plt.gca().coastlines()
    title_str = f"Multimodel mean of {cfg['quickplot']['title_var']} at {gwl} $^\\circ$ C"
    plt.title(title_str)
    plt.tight_layout()
    plt.savefig(mean_plot_file)
    plt.close()


def plot_stdev_stats(cfg, project, stdev_cube_file, gwl):
    """Plot multimodel stdev."""
    cmap_stdev_str = cfg["quickplot"]["cmap_stdev"]
    stdev_level_params = cfg["quickplot"]["stdev_level_params"]
    filename = SEP.join([project, "mm_stdev", str(gwl)]) + ".png"
    stdev_plot_file = os.path.join(cfg["plot_dir"], filename)

    cube = iris.load_cube(stdev_cube_file)
    plt.axes(projection=ccrs.Robinson())
    levels = np.arange(
        float(stdev_level_params[0]),
        float(stdev_level_params[1]),
        float(stdev_level_params[2]),
    )
    qplt.contourf(cube, cmap=cmap_stdev_str, levels=levels, extend="max")
    title_str = f"Multimodel standard deviation of {cfg['quickplot']['title_var']} at {gwl} $^\\circ$ C"

    plt.gca().coastlines()
    plt.title(title_str)
    plt.tight_layout()
    plt.savefig(stdev_plot_file)
    plt.close()


def main(cfg):
    """Execute GWL statistics and plot."""
    gwl_filename = io.get_ancestor_file(cfg, cfg["pattern"])
    logger.info("GWL exceedance years file is %s", gwl_filename)
    gwl_df = pd.read_csv(gwl_filename)
    projects = []
    for data in cfg["input_data"].values():
        # select by by project
        if data["project"] not in projects:
            projects.append(data["project"])

    logger.info("List of Projects: %s", projects)
    for project in projects:
        # get all data sets for project across experiments
        project_data = select_metadata(
            cfg["input_data"].values(), project=project
        )
        for gwl in cfg["gwls"]:
            gwl_subset_df = gwl_df[gwl_df["GWL"] == gwl]
            if gwl_subset_df.shape[0] > 0:
                logger.info(
                    "Calculating means and standard deviations for GWL: %s",
                    str(gwl),
                )
                filename = SEP.join([project, "mm_mean", str(gwl)]) + ".nc"
                mean_file = os.path.join(cfg["work_dir"], filename)
                filename = SEP.join([project, "mm_stdev", str(gwl)]) + ".nc"
                stdev_file = os.path.join(cfg["work_dir"], filename)

                cubes = calculate_gwl_mm_cube(
                    project_data,
                    gwl_subset_df,
                    cfg["window_size"],
                )
                calculate_mm_stats(cubes, mean_file, stdev_file)

                ancestors_file = [gwl_filename]
                ancestors_file = ancestors_file + [
                    d["filename"] for d in project_data
                ]

                if os.path.isfile(mean_file):
                    plot_mean_stats(cfg, project, mean_file, gwl)
                    # write provenance for the NetCDF and png files
                    provenance_dict = {
                        "caption": f"Multimodel mean of {cfg['quickplot']['title_var']} at {gwl} $^\\circ$ C",
                        "ancestors": ancestors_file,
                        "authors": ["swaminathan_ranjini"],
                        "references": ["swaminathan22jclim"],
                        "projects": ["ukesm"],
                        "domains": ["global"],
                    }
                    log_provenance(provenance_dict, mean_file, cfg)
                    provenance_dict.update({"plot_type": ["map"]})
                    filename = (
                        SEP.join([project, "mm_mean", str(gwl)]) + ".png"
                    )
                    filename = os.path.join(cfg["plot_dir"], filename)
                    log_provenance(provenance_dict, filename, cfg)

                if os.path.isfile(stdev_file):
                    plot_stdev_stats(cfg, project, stdev_file, gwl)
                    provenance_dict = {
                        "caption": f"Multimodel standard deviation of {cfg['quickplot']['title_var']} at {gwl} $^\\circ$ C",
                        "ancestors": ancestors_file,
                        "authors": ["swaminathan_ranjini"],
                        "references": ["swaminathan22jclim"],
                        "projects": ["ukesm"],
                        "domains": ["global"],
                        "plot_type": ["map"],
                    }
                    log_provenance(provenance_dict, stdev_file, cfg)
                    provenance_dict.update({"plot_type": ["map"]})

                    filename = (
                        SEP.join([project, "mm_stdev", str(gwl)]) + ".png"
                    )
                    filename = os.path.join(cfg["plot_dir"], filename)
                    log_provenance(provenance_dict, filename, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
