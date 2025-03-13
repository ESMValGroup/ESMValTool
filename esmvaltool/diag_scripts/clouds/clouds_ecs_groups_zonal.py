"""Python diagnostic for plotting zonal averages."""

import logging
from copy import deepcopy
from pathlib import Path

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)

VAR_NAMES = {
    "clt": "total_cloud_fraction",
    "clivi": "ice_water_path",
    "lwp": "liquid_water_path",
    "swcre": "shortwave_cloud_radiative_effect",
    "lwcre": "longwave_cloud_radiative_effect",
    "netcre": "net_cloud_radiative_effect",
}
LINE_LEGEND = {
    "ECS_high_hist": "ECS_high",
    "ECS_med_hist": "ECS_med",
    "ECS_low_hist": "ECS_low",
}
LINE_COLOR = {
    "ECS_high_hist": "royalblue",
    "ECS_high_scen": "royalblue",
    "ECS_med_hist": "green",
    "ECS_med_scen": "green",
    "ECS_low_hist": "orange",
    "ECS_low_scen": "orange",
    "CMIP6": "firebrick",
    "CMIP5": "royalblue",
    "CMIP3": "darkcyan",
    "OBS": "black",
}
LINE_DASH = {
    "ECS_high_hist": "solid",
    "ECS_high_scen": "dashed",
    "ECS_med_hist": "solid",
    "ECS_med_scen": "dashed",
    "ECS_low_hist": "solid",
    "ECS_low_scen": "dashed",
    "CMIP6": "solid",
    "CMIP5": "solid",
    "CMIP3": "solid",
    "OBS": "solid",
}


def get_provenance_record(short_name, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = (
        f"Zonally averaged group means of {short_name} in the upper"
        "panel and the corresponding relative differences in lower"
        "panel."
    )

    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_types": ["zonal"],
        "authors": [
            "bock_lisa",
        ],
        "references": [
            "bock24acp",
        ],
        "ancestors": ancestor_files,
    }
    return record


def _get_multi_model_mean(cubes, var):
    """Compute multi-model mean."""
    logger.debug("Calculating multi-model mean")
    datasets = []
    mmm = []
    for dataset_name, cube in cubes.items():
        datasets.append(dataset_name)
        mmm.append(cube.data)
    mmm = np.ma.array(mmm)
    dataset_0 = list(cubes.keys())[0]
    mmm_cube = cubes[dataset_0].copy(data=np.ma.mean(mmm, axis=0))
    attributes = {
        "dataset": "MultiModelMean",
        "short_name": var,
        "datasets": "|".join(datasets),
    }
    mmm_cube.attributes = attributes
    return mmm_cube


def _get_multi_model_quantile(cubes, var, quantile):
    """Compute multi-model quantile."""
    logger.debug("Calculating multi-model %s quantile", quantile)
    datasets = []
    mmq = []
    for dataset_name, cube in cubes.items():
        datasets.append(dataset_name)
        mmq.append(cube.data)
    mmq = np.ma.array(mmq)
    dataset_0 = list(cubes.keys())[0]
    mmq_cube = cubes[dataset_0].copy(data=np.quantile(mmq, quantile, axis=0))
    attributes = {
        "dataset": "MultiModel" + str(quantile),
        "short_name": var,
        "datasets": "|".join(datasets),
    }
    mmq_cube.attributes = attributes
    return mmq_cube


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    if cube.var_name == "cli":
        cube.convert_units("g/kg")
    elif cube.var_name == "clw":
        cube.convert_units("g/kg")

    logger.debug("Reading %s", filename)
    cube = iris.util.squeeze(cube)
    return cube


def compute_diff(filename1, filename2):
    """Compute difference between two cubes."""
    logger.debug("Loading %s", filename1)
    cube1 = iris.load_cube(filename1)
    cube2 = iris.load_cube(filename2)

    if cube1.var_name == "cli":
        cube1.convert_units("g/kg")
        cube2.convert_units("g/kg")
    elif cube1.var_name == "clw":
        cube1.convert_units("g/kg")
        cube2.convert_units("g/kg")

    cube = cube2 - cube1
    cube.metadata = cube1.metadata
    cube = iris.util.squeeze(cube)
    return cube


def compute_diff_temp(input_data, group, dataset, plot_type):
    """Compute relative change per temperture change."""
    dataset_name = dataset["dataset"]
    var = dataset["short_name"]

    input_file_1 = dataset["filename"]

    var_data_2 = select_metadata(
        input_data,
        short_name=var,
        dataset=dataset_name,
        variable_group=group[1],
    )
    if not var_data_2:
        raise ValueError(
            f"No '{var}' data for '{dataset_name}' in '{group[1]}' available"
        )

    input_file_2 = var_data_2[0]["filename"]

    if plot_type == "zonal":
        ta_data_1 = select_metadata(
            input_data,
            short_name="tas",
            dataset=dataset_name,
            variable_group="tas_" + group[0],
        )
        ta_data_2 = select_metadata(
            input_data,
            short_name="tas",
            dataset=dataset_name,
            variable_group="tas_" + group[1],
        )
    elif plot_type == "height":
        ta_data_1 = select_metadata(
            input_data,
            short_name="ta",
            dataset=dataset_name,
            variable_group="ta_" + group[0],
        )
        ta_data_2 = select_metadata(
            input_data,
            short_name="ta",
            dataset=dataset_name,
            variable_group="ta_" + group[1],
        )
    else:
        raise ValueError(f"The plot_type '{var}' is not implemented.")

    if not ta_data_1:
        raise ValueError(
            f"No temperature data for '{dataset_name}' "
            f"in '{group[0]}' available"
        )
    if not ta_data_2:
        raise ValueError(
            f"No temperature data for '{dataset_name}' "
            f"in '{group[1]}' available"
        )
    input_file_ta_1 = ta_data_1[0]["filename"]
    input_file_ta_2 = ta_data_2[0]["filename"]

    cube = compute_diagnostic(input_file_1)
    if var in ["lwp", "clivi", "clw", "cli"]:
        cube.data[cube.data < 0.001] = np.nan
    elif var in ["cl"]:
        cube.data[cube.data < 0.1] = np.nan
    elif var in ["netcre", "swcre", "lwcre"]:
        cube.data[abs(cube.data) < 1.0] = np.nan

    cube_diff = compute_diff(input_file_1, input_file_2)
    cube_ta_diff = compute_diff(input_file_ta_1, input_file_ta_2)

    cube_ta_diff.data[cube_ta_diff.data < 1.0] = np.nan

    cube_diff = (
        100.0 * (cube_diff / iris.analysis.maths.abs(cube)) / cube_ta_diff
    )

    cube_diff.metadata = cube.metadata

    if plot_type == "zonal":
        logger.debug("Computing zonal mean")
        cube_diff = cube_diff.collapsed("longitude", iris.analysis.MEAN)
    elif plot_type == "height":
        logger.debug("Computing field mean")
        grid_areas = iris.analysis.cartography.area_weights(cube_diff)
        cube_diff = cube_diff.collapsed(
            ["longitude", "latitude"], iris.analysis.MEAN, weights=grid_areas
        )
    else:
        raise ValueError(f"Plot type {plot_type} is not implemented.")

    cube_diff.units = "%/K"

    return cube_diff


def plot_diagnostic(cube, legend, plot_type):
    """Create diagnostic data and plot it."""
    cube_label = legend
    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.subplot(211)

    if plot_type == "height":
        cube.coord("air_pressure").convert_units("hPa")
        y_axis = cube.coord("air_pressure")
        qplt.plot(
            cube,
            y_axis,
            label=cube_label,
            color=line_color,
            linestyle=line_dash,
        )
    else:
        lat = cube.coord("latitude")
        qplt.plot(
            lat, cube, label=cube_label, color=line_color, linestyle=line_dash
        )

    logger.info("Plotting %s", legend)


def plot_diagnostic_diff(cube, legend, plot_type):
    """Create diagnostic data and plot it."""
    cube_label = LINE_LEGEND.get(legend, legend)
    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.subplot(212)

    if cube.var_name == "pr":
        cube.units = cube.units / "kg m-3"
        cube.data = cube.core_data() / 1000.0
        cube.convert_units("mm day-1")
    elif cube.var_name == "cli":
        cube.convert_units("g/kg")
    elif cube.var_name == "clw":
        cube.convert_units("g/kg")

    if plot_type == "height":
        cube.coord("air_pressure").convert_units("hPa")
        y_axis = cube.coord("air_pressure")
        qplt.plot(
            cube,
            y_axis,
            label=cube_label,
            color=line_color,
            linestyle=line_dash,
        )
    else:
        lat = cube.coord("latitude")
        qplt.plot(
            lat, cube, label=cube_label, color=line_color, linestyle=line_dash
        )

    logger.info("Plotting %s", legend)


def plot_errorband(cube1, cube2, legend, plot_type):
    """Create diagnostic data and plot it."""
    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.subplot(211)

    if cube1.var_name == "pr":
        cube1.units = cube1.units / "kg m-3"
        cube1.data = cube1.core_data() / 1000.0
        cube1.convert_units("mm day-1")
        cube2.units = cube2.units / "kg m-3"
        cube2.data = cube2.core_data() / 1000.0
        cube2.convert_units("mm day-1")
    elif cube1.var_name == "cli":
        cube1.convert_units("g/kg")
        cube2.convert_units("g/kg")
    elif cube1.var_name == "clw":
        cube1.convert_units("g/kg")
        cube2.convert_units("g/kg")

    if plot_type == "height":
        cube1.coord("air_pressure").convert_units("hPa")
        cube2.coord("air_pressure").convert_units("hPa")
        y_axis = cube1.coord("air_pressure").points
        plt.fill_betweenx(
            y_axis,
            cube1.data,
            cube2.data,
            color=line_color,
            linestyle=line_dash,
            alpha=0.1,
        )
    else:
        lat = cube1.coord("latitude").points
        plt.fill_between(
            lat,
            cube1.data,
            cube2.data,
            color=line_color,
            linestyle=line_dash,
            alpha=0.1,
        )
    logger.info("Plotting %s", legend)


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault("filename_attach", "base")

    plot_type = cfg["plot_type"]

    input_data = list(cfg["input_data"].values())

    groups = group_metadata(input_data, "variable_group", sort="dataset")

    plt.figure(figsize=(8, 12))

    for group_name in groups:
        if ("tas_" not in group_name) and ("ta_" not in group_name):
            logger.info("Processing variable %s", group_name)

            dataset_names = []
            cubes = {}

            for dataset in groups[group_name]:
                dataset_name = dataset["dataset"]
                var = dataset["short_name"]

                if dataset_name not in [
                    "MultiModelMean",
                    "MultiModelP5",
                    "MultiModelP95",
                ]:
                    logger.info("Loop dataset %s", dataset_name)

                    input_file = dataset["filename"]
                    cube = compute_diagnostic(input_file)
                    logger.debug("Computing zonal mean")
                    if plot_type == "zonal":
                        cube = cube.collapsed("longitude", iris.analysis.MEAN)
                    elif plot_type == "height":
                        grid_areas = iris.analysis.cartography.area_weights(
                            cube
                        )
                        cube = cube.collapsed(
                            ["longitude", "latitude"],
                            iris.analysis.MEAN,
                            weights=grid_areas,
                        )
                    else:
                        raise ValueError(
                            f"Plot type {plot_type} is not implemented."
                        )

                    cubes[dataset_name] = cube

            cube_mmm = _get_multi_model_mean(cubes, var)

            plot_diagnostic(cube_mmm, group_name, plot_type)

            cube_p5 = _get_multi_model_quantile(cubes, var, 0.05)
            cube_p95 = _get_multi_model_quantile(cubes, var, 0.95)

            plot_errorband(cube_p5, cube_p95, group_name, plot_type)

    if plot_type == "height":
        plt.ylim(1000.0, 100.0)
        plt.yscale("log")
        plt.yticks(
            [1000.0, 800.0, 600.0, 400.0, 300.0, 200.0, 100.0],
            [1000, 800, 600, 400, 300, 200, 100],
        )

    long_name = input_data[0]["long_name"]
    if plot_type == "height":
        title = "Vertical mean of " + long_name
    elif plot_type == "zonal":
        if long_name == "Total Cloud Cover Percentage":
            title = "Zonal mean of Total Cloud Fraction"
        else:
            title = "Zonal mean of " + long_name
    else:
        title = long_name

    plt.title(title)
    plt.legend(ncol=1)
    plt.grid(True)

    for group_name in cfg["group_by"]:
        logger.info("Processing group %s", group_name[0])

        dataset_names = []
        cubes_diff = {}

        for dataset in groups[group_name[0]]:
            dataset_name = dataset["dataset"]
            var = dataset["short_name"]

            if dataset_name not in [
                "MultiModelMean",
                "MultiModelP5",
                "MultiModelP95",
            ]:
                logger.info("Loop dataset %s", dataset_name)
                dataset_names.append(dataset_name)

                cube_diff = compute_diff_temp(
                    input_data, group_name, dataset, plot_type
                )

                cubes_diff[dataset_name] = cube_diff

        cube_mmm = _get_multi_model_mean(cubes_diff, var)

        plot_diagnostic_diff(cube_mmm, group_name[0], plot_type)

    if plot_type == "height":
        plt.xlim(0.0, 1.0)
        plt.ylim(1000.0, 100.0)
        plt.yscale("log")
        plt.yticks(
            [1000.0, 800.0, 600.0, 400.0, 300.0, 200.0, 100.0],
            [1000, 800, 600, 400, 300, 200, 100],
        )
        plt.axvline(x=0, ymin=0.0, ymax=1.0, color="black", linewidth=3)
        title = "Difference of vertical mean of " + long_name
    elif plot_type == "zonal":
        plt.axhline(y=0, xmin=-90.0, xmax=90.0, color="black", linewidth=3)
        title = "Difference of zonal mean of " + long_name
    else:
        title = long_name

    plt.title(title)
    plt.legend(ncol=1)
    plt.grid(True)

    short_name = input_data[0]["short_name"]
    provenance_record = get_provenance_record(
        short_name, ancestor_files=[d["filename"] for d in input_data]
    )

    if plot_type == "height":
        basename = "level_diff_" + short_name + "_" + cfg["filename_attach"]
    else:
        basename = "zonal_diff_" + short_name + "_" + cfg["filename_attach"]

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cube_mmm)

    # And save the plot
    save_figure(basename, provenance_record, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
