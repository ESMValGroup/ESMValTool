"""Utility functions for calculations."""

import json
import logging
import warnings
from pathlib import Path
from typing import Final

import dask.array as da
import iris
import iris.analysis.cartography
import numpy as np
import pandas as pd
from plot_utils import plot_corr_rmse_heatmaps, plot_maps
from wt_utils import (
    check_mapping_dict_format,
    get_driver,
    get_mapping_dict,
    get_provenance_record,
    log_provenance,
    map_lwt_to_slwt,
    reverse_convert_dict,
    write_corr_rmse_to_csv,
    write_mapping_dict,
)

iris.FUTURE.datum_support = True

logger = logging.getLogger(Path(__file__).name)

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")


# defining constants for direction calculations
DIR_0: Final = 0.0
DIR_22_5: Final = 22.5
DIR_67_5: Final = 67.5
DIR_112_5: Final = 112.5
DIR_157_5: Final = 157.5
DIR_202_5: Final = 202.5
DIR_247_5: Final = 247.5
DIR_292_5: Final = 292.5
DIR_337_5: Final = 337.5

# defining constant for flow limit
FLOW_LIMIT: Final = 6


def calc_slwt_obs(
    cfg: dict,
    lwt: np.array,
    cube: iris.cube.Cube,
    dataset: str,
    ancestors: list,
    timerange: str,
) -> np.array:
    """Calculate simplified weathertypes for observational data.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    lwt
        Array of Lamb WT
    cube
        Cube of psl data
    dataset
        Name of dataset
    ancestors
        List of ancestors
    timerange
        Time range for the calculation

    Returns
    -------
        np.array: Simplified Lamb Weathertypes

    """
    logger.info("Calculating simplified Lamb Weathertypes for %s", dataset)

    wt_data_prcp = []
    for wt_ in range(1, 28):
        target_indices = da.where(lwt == wt_)[0].compute()
        if target_indices.size == 0:
            logger.info(
                "calc_slwt_obs - CAUTION: Skipped wt %s \
                for dataset %s!",
                wt_,
                dataset,
            )
            continue
        dates = [
            cube.coord("time").units.num2date(cube.coord("time").points[i])
            for i in target_indices
        ]
        if dataset == "E-OBS":
            extracted_cube = cube[target_indices]
        else:
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d),
            )
        wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)

        data = wt_cube_mean.data

        if isinstance(data, np.ma.MaskedArray):
            wt_data_prcp.append(data.compressed())
        else:
            # if not MaskedArray, flatten should be fine
            wt_data_prcp.append(data.flatten())

    selected_pairs = process_prcp_mean(
        cfg,
        wt_data_prcp,
        dataset,
        timerange,
    )

    with Path(
        f"{cfg.get('work_dir')}/wt_selected_pairs_{dataset}.json",
    ).open("w", encoding="utf-8") as file:
        json.dump(selected_pairs, file)

    mapping_dict = get_mapping_dict(selected_pairs)

    write_mapping_dict(cfg.get("work_dir"), dataset, mapping_dict)

    provenance_record = get_provenance_record(
        "Lamb Weathertypes",
        ancestors,
        ["Lamb Weathertypes"],
    )

    log_provenance(
        f"{cfg.get('work_dir')}/wt_selected_pairs_{dataset}",
        cfg,
        provenance_record,
    )

    return map_lwt_to_slwt(lwt, mapping_dict)


def calc_const():
    """Calculate constants for weathertyping algorithm.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Returns
    -------
    tuple
        The four constants needed for WT calculation.
    """
    const1 = 1 / np.cos(45 * np.pi / 180)
    const2 = np.sin(45 * np.pi / 180) / np.sin(40 * np.pi / 180)
    const3 = np.sin(45 * np.pi / 180) / np.sin(50 * np.pi / 180)
    const4 = 1 / (2 * np.cos(45 * np.pi / 180) ** 2)

    return const1, const2, const3, const4


def calc_westerly_flow(cube: iris.cube.Cube) -> np.array:
    """Calculate the westerly flow over area.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    cube
        Cube of psl data.

    Returns
    -------
    np.array
        westerly flow
    """
    return 1 / 2 * (cube.data[:, 1, 2] + cube.data[:, 1, 4]) - 1 / 2 * (
        cube.data[:, 3, 2] + cube.data[:, 3, 4]
    )


def calc_southerly_flow(cube: iris.cube.Cube, const1: float) -> np.array:
    """Calculate the southerly flow over area.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    cube
        cube of psl data
    const1
        const1

    Returns
    -------
    np.array
        southerly flow
    """
    return const1 * (
        1
        / 4
        * (cube.data[:, 3, 4] + 2 * cube.data[:, 2, 4] + cube.data[:, 1, 4])
        - 1
        / 4
        * (cube.data[:, 3, 2] + 2 * cube.data[:, 2, 2] + cube.data[:, 1, 2])
    )


def calc_resultant_flow(
    westerly_flow: np.array,
    southerly_flow: np.array,
) -> np.array:
    """Calculate the resultant flow.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    westerly_flow
        westerly flow
    southerly_flow
        southerly flow

    Returns
    -------
    np.array
        resultant flow
    """
    return (southerly_flow**2 + westerly_flow**2) ** (1 / 2)


def calc_westerly_shear_velocity(
    cube: iris.cube.Cube,
    const2: float,
    const3: float,
) -> np.array:
    """Calculate westerly shear velocity.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    cube
        cube of psl data
    const2
        const2
    const3
        const3

    Returns
    -------
    np.array
        westerly shear velocity
    """
    return const2 * (
        1 / 2 * (cube.data[:, 0, 2] + cube.data[:, 0, 4])
        - 1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
    ) - const3 * (
        1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
        - 1 / 2 * (cube.data[:, 4, 2] + cube.data[:, 4, 4])
    )


def calc_southerly_shear_velocity(
    cube: iris.cube.Cube,
    const4: float,
) -> np.array:
    """Calculate southerly shear velocity.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    cube
        cube of psl data
    const4
        const

    Returns
    -------
    np.array
        southerly shear velocity
    """
    return const4 * (
        1
        / 4
        * (cube.data[:, 3, 6] + 2 * cube.data[:, 2, 6] + cube.data[:, 1, 6])
        - 1
        / 4
        * (cube.data[:, 3, 4] + 2 * cube.data[:, 2, 4] + cube.data[:, 1, 4])
        - 1
        / 4
        * (cube.data[:, 3, 2] + 2 * cube.data[:, 2, 2] + cube.data[:, 1, 2])
        + 1
        / 4
        * (cube.data[:, 3, 0] + 2 * cube.data[:, 2, 0] + cube.data[:, 1, 0])
    )


def calc_total_shear_velocity(
    westerly_shear_velocity: np.array,
    southerly_shear_velocity: np.array,
) -> np.array:
    """Calculate total shear velocity.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    westerly_shear_velocity
        westerly shear velocity
    southerly_shear_velocity
        southerly shear velocity

    Returns
    -------
    np.array
        total shear velocity
    """
    return westerly_shear_velocity + southerly_shear_velocity


def lamp_pure_directional_type(
    direction: float,
    weathertypes: np.array,
    i: int,
) -> int:
    """Calculate Lamp pure directional weathertype.

    Parameters
    ----------
    direction
        direction in degrees

    Returns
    -------
    int
        Lamb weathertype
    """
    if direction >= DIR_337_5 or direction < DIR_22_5:
        weathertypes[i] = 1
    if DIR_22_5 <= direction < DIR_67_5:
        weathertypes[i] = 2
    if DIR_67_5 <= direction < DIR_112_5:
        weathertypes[i] = 3
    if DIR_112_5 <= direction < DIR_157_5:
        weathertypes[i] = 4
    if DIR_157_5 <= direction < DIR_202_5:
        weathertypes[i] = 5
    if DIR_202_5 <= direction < DIR_247_5:
        weathertypes[i] = 6
    if DIR_247_5 <= direction < DIR_292_5:
        weathertypes[i] = 7
    if DIR_292_5 <= direction < DIR_337_5:
        weathertypes[i] = 8


def lamp_synoptic_directional_type(
    direction: float,
    weathertypes: np.array,
    i: int,
) -> int:
    """Calculate Lamb synoptic/directional hybrid weathertype.

    Parameters
    ----------
    direction
        direction in degrees

    Returns
    -------
    int
        Lamb weathertype
    """
    if direction >= DIR_337_5 or direction < DIR_22_5:
        weathertypes[i] = 11
    if DIR_22_5 <= direction < DIR_67_5:
        weathertypes[i] = 12
    if DIR_67_5 <= direction < DIR_112_5:
        weathertypes[i] = 13
    if DIR_112_5 <= direction < DIR_157_5:
        weathertypes[i] = 14
    if DIR_157_5 <= direction < DIR_202_5:
        weathertypes[i] = 15
    if DIR_202_5 <= direction < DIR_247_5:
        weathertypes[i] = 16
    if DIR_247_5 <= direction < DIR_292_5:
        weathertypes[i] = 17
    if DIR_292_5 <= direction < DIR_337_5:
        weathertypes[i] = 18


def lamb_synoptic_directional_type_zlt0(
    direction: float,
    weathertypes: np.array,
    i: int,
) -> int:
    """Calculate Lamb synoptic/directional hybrid weathertype with z <0.

    Parameters
    ----------
    direction
        direction in degrees

    Returns
    -------
    int
        Lamb weathertype
    """
    if direction >= DIR_337_5 or direction < DIR_22_5:
        weathertypes[i] = 19
    if DIR_22_5 <= direction < DIR_67_5:
        weathertypes[i] = 20
    if DIR_67_5 <= direction < DIR_112_5:
        weathertypes[i] = 21
    if DIR_112_5 <= direction < DIR_157_5:
        weathertypes[i] = 22
    if DIR_157_5 <= direction < DIR_202_5:
        weathertypes[i] = 23
    if DIR_202_5 <= direction < DIR_247_5:
        weathertypes[i] = 24
    if DIR_247_5 <= direction < DIR_292_5:
        weathertypes[i] = 25
    if DIR_292_5 <= direction < DIR_337_5:
        weathertypes[i] = 26


def wt_algorithm(cube: iris.cube.Cube, dataset: str) -> np.array:
    """Algorithm to calculate Lamb weathertypes.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Parameters
    ----------
    cube
        Cube of psl data
    dataset
        Name of dataset

    Returns
    -------
    np.array
        Lamb weathertypes
    """
    # lats and lons corresponding to datapoints
    # 55, 5 -> 1
    # 55, 15 -> 2
    # 50, -5 -> 3
    # 50, 5 -> 4
    # 50, 15 -> 5
    # 50, 25 -> 6
    # 45, -5 -> 7
    # 45, 5 -> 8
    # 45, 15 -> 9
    # 45, 25 -> 10
    # 40, -5 -> 11
    # 40, 5 -> 12
    # 40, 15 -> 13
    # 40, 25 -> 14
    # 35, 5 -> 15
    # 35, 15 -> 16

    # lons: -5, 0, 5, 10, 15, 20, 25
    # lats: 35, 40, 45, 50, 55
    logger.info("Calculating Lamb Weathertypes for %s", dataset)

    const1, const2, const3, const4 = calc_const()

    # westerly flow
    westerly_flow = calc_westerly_flow(cube)
    # southerly flow
    southerly_flow = calc_southerly_flow(cube, const1)
    # resultant flow
    total_flow = calc_resultant_flow(westerly_flow, southerly_flow)
    # westerly shear vorticity
    westerly_shear_velocity = calc_westerly_shear_velocity(
        cube,
        const2,
        const3,
    )
    # southerly shear vorticity
    southerly_shear_velocity = calc_southerly_shear_velocity(cube, const4)
    # total shear vorticity
    total_shear_velocity = calc_total_shear_velocity(
        westerly_shear_velocity,
        southerly_shear_velocity,
    )

    weathertypes = np.zeros(len(total_shear_velocity))

    for i, z_i in enumerate(total_shear_velocity):
        direction = (
            np.arctan(westerly_flow[i] / southerly_flow[i]) * 180 / np.pi
        )  # deg
        if southerly_flow[i] >= 0:
            direction += 180  # deg

        if direction < DIR_0:
            direction += 360  # deg

        # Lamb pure directional type
        if abs(z_i) < total_flow[i]:
            # writes directly to weathertypes array
            lamp_pure_directional_type(direction, weathertypes, i)
        # Lamb pure cyclonic and anticyclonic type
        elif (2 * total_flow[i]) < abs(z_i):
            if z_i > 0:
                weathertypes[i] = 9

            elif z_i < 0:
                weathertypes[i] = 10
        # Lamb synoptic/direction hybrid types
        elif total_flow[i] < abs(z_i) < (2 * total_flow[i]):
            if z_i > 0:
                # writes directly to weathertypes array
                lamp_synoptic_directional_type(direction, weathertypes, i)

            elif z_i < 0:
                # writes directly to weathertypes array
                lamb_synoptic_directional_type_zlt0(direction, weathertypes, i)
        # light indeterminate flow, corresponding to Lamb unclassified type U
        elif abs(z_i) < FLOW_LIMIT and total_flow[i] < FLOW_LIMIT:
            weathertypes[i] = 27

    return weathertypes


def calc_lwt_slwt_model(
    cfg: dict,
    cube: iris.cube.Cube,
    data_info: dict,
    predefined_slwt: dict | None = None,
):
    """Calculate Lamb as well as simplified weathertypes for model and write them to file.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    cube
        PSL field of dataset
    data_info
        Dictionary with info to dataset
    predefined_slwt
        If None, automatic_slwt will be used.
        If dict, this mapping dict will be used.
            (see recipe option predefined_slwt)
    """
    driver = get_driver(data_info)

    if not Path(
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}",
    ).exists():
        Path(
            f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}",
        ).mkdir(
            parents=True,
            exist_ok=True,
        )

    lwt = wt_algorithm(cube, data_info.get("dataset"))

    time_points = cube.coord("time").units.num2date(cube.coord("time").points)

    logger.info(
        "Calculating simplified Lamb Weathertypes for %s %s",
        data_info.get("dataset"),
        data_info.get("ensemble", ""),
    )

    if not predefined_slwt:
        with Path(
            f"{cfg.get('work_dir')}/wt_mapping_dict_ERA5.json",
        ).open("r", encoding="utf-8") as file:
            mapping_dict_era5_f = json.load(file)
            mapping_dict_era5 = reverse_convert_dict(mapping_dict_era5_f)

        with Path(
            f"{cfg.get('work_dir')}/wt_mapping_dict_E-OBS.json",
        ).open("r", encoding="utf-8") as file:
            mapping_dict_eobs_f = json.load(file)
            mapping_dict_eobs = reverse_convert_dict(mapping_dict_eobs_f)

        slwt_era5 = map_lwt_to_slwt(lwt, mapping_dict_era5)
        slwt_eobs = map_lwt_to_slwt(lwt, mapping_dict_eobs)

        print("ERA5: ", slwt_era5)
        print("EOBS: ", slwt_eobs)
    else:
        predefined_slwt = check_mapping_dict_format(predefined_slwt)
        write_mapping_dict(cfg.get("work_dir"), "ERA5", predefined_slwt)
        write_mapping_dict(cfg.get("work_dir"), "E-OBS", predefined_slwt)
        slwt_era5 = map_lwt_to_slwt(lwt, predefined_slwt)
        slwt_eobs = map_lwt_to_slwt(lwt, predefined_slwt)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(lwt, long_name="lwt"))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name="slwt_era5"))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name="slwt_eobs"))

    wt_cube[0].add_dim_coord(cube.coord("time"), 0)
    wt_cube[1].add_dim_coord(cube.coord("time"), 0)
    wt_cube[2].add_dim_coord(cube.coord("time"), 0)

    iris.save(
        wt_cube,
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}/{data_info.get('dataset')}{driver}_"
        f"{data_info.get('ensemble', '')}_{data_info.get('timerange')}.nc",
    )

    # write to csv file
    d = {
        "date": time_points[:],
        "lwt": np.int8(lwt),
        "slwt_ERA5": np.int8(slwt_era5),
        "slwt_EOBS": np.int8(slwt_eobs),
    }
    df = pd.DataFrame(data=d)
    df.to_csv(
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}/{data_info.get('dataset')}{driver}_{data_info.get('ensemble', '')}_"
        f"{data_info.get('timerange')}.csv",
        index=False,
    )

    ancestors = [
        data_info.get("preproc_path"),
        f"{cfg.get('work_dir')}/wt_mapping_dict_ERA5.json",
        f"{cfg.get('work_dir')}/wt_mapping_dict_E-OBS.json",
    ]
    provenance_record = get_provenance_record(
        "Lamb Weathertypes",
        ancestors,
        ["Lamb Weathertypes"],
    )

    log_provenance(
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}/{data_info.get('dataset')}_{data_info.get('ensemble', '')}",
        cfg,
        provenance_record,
    )


def rmse(subarray1: np.array, subarray2: np.array) -> np.array:
    """Calculate root mean square error between two arrays."""
    return np.sqrt(np.mean((subarray1 - subarray2) ** 2))


def process_prcp_mean(
    cfg: dict,
    data: np.array,
    dataset: str,
    timerange: str,
) -> list:
    """Return which weathertypes can be grouped together for a certain precipitation pattern.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    data
        Array of precipitation means for each WT
    dataset
        Name of dataset
    timerange
        Time range of dataset

    Returns
    -------
    list
        Selected pairs of WT. This is passed to get_mapping_dict
    """
    logger.info("Calculating corr and rsme matrices for %s", dataset)

    selected_pairs = []
    pattern_correlation_matrix = np.ma.corrcoef(data)
    n = len(data)
    rmse_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            rmse_matrix[i][j] = rmse(data[i], data[j])
            rmse_matrix[j][i] = rmse_matrix[i][j]
            if pattern_correlation_matrix[i][j] >= cfg.get(
                "correlation_threshold",
            ) and rmse_matrix[i][j] <= cfg.get("rmse_threshold"):
                selected_pairs.append(
                    (i + 1, j + 1),
                )

    # write matrices to csv
    write_corr_rmse_to_csv(
        cfg,
        pattern_correlation_matrix,
        rmse_matrix,
        dataset,
    )
    # plot heatmaps for matrices
    plot_corr_rmse_heatmaps(
        cfg,
        pattern_correlation_matrix,
        rmse_matrix,
        dataset,
        timerange,
    )

    return selected_pairs


def calc_wt_means(
    cfg: dict,
    cube: iris.cube.Cube,
    wt_cubes: iris.cube.CubeList,
    data_info: dict,
):
    """Calculate means for psl, tas or pr for weathertypes.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    cube
        Cube with variable data
    wt_cubes
        List of cubes of lwt, slwt_ERA5 and slwt_EOBS
    data_info
        Dictionary with info to dataset
    """
    var_name = data_info.get("var")
    wt_string = data_info.get("wt_string")
    driver = get_driver(data_info)

    logger.info(
        "Calculating %s %s means for %s",
        data_info.get("dataset"),
        var_name,
        wt_string,
    )

    wt_array, tcoord = get_wt_array(wt_string, wt_cubes)
    counts = da.bincount(wt_array.flatten().astype(int)).compute()

    for wt in range(1, len(counts)):
        target_indices = da.where(wt_array == wt)[0].compute()
        if target_indices.size == 0:
            logger.info(
                "calc_wt_means - CAUTION: Skipped %s %s \
                for dataset %s!",
                wt_string,
                wt,
                data_info.get("dataset"),
            )
            continue
        dates = [
            tcoord.units.num2date(tcoord.points[i]) for i in target_indices
        ]
        extracted_cube = cube.extract(
            iris.Constraint(time=lambda t, d=dates: t.point in d),
        )
        wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
        plot_maps(wt, cfg, wt_cube_mean, data_info, "mean")

        ancestors = [
            f"{data_info.get('preproc_path')}",
            f"{cfg.get('work_dir')}/ERA5.nc",
        ]
        provenance_record = get_provenance_record(
            f"{var_name} means for {wt_string}, wt: {wt}, ",
            ancestors,
            [var_name],
            ["map"],
            ["mean"],
        )

        local_path = f"{cfg.get('plot_dir')}/mean"

        log_provenance(
            f"{local_path}/{wt_string}_{wt}{driver}_"
            f"{data_info.get('dataset')}_{data_info.get('ensemble', '')}"
            f"_{var_name}_mean_{data_info.get('timerange')}",
            cfg,
            provenance_record,
        )


def get_wt_array(wt_string: str, wt_cubes: iris.cube.CubeList) -> tuple:
    """Get weathertype array and time coordinate based on wt_string.

    Parameters
    ----------
    wt_string
        string for weathertype selection
    wt_cubes
        list of weathertype cubes

    Raises
    ------
    NameError
        if wt_array does not exist for the given wt_string

    Returns
    -------
    tuple
        weathertype array and time coordinate
    """
    if wt_string == "slwt_ERA5":
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord("time")
        wt_array = slwt_era5_cube.core_data()[:]
    elif wt_string == "slwt_EOBS":
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord("time")
        wt_array = slwt_eobs_cube.core_data()[:]
    elif wt_string == "lwt":
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord("time")
        wt_array = lwt_cube.core_data()[:]
    else:
        error_str = "wt_array does not exist for calc_wt_means, calc_wt_anomalies or calc_wt_std."
        raise NameError(
            error_str,
        )

    return wt_array, tcoord


def calc_wt_anomalies(
    cfg: dict,
    cube: iris.cube.Cube,
    wt_cubes: iris.cube.CubeList,
    data_info: dict,
):
    """Calculate anomalies for psl, tas and pr for weathertypes.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    cube
        Cube with variable data
    wt_cubes
        List of cubes of lwt, slwt_ERA5 and slwt_EOBS
    data_info
        Dictionary with info to dataset
    """
    var_name = data_info.get("var_name")
    wt_string = data_info.get("wt_string")

    logger.info(
        "Calculating %s %s anomalies for %s",
        data_info.get("dataset"),
        var_name,
        wt_string,
    )

    wt_array, tcoord = get_wt_array(wt_string, wt_cubes)

    for wt in range(1, len(np.unique(wt_array)) + 1):
        target_indices = da.where(wt_array == wt).compute()
        if target_indices.size == 0:
            logger.info(
                "calc_wt_anomalies - CAUTION: Skipped wt %s \
                for dataset %s!",
                wt,
                data_info.get("dataset"),
            )
            continue
        dates = [
            tcoord.units.num2date(tcoord.points[i]) for i in target_indices
        ]
        extracted_cube = cube.extract(
            iris.Constraint(time=lambda t, d=dates: t.point in d[0]),
        )
        wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
        plot_maps(
            wt,
            cfg,
            cube.collapsed("time", iris.analysis.MEAN) - wt_cube_mean,
            data_info,
            "anomaly",
        )

        ancestors = [
            f"{data_info.get('preproc_path')}",
            f"{cfg.get('work_dir')}/ERA5.nc",
        ]
        provenance_record = get_provenance_record(
            f"{var_name} anomaly for, wt: {wt} \
                                                {wt_string}",
            ancestors,
            [var_name],
            ["map"],
            ["anomaly"],
        )

        log_provenance(
            f"{cfg.get('plot_dir')}/anomaly/{wt_string}_{wt}_{data_info.get('dataset')}_{data_info.get('ensemble', '')}"
            f"_{var_name}_anomaly__{data_info.get('timerange')}",
            cfg,
            provenance_record,
        )


def calc_wt_std(
    cfg: dict,
    cube: iris.cube.Cube,
    wt_cubes: iris.cube.CubeList,
    data_info: dict,
):
    """Calculate standard deviation for psl, tas and pr for weathertypes.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    cube
        Cube with variable data
    wt_cubes
        List of cubes of lwt, slwt_ERA5 and slwt_EOBS
    data_info
        Dictionary with info to dataset
    """
    var_name = data_info.get("var_name")
    wt_string = data_info.get("wt_string")

    logger.info(
        "Calculating %s %s standard deviation for %s",
        data_info.get("dataset"),
        var_name,
        wt_string,
    )

    wt_array, tcoord = get_wt_array(wt_string, wt_cubes)

    for wt in range(1, len(np.unique(wt_array)) + 1):
        target_indices = da.where(wt_array == wt).compute()
        if target_indices.size == 0:
            logger.info(
                "calc_slwt_obs - CAUTION: Skipped wt %s \
                for dataset %s!",
                wt,
                data_info.get("dataset"),
            )
            continue
        dates = [
            tcoord.units.num2date(tcoord.points[i]) for i in target_indices
        ]
        extracted_cube = cube.extract(
            iris.Constraint(time=lambda t, d=dates: t.point in d[0]),
        )
        wt_cube_std = extracted_cube.collapsed("time", iris.analysis.STD_DEV)
        plot_maps(wt, cfg, wt_cube_std, data_info, "stddev")

        ancestors = [
            f"{data_info.get('preproc_path')}",
            f"{cfg.get('work_dir')}/ERA5.nc",
        ]
        provenance_record = get_provenance_record(
            f"{var_name} standard, wt: {wt}\
                                                deviation for \
                                                {wt_string}",
            ancestors,
            [var_name],
            ["map"],
            ["stddev"],
        )

        local_path = f"{cfg.get('plot_dir')}/stddev"

        log_provenance(
            f"{local_path}/{wt_string}_{wt}_{data_info.get('dataset')}_{data_info.get('ensemble', '')}"
            f"_{var_name}_stddev_{data_info.get('timerange')}",
            cfg,
            provenance_record,
        )


def calc_lwt_model(cfg: dict, cube: iris.cube.Cube, data_info: dict):
    """Calculate Lamb weathertypes for models.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    cube
        Cube with variable data
    data_info
        Dictionary with info to dataset
    """
    driver = get_driver(data_info)

    if not Path(
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}",
    ).exists():
        Path(
            f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}",
        ).mkdir(parents=True, exist_ok=True)

    wt = wt_algorithm(cube, data_info.get("dataset"))

    time_points = cube.coord("time").units.num2date(cube.coord("time").points)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(wt, long_name="lwt"))

    logger.info(
        "Writing Lamb Weathertype for %s \
                to file %s.nc",
        data_info.get("dataset"),
        data_info.get("dataset"),
    )

    wt_cube[0].add_dim_coord(cube.coord("time"), 0)

    iris.save(
        wt_cube,
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}/{data_info.get('dataset')}{driver}"
        f"_{data_info.get('ensemble', '')}_{data_info.get('timerange')}.nc",
    )

    # write to csv file
    d = {"date": time_points[:], "lwt": np.int8(wt)}
    df = pd.DataFrame(data=d)
    df.to_csv(
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}/{data_info.get('dataset')}{driver}_{data_info.get('ensemble', '')}_"
        f"{data_info.get('timerange')}.csv",
        index=False,
    )

    ancestors = [
        f"{cfg.get('work_dir')}/wt_mapping_dict_ERA5.json",
        f"{cfg.get('work_dir')}/wt_mapping_dict_E-OBS.json",
    ]
    provenance_record = get_provenance_record(
        "Lamb Weathertypes",
        ancestors,
        ["Lamb Weathertypes"],
    )

    log_provenance(
        f"{cfg.get('work_dir')}/{data_info.get('output_file_path')}/{data_info.get('dataset')}_{data_info.get('ensemble', '')}",
        cfg,
        provenance_record,
    )


def plot_means(
    cfg: dict,
    preproc_var: np.array,
    wt_cubes: iris.cube.Cube,
    data_info: dict,
    mode: str = "slwt",
):
    """Plot means, anomalies and standard deviations.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    preproc_var
        Preprocessed variable cube
    wt_cubes
        List of cubes of lwt, slwt_ERA5 and slwt_EOBS
    data_info
        Dictionary with info to dataset
    mode
        Mode of weathertype calculation, either "slwt" or "lwt". Defaults to "slwt".
    """
    if mode == "slwt":
        data_info["wt_string"] = "lwt"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
        data_info["wt_string"] = "slwt_ERA5"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
        data_info["wt_string"] = "slwt_EOBS"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
    elif mode == "lwt":
        data_info["wt_string"] = "lwt"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
    else:
        e = "mode must be either 'slwt' or 'lwt'"
        raise ValueError(e)
