"""Utility functions for calculations."""

# operating system manipulations (e.g. path constructions)
import json
import logging
import os
import warnings

# to manipulate iris cubes
import iris
import iris.analysis.cartography

# general imports
import numpy as np
import pandas as pd

# local imports
from plot_utils import plot_corr_rmse_heatmaps, plot_maps
from wt_utils import (
    check_mapping_dict_format,
    get_mapping_dict,
    get_provenance_record,
    log_provenance,
    map_lwt_to_slwt,
    reverse_convert_dict,
    write_corr_rmse_to_csv,
    write_mapping_dict,
)

iris.FUTURE.datum_support = True

logger = logging.getLogger(os.path.basename(__file__))

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")


def calc_slwt_obs(
    cfg: dict,
    lwt: np.array,
    cube: iris.cube.Cube,
    dataset: str,
    ancestors: list,
    timerange: str,
) -> np.array:
    """Calculate simplified weathertypes for observation datasets based on \
        precipitation patterns over specified area.

    Args:
    ----
        cfg (dict): Configuration dictionary from recipe
        lwt (np.array): Array of Lamb WT
        cube (iris.cube.Cube): preprocessor cube to keep time coordinate
        dataset (str): Name of dataset
        correlation_thresold (float): correlation_threshold
        rmse_threshold (float): rsme_threshold
        ancestors (list): list of ancestors

    Returns
    -------
        np.array: _description_
    """
    logger.info("Calculating simplified Lamb Weathertypes for %s", dataset)

    work_dir = cfg.get("work_dir")
    correlation_threshold = cfg.get("correlation_threshold")
    rmse_threshold = cfg.get("rmse_threshold")
    tcoord = cube.coord("time")

    wt_data_prcp = []
    for wt_ in range(1, 28):
        target_indices = np.where(lwt == wt_)
        if len(target_indices[0]) < 1:
            logger.info(
                "calc_slwt_obs - CAUTION: Skipped wt %s \
                for dataset %s!",
                wt_,
                dataset,
            )
            continue
        dates = [
            tcoord.units.num2date(tcoord.points[i]) for i in target_indices
        ]
        if dataset == "E-OBS":
            extracted_cube = cube[target_indices]
        else:
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
        wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
        wt_data_prcp.append(wt_cube_mean.data.compressed())
    selected_pairs = process_prcp_mean(
        cfg,
        wt_data_prcp,
        correlation_threshold,
        rmse_threshold,
        dataset,
        timerange,
    )

    with open(
        f"{work_dir}/wt_selected_pairs_{dataset}.json", "w", encoding="utf-8"
    ) as file:
        json.dump(selected_pairs, file)

    mapping_dict = get_mapping_dict(selected_pairs)

    write_mapping_dict(work_dir, dataset, mapping_dict)

    provenance_record = get_provenance_record(
        "Lamb Weathertypes", ancestors, ["Lamb Weathertypes"], False, False
    )

    log_provenance(
        f"{work_dir}/wt_selected_pairs_{dataset}", cfg, provenance_record
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
        tuple: The four constants needed for WT calculation.
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

    Args:
    ----
        cube (iris.cube.Cube): Cube of psl data.

    Returns
    -------
        np.array: westerly flow
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

    Args:
    ----
        cube (iris.cube.Cube): Cube of psl data.
        const1 (float): const1

    Returns
    -------
        np.array: southerly flow
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
    westerly_flow: np.array, southerly_flow: np.array
) -> np.array:
    """Calculate the resultant flow.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
    ----
        westerly_flow (np.array): westerly flow.
        southerly_flow (np.array): southerly flow

    Returns
    -------
        np.array: resultant flow
    """
    return (southerly_flow**2 + westerly_flow**2) ** (1 / 2)


def calc_westerly_shear_velocity(
    cube: iris.cube.Cube, const2: float, const3: float
) -> np.array:
    """Calculate westerly shear velocity.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
    ----
        cube (iris.cube.Cube): cube of psl data
        const2 (float): const2
        const3 (float): const3

    Returns
    -------
        np.array: westerly shear velocity
    """
    return const2 * (
        1 / 2 * (cube.data[:, 0, 2] + cube.data[:, 0, 4])
        - 1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
    ) - const3 * (
        1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
        - 1 / 2 * (cube.data[:, 4, 2] + cube.data[:, 4, 4])
    )


def calc_southerly_shear_velocity(
    cube: iris.cube.Cube, const4: float
) -> np.array:
    """Calculate southerly shear velocity.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
    ----
        cube (iris.cube.Cube): cube of psl data
        const4 (float): const4

    Returns
    -------
        np.array: southerly shear velocity
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
    westerly_shear_velocity: np.array, southerly_shear_velocity: np.array
) -> np.array:
    """Calculate total shear velocity.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
    ----
        westerly_shear_velocity (np.array): westerly shear velocity
        southerly_shear_velocity (np.array): southerly shear velocity

    Returns
    -------
        np.array: total shear velocity
    """
    return westerly_shear_velocity + southerly_shear_velocity


def wt_algorithm(cube: iris.cube.Cube, dataset: str) -> np.array:
    """Algorithm to calculate Lamb weathertypes.

    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): PSL field of dataset
        dataset (str): Name of dataset

    Returns
    -------
        np.array: Array of Lamb WT for each day
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
        cube, const2, const3
    )
    # southerly shear vorticity
    southerly_shear_velocity = calc_southerly_shear_velocity(cube, const4)
    # total shear vorticity
    total_shear_velocity = calc_total_shear_velocity(
        westerly_shear_velocity, southerly_shear_velocity
    )

    weathertypes = np.zeros(len(total_shear_velocity))

    for i, z_i in enumerate(total_shear_velocity):
        direction = (
            np.arctan(westerly_flow[i] / southerly_flow[i]) * 180 / np.pi
        )  # deg
        if southerly_flow[i] >= 0:
            direction += 180  # deg

        if direction < 0:
            direction += 360  # deg

        # Lamb pure directional type
        if abs(z_i) < total_flow[i]:
            if 337.5 <= direction or direction < 22.5:
                weathertypes[i] = 1
            elif 22.5 <= direction < 67.5:
                weathertypes[i] = 2
            elif 67.5 <= direction < 112.5:
                weathertypes[i] = 3
            elif 112.5 <= direction < 157.5:
                weathertypes[i] = 4
            elif 157.5 <= direction < 202.5:
                weathertypes[i] = 5
            elif 202.5 <= direction < 247.5:
                weathertypes[i] = 6
            elif 247.5 <= direction < 292.5:
                weathertypes[i] = 7
            elif 292.5 <= direction < 337.5:
                weathertypes[i] = 8
        # Lamb’s pure cyclonic and anticyclonic type
        elif (2 * total_flow[i]) < abs(z_i):
            if z_i > 0:
                weathertypes[i] = 9

            elif z_i < 0:
                weathertypes[i] = 10
        # Lambs’s synoptic/direction hybrid types
        elif total_flow[i] < abs(z_i) < (2 * total_flow[i]):
            if z_i > 0:
                if 337.5 <= direction or direction < 22.5:
                    weathertypes[i] = 11
                elif 22.5 <= direction < 67.5:
                    weathertypes[i] = 12
                elif 67.5 <= direction < 112.5:
                    weathertypes[i] = 13
                elif 112.5 <= direction < 157.5:
                    weathertypes[i] = 14
                elif 157.5 <= direction < 202.5:
                    weathertypes[i] = 15
                elif 202.5 <= direction < 247.5:
                    weathertypes[i] = 16
                elif 247.5 <= direction < 292.5:
                    weathertypes[i] = 17
                elif 292.5 <= direction < 337.5:
                    weathertypes[i] = 18

            elif z_i < 0:
                if 337.5 <= direction or direction < 22.5:
                    weathertypes[i] = 19
                elif 22.5 <= direction < 67.5:
                    weathertypes[i] = 20
                elif 67.5 <= direction < 112.5:
                    weathertypes[i] = 21
                elif 112.5 <= direction < 157.5:
                    weathertypes[i] = 22
                elif 157.5 <= direction < 202.5:
                    weathertypes[i] = 23
                elif 202.5 <= direction < 247.5:
                    weathertypes[i] = 24
                elif 247.5 <= direction < 292.5:
                    weathertypes[i] = 25
                elif 292.5 <= direction < 337.5:
                    weathertypes[i] = 26
        # light indeterminate flow, corresponding to Lamb’s unclassified type U
        elif abs(z_i) < 6 and total_flow[i] < 6:
            weathertypes[i] = 27

    return weathertypes


def calc_lwt_slwt_model(
    cfg: dict,
    cube: iris.cube.Cube,
    data_info: dict,
    predefined_slwt: bool | dict,
):
    """Calculate Lamb WT and simplified WT for model data.

    Args:
        cfg (dict): Configuration dicitonary from recipe
        cube (iris.cube.Cube): PSL field of dataset
        data_info (dict): Dictionary with info to dataset
        predefined_slwt (bool | dict): If False, automatic_slwt will be used.
        If dict, this mapping dict will be used.
        (see recipe option predefined_slwt)
    """
    work_dir = cfg.get("work_dir")
    dataset = data_info.get("dataset")
    preproc_path = data_info.get("preproc_path")
    output_file_path = data_info.get("output_file_path")
    ensemble = data_info.get("ensemble", "")
    timerange = data_info.get("timerange")
    driver = data_info.get("driver", "")
    if driver != "":
        driver = f"_{driver}"

    if not os.path.exists(f"{work_dir}/{output_file_path}"):
        os.makedirs(f"{work_dir}/{output_file_path}")

    lwt = wt_algorithm(cube, dataset)

    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    logger.info(
        "Calculating simplified Lamb Weathertypes for %s %s", dataset, ensemble
    )

    if not predefined_slwt:
        with open(
            f"{work_dir}/wt_mapping_dict_ERA5.json", encoding="utf-8"
        ) as file:
            mapping_dict_era5_f = json.load(file)

        with open(
            f"{work_dir}/wt_mapping_dict_E-OBS.json", encoding="utf-8"
        ) as file:
            mapping_dict_eobs_f = json.load(file)
        mapping_dict_era5 = reverse_convert_dict(mapping_dict_era5_f)
        mapping_dict_eobs = reverse_convert_dict(mapping_dict_eobs_f)

        slwt_era5 = map_lwt_to_slwt(lwt, mapping_dict_era5)
        slwt_eobs = map_lwt_to_slwt(lwt, mapping_dict_eobs)
    else:
        predefined_slwt = check_mapping_dict_format(predefined_slwt)
        write_mapping_dict(work_dir, "ERA5", predefined_slwt)
        write_mapping_dict(work_dir, "E-OBS", predefined_slwt)
        slwt_era5 = map_lwt_to_slwt(lwt, predefined_slwt)
        slwt_eobs = map_lwt_to_slwt(lwt, predefined_slwt)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(lwt, long_name="lwt"))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name="slwt_era5"))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name="slwt_eobs"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(
        wt_cube,
        f"{work_dir}/{output_file_path}/{dataset}{driver}_"
        f"{ensemble}_{timerange}.nc",
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
        f"{work_dir}/{output_file_path}/{dataset}{driver}_{ensemble}_"
        f"{timerange}.csv",
        index=False,
    )

    ancestors = [
        preproc_path,
        f"{work_dir}/wt_mapping_dict_ERA5.json",
        f"{work_dir}/wt_mapping_dict_E-OBS.json",
    ]
    provenance_record = get_provenance_record(
        "Lamb Weathertypes", ancestors, ["Lamb Weathertypes"], False, False
    )

    log_provenance(
        f"{work_dir}/{output_file_path}/{dataset}_{ensemble}",
        cfg,
        provenance_record,
    )


def rmse(subarray1: np.array, subarray2: np.array) -> np.array:
    """Calculate rsme.

    Args:
        subarray1 (np.array): array1
        subarray2 (np.array): array2

    Returns
    -------
        np.array: rsme array
    """
    return np.sqrt(np.mean((subarray1 - subarray2) ** 2))


def process_prcp_mean(
    cfg: dict,
    data: np.array,
    correlation_threshold: float,
    rmse_threshold: float,
    dataset: str,
    timerange: str,
) -> list:
    """Process precipitation fields.

    Specified area to get a list of selected pairs of weathertypes
    with the highest correlation (higher than correlation_threshold)
    and smallest RSME (smaller than rsme_threshold) for
    further processing and simplifying the WT.

    Args:
        cfg (dict): Configuration dictionary from recipe
        data (np.array): Precipitation data
        correlation_threshold (float): Correlation threshold
        rmse_threshold (float): RMSE threshold
        dataset (str): Name of dataset

    Returns
    -------
        list: Selected pairs of WT. This is passed to get_mapping_dict
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
            if (
                pattern_correlation_matrix[i][j] >= correlation_threshold
                and rmse_matrix[i][j] <= rmse_threshold
            ):
                selected_pairs.append(
                    (
                        (i + 1, j + 1),
                        pattern_correlation_matrix[i][j],
                        rmse_matrix[i][j],
                    )
                )

    # write matrices to csv
    write_corr_rmse_to_csv(
        cfg, pattern_correlation_matrix, rmse_matrix, dataset
    )
    # plot heatmaps for matrices
    plot_corr_rmse_heatmaps(
        cfg, pattern_correlation_matrix, rmse_matrix, dataset, timerange
    )

    return selected_pairs


def calc_wt_means(
    cfg: dict,
    cube: iris.cube.Cube,
    wt_cubes: iris.cube.CubeList,
    data_info: dict,
):
    """Calculate means of fields of each weathertype.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube with variable data
        wt_cubes (iris.cube.CubeList): List of cubes of lwt, slwt_ERA5
                                        and slwt_EOBS
        data_info (dict): Dictionary with info to dataset
        preproc_path (str): Ancestor path
    """
    dataset = data_info.get("dataset")
    var_name = data_info.get("var")
    wt_string = data_info.get("wt_string")
    preproc_path = data_info.get("preproc_path")
    ensemble = data_info.get("ensemble")
    timerange = data_info.get("timerange")
    driver = data_info.get("driver", "")
    if driver != "":
        driver = f"_{driver}"

    logger.info("Calculating %s %s means for %s", dataset, var_name, wt_string)

    work_dir = cfg.get("work_dir")

    num_slwt = 0
    target_indices = []
    lwt, slwt_eobs, slwt_era5 = [], [], []

    if wt_string == "slwt_ERA5":
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord("time")
        slwt_era5 = slwt_era5_cube.data[:]
        num_slwt = len(np.unique(slwt_era5))
    elif wt_string == "slwt_EOBS":
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord("time")
        slwt_eobs = slwt_eobs_cube.data[:]
        num_slwt = len(np.unique(slwt_eobs))
    elif wt_string == "lwt":
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord("time")
        lwt = lwt_cube.data[:]

    if "slwt" in wt_string:
        for wt in range(1, num_slwt + 1):
            if wt_string == "slwt_ERA5":
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == "slwt_EOBS":
                target_indices = np.where(slwt_eobs == wt)
            else:
                logger.info("WT_STRING not supported!")
            if len(target_indices[0]) < 1:
                logger.info(
                    "calc_wt_means - CAUTION: Skipped %s %s \
                    for dataset %s!",
                    wt_string,
                    wt,
                    dataset,
                )
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
            wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
            plot_maps(wt, cfg, wt_cube_mean, data_info, "mean")
    elif wt_string == "lwt":
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    "calc_wt_means - CAUTION: Skipped lwt %s \
                    for dataset %s!",
                    wt,
                    dataset,
                )
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
            wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
            plot_maps(wt, cfg, wt_cube_mean, data_info, "mean")
    else:
        logger.info("WT_STRING NOT SUPPORTED.")

    ancestors = [f"{preproc_path}", f"{work_dir}/ERA5.nc"]
    provenance_record = get_provenance_record(
        f"{var_name} means for \
                                              {wt_string}",
        ancestors,
        [var_name],
        ["map"],
        ["mean"],
    )

    local_path = f"{cfg.get('plot_dir')}/mean"

    log_provenance(
        f"{local_path}/{wt_string}_{wt}{driver}_"
        f"{dataset}_{ensemble}"
        f"_{var_name}_mean_{timerange}",
        cfg,
        provenance_record,
    )


def calc_wt_anomalies(
    cfg: dict,
    cube: iris.cube.Cube,
    wt_cubes: iris.cube.CubeList,
    data_info: dict,
):
    """Calculate anomalies of fields of each weathertype.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube with variable data
        wt_cubes (iris.cube.CubeList): List of cubes of lwt, slwt_ERA5
                                        and slwt_EOBS
        data_info (dict): Dictionary with info to dataset
        preproc_path (str): Ancestor path
    """
    work_dir = cfg.get("work_dir")
    dataset = data_info.get("dataset")
    var_name = data_info.get("var_name")
    wt_string = data_info.get("wt_string")
    preproc_path = data_info.get("preproc_path")
    ensemble = data_info.get("ensemble")
    timerange = data_info.get("timerange")

    logger.info(
        "Calculating %s %s anomalies for %s", dataset, var_name, wt_string
    )

    target_indices = []
    lwt, slwt_eobs, slwt_era5 = [], [], []

    if wt_string == "slwt_ERA5":
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord("time")
        slwt_era5 = slwt_era5_cube.data[:]
    elif wt_string == "slwt_EOBS":
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord("time")
        slwt_eobs = slwt_eobs_cube.data[:]
    elif wt_string == "lwt":
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord("time")
        lwt = lwt_cube.data[:]

    num_slwt_era5 = len(np.unique(slwt_era5))
    num_slwt_eobs = len(np.unique(slwt_eobs))

    if num_slwt_eobs != num_slwt_era5:
        logger.info(
            "calc_wt_anomalies - CAUTION: unequal number of \
                    slwt_era5 (%s) and slwt_eobs (%s)!",
            num_slwt_era5,
            num_slwt_eobs,
        )

    if "slwt" in wt_string:
        for wt in range(1, max(num_slwt_era5, num_slwt_eobs)):
            if wt_string == "slwt_ERA5":
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == "slwt_EOBS":
                target_indices = np.where(slwt_eobs == wt)
            else:
                logger.info("WT_STRING not supported!")
            if len(target_indices[0]) < 1:
                logger.info(
                    "calc_wt_anomalies - CAUTION: Skipped wt %s \
                    for dataset %s!",
                    wt,
                    dataset,
                )
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
            wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
            plot_maps(
                wt,
                cfg,
                cube.collapsed("time", iris.analysis.MEAN) - wt_cube_mean,
                data_info,
                "anomaly",
            )
    elif wt_string == "lwt":
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    "calc_wt_anomalies - CAUTION: Skipped wt %s \
                    for dataset %s!",
                    wt,
                    dataset,
                )
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
            wt_cube_mean = extracted_cube.collapsed("time", iris.analysis.MEAN)
            plot_maps(
                wt,
                cfg,
                cube.collapsed("time", iris.analysis.MEAN) - wt_cube_mean,
                data_info,
                "anomaly",
            )
    else:
        logger.info("WT_STRING NOT SUPPORTED.")

    ancestors = [f"{preproc_path}", f"{work_dir}/ERA5.nc"]
    provenance_record = get_provenance_record(
        f"{var_name} anomaly for \
                                              {wt_string}",
        ancestors,
        [var_name],
        ["map"],
        ["anomaly"],
    )

    local_path = f"{cfg.get('plot_dir')}/anomaly"

    log_provenance(
        f"{local_path}/{wt_string}_{wt}_{dataset}_{ensemble}"
        f"_{var_name}_anomaly__{timerange}",
        cfg,
        provenance_record,
    )


def calc_wt_std(
    cfg: dict,
    cube: iris.cube.Cube,
    wt_cubes: iris.cube.CubeList,
    data_info: dict,
):
    """Calculate standard deviation of fields of each weathertype.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube with variable data
        wt_cubes (iris.cube.CubeList): List of cubes of lwt, slwt_ERA5
                                        and slwt_EOBS
        data_info (dict): Dictionary with info to dataset
        preproc_path (str): Ancestor path
    """
    work_dir = cfg.get("work_dir")
    dataset = data_info.get("dataset")
    var_name = data_info.get("var_name")
    wt_string = data_info.get("wt_string")
    preproc_path = data_info.get("preproc_path")
    ensemble = data_info.get("ensemble")
    timerange = data_info.get("timerange")

    logger.info(
        "Calculating %s %s standard deviation for %s",
        dataset,
        var_name,
        wt_string,
    )

    target_indices = []
    lwt, slwt_eobs, slwt_era5 = [], [], []

    if wt_string == "slwt_ERA5":
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord("time")
        slwt_era5 = slwt_era5_cube.data[:]
    elif wt_string == "slwt_EOBS":
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord("time")
        slwt_eobs = slwt_eobs_cube.data[:]
    elif wt_string == "lwt":
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord("time")
        lwt = lwt_cube.data[:]

    num_slwt_era5 = len(np.unique(slwt_era5))
    num_slwt_eobs = len(np.unique(slwt_eobs))

    if num_slwt_eobs != num_slwt_era5:
        logger.info(
            "calc_wt_std - CAUTION: unequal number of \
                    slwt_era5 (%s) and slwt_eobs (%s)!",
            num_slwt_era5,
            num_slwt_eobs,
        )

    if "slwt" in wt_string:
        for wt in range(1, max(num_slwt_era5, num_slwt_eobs)):
            if wt_string == "slwt_ERA5":
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == "slwt_EOBS":
                target_indices = np.where(slwt_eobs == wt)
            else:
                logger.info("WT_STRING not supported!")
            if len(target_indices[0]) < 1:
                logger.info(
                    "calc_slwt_obs - CAUTION: Skipped wt %s \
                    for dataset %s!",
                    wt,
                    dataset,
                )
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
            wt_cube_std = extracted_cube.collapsed(
                "time", iris.analysis.STD_DEV
            )
            plot_maps(wt, cfg, wt_cube_std, data_info, "stddev")
    elif wt_string == "lwt":
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    "calc_wt_std - CAUTION: Skipped wt %s \
                    for dataset %s!",
                    wt,
                    dataset,
                )
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t, d=dates: t.point in d[0])
            )
            wt_cube_std = extracted_cube.collapsed(
                "time", iris.analysis.STD_DEV
            )
            plot_maps(wt, cfg, wt_cube_std, data_info, "stddev")
    else:
        logger.info("WT_STRING NOT SUPPORTED.")

    ancestors = [f"{preproc_path}", f"{work_dir}/ERA5.nc"]
    provenance_record = get_provenance_record(
        f"{var_name} standard \
                                              deviation for \
                                              {wt_string}",
        ancestors,
        [var_name],
        ["map"],
        ["stddev"],
    )

    local_path = f"{cfg.get('plot_dir')}/stddev"

    log_provenance(
        f"{local_path}/{wt_string}_{wt}_{dataset}_{ensemble}"
        f"_{var_name}_stddev_{timerange}",
        cfg,
        provenance_record,
    )


def calc_lwt_model(
    cfg: dict, cube: iris.cube.Cube, dataset: str, data_info: dict
):
    """Calculate lwt for model data.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube to keep time coordinate
        dataset (str): Name of dataset
        data_info (dict): Dictionary with info to dataset
    """
    work_dir = cfg.get("work_dir")
    dataset = data_info.get("dataset")
    output_file_path = data_info.get("output_file_path")
    ensemble = data_info.get("ensemble", "")
    timerange = data_info.get("timerange")
    driver = data_info.get("driver", "")
    if driver != "":
        driver = f"_{driver}"

    if not os.path.exists(f"{work_dir}/{output_file_path}"):
        os.makedirs(f"{work_dir}/{output_file_path}")

    wt = wt_algorithm(cube, dataset)

    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(wt, long_name="lwt"))

    logger.info(
        "Writing Lamb Weathertype for %s \
                to file %s.nc",
        dataset,
        dataset,
    )

    wt_cube[0].add_dim_coord(tcoord, 0)

    iris.save(
        wt_cube,
        f"{work_dir}/{output_file_path}/{dataset}{driver}"
        f"_{ensemble}_{timerange}.nc",
    )

    # write to csv file
    d = {"date": time_points[:], "lwt": np.int8(wt)}
    df = pd.DataFrame(data=d)
    df.to_csv(
        f"{work_dir}/{output_file_path}/{dataset}{driver}_{ensemble}_"
        f"{timerange}.csv",
        index=False,
    )

    ancestors = [
        f"{work_dir}/wt_mapping_dict_ERA5.json",
        f"{work_dir}/wt_mapping_dict_E-OBS.json",
    ]
    provenance_record = get_provenance_record(
        "Lamb Weathertypes", ancestors, ["Lamb Weathertypes"], False, False
    )

    log_provenance(
        f"{work_dir}/{output_file_path}/{dataset}_{ensemble}",
        cfg,
        provenance_record,
    )


def plot_means(
    cfg: dict,
    preproc_var: np.array,
    wt_cubes: iris.cube.Cube,
    data_info: dict,
    only_lwt=False,
):
    """Wrapper function to plot various means/std/anomalies.

    Args:
        cfg (dict): cfg dictionary provided by recipe
        preproc_var (np.array): variable to be plotted
        wt_cubes (iris.cube.Cube): list of wt cubes
        data_info (dict): dictionary with info to dataset
        only_lwt (bool, optional): If True,
        only Lamb weathertypes will be loaded. Defaults to False.
        (useful for automatic_slwt = False)
    """
    if not only_lwt:
        data_info["wt_string"] = "lwt"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
        data_info["wt_string"] = "slwt_ERA5"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
        data_info["wt_string"] = "slwt_EOBS"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
    else:
        data_info["wt_string"] = "lwt"
        calc_wt_means(cfg, preproc_var, wt_cubes, data_info)
