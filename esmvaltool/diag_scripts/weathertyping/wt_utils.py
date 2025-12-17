"""Utility functions for weathertyping script."""

import json
import logging
import warnings
from pathlib import Path

import iris
import iris.analysis.cartography
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import ProvenanceLogger

iris.FUTURE.datum_support = True

logger = logging.getLogger(Path(__file__).name)

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")


def get_driver(data_info: dict) -> str:
    """Get driving model name and string for further use.

    Parameters
    ----------
    data_info
        Data information dictionary.

    Returns
    -------
    str:
        Driver string with leading underscore or empty string.
    """
    return data_info.get("driver", "")


def load_wt_preprocessors(dataset: str, preproc_variables_dict: dict) -> tuple:
    """Load preprocessor cubes for calculating Lamb weathertypes.

    Parameters
    ----------
    dataset
        Name of dataset
    preproc_variables_dict
        Dictionary with info on preprocessor variables

    Returns
    -------
    tuple
        Preprocessor cubes for weathertyping
    """
    wt_preproc = iris.load_cube(
        preproc_variables_dict.get(dataset)[0].get("filename"),
    )
    # try and except block to catch missing precipitation preprocessor
    # if it is not supplied, predefined_slwt has to be given by user
    try:
        wt_preproc_prcp = iris.load_cube(
            preproc_variables_dict.get(dataset)[1].get("filename"),
        )
        wt_preproc_prcp_eobs = iris.load_cube(
            preproc_variables_dict.get("E-OBS")[0].get("filename"),
        )
    except (IndexError, KeyError, TypeError):
        logger.info(
            "ERA5 precipitation preprocessor not found for automatic slwt.",
        )
        wt_preproc_prcp = None
        wt_preproc_prcp_eobs = None

    return wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs


def get_ancestors_era5_eobs(
    dataset: str,
    preproc_variables_dict: dict,
) -> tuple:
    """Get ancestors for observational data.

    Parameters
    ----------
    dataset
        Name of dataset
    preproc_variables_dict
        Dictionary with info on preprocessor variables

    Returns
    -------
    tuple
        Lists of ERA5 and E-OBS ancestors
    """
    era5_ancestors = [
        preproc_variables_dict.get(dataset)[0].get("filename"),
    ]
    try:
        era5_ancestors.append(
            preproc_variables_dict.get(dataset)[1].get("filename"),
        )
        eobs_ancestors = [
            preproc_variables_dict.get(dataset)[0].get("filename"),
            preproc_variables_dict.get("E-OBS")[0].get("filename"),
        ]
    except (IndexError, KeyError, TypeError):
        logger.info(
            "No ancestors for ERA5 and E-OBS precipitation preprocessor.",
        )
        eobs_ancestors = []

    return era5_ancestors, eobs_ancestors


def get_model_output_filepath(dataset: str, data_info: list) -> tuple:
    """Get output filepaths for models.

    Parameters
    ----------
    dataset
        Name of dataset
    data_info
        Model variables

    Returns
    -------
    tuple
        Output filepath and preprocessor path for future referencing.
    """
    timerange = data_info.get("timerange").replace("/", "-")
    experiment = data_info.get("exp")
    ensemble = data_info.get("ensemble")
    driver = data_info.get("driver", "")

    out_path = f"{dataset}/{driver}/{experiment}/{ensemble}/{timerange}"
    preproc_path = data_info.get("filename")

    return out_path, preproc_path


def get_preproc_lists(preproc_vars: list):
    """Put preprocessors and paths into lists for further use.

    Parameters
    ----------
    preproc_vars
        List of preprocessor variables.

    Returns
    -------
    tuple
        Preprocessor cubes and paths.
    """
    preproc_path_psl = preproc_vars[-3].get("filename")
    preproc_path_prcp = preproc_vars[-2].get("filename")
    preproc_path_tas = preproc_vars[-1].get("filename")

    mean_preproc_psl = iris.load_cube(preproc_vars[-3].get("filename"))
    mean_preproc_prcp = iris.load_cube(preproc_vars[-2].get("filename"))
    mean_preproc_tas = iris.load_cube(preproc_vars[-1].get("filename"))

    preproc_list = [mean_preproc_psl, mean_preproc_prcp, mean_preproc_tas]
    preproc_path_list = [preproc_path_psl, preproc_path_prcp, preproc_path_tas]

    return preproc_list, preproc_path_list


def get_preproc_lists_ensemble(preproc_vars: list):
    """Put preprocessors and paths into lists for further use.

    Parameters
    ----------
    preproc_vars
        List of preprocessor variables.

    Returns
    -------
    tuple
        Preprocessor cubes and paths.
    """
    preproc_path = preproc_vars.get("filename")

    preproc = iris.load_cube(preproc_vars.get("filename"))

    return preproc, preproc_path


def get_looping_dict(preproc_vars: list):
    """Put variable preprocessors into dict for looping.

    Parameters
    ----------
    preproc_vars
        List of preprocessor variables.

    Returns
    -------
    dict
        Dictionary of preprocessor cubes and paths.
    """
    preproc, preproc_path = get_preproc_lists(preproc_vars)

    return {
        "psl": [preproc[0], preproc_path[0]],
        "pr": [preproc[1], preproc_path[1]],
        "tas": [preproc[2], preproc_path[2]],
    }


def load_wt_files(path: str, mode="slwt"):
    """Load wt files.

    Parameters
    ----------
    path
        Path of wt file
    mode
        Type of weathertype to load. Defaults to "slwt".

    Returns
    -------
    list
        List of weathertype cubes.
    """
    if mode == "slwt":
        lwt_cube = iris.load_cube(path, "lwt")
        slwt_era5_cube = iris.load_cube(path, "slwt_era5")
        slwt_eobs_cube = iris.load_cube(path, "slwt_eobs")
        wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]
    elif mode == "lwt":
        lwt_cube = iris.load_cube(path, "lwt")
        wt_cubes = [lwt_cube]
    else:
        e = "Mode not recognized. Use 'slwt' or 'lwt'."
        raise ValueError(e)

    return wt_cubes


def get_provenance_record(
    caption: str,
    ancestors: list,
    long_names: list,
    plot_types: list | None = None,
    statistics: list | None = None,
) -> dict:
    """Get provenance record.

    Parameters
    ----------
    caption
        Caption of plot
    ancestors
        List of ancestor plots
    long_names
        List of variable long names
    plot_types
        Type of plot
    statistics
        Types of statistics used

    Returns
    -------
    dict
        Provenance record
    """
    record = {
        "caption": caption,
        "domains": ["reg"],
        "authors": ["jury_martin", "kroissenbrunner_thomas"],
        "references": ["maraun21jgr", "jones93ijc"],
        "projects": ["preval"],
        "long_names": long_names,
        "ancestors": ancestors,
    }
    if plot_types:
        record["plot_types"] = plot_types
    if statistics:
        record["statistics"] = statistics
    return record


def log_provenance(filename: str, cfg: dict, provenance_record: dict):
    """Log provenance.

    Parameters
    ----------
    filename:
        Filename of xml file
    cfg
        Configuration dictionary
    provenance_record
        Provenance record dictionary
    """
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info(
        "Provenance added to %s",
        f"{cfg['run_dir']}/diagnostic_provenance.yml",
    )


def turn_list_to_mapping_dict(list_: list) -> dict:
    """Turn list of combined WT to a dictionary for further use.

    Parameters
    ----------
    list_
        List of combined WTs

    Returns
    -------
    dict
        Dictionary of combined WTs
    """
    result_dict = {}

    for i, s in enumerate(list_):
        for elem in s:
            if elem not in result_dict:
                result_dict[elem] = i + 1

    return result_dict


def get_mapping_dict(selected_pairs: list) -> dict:
    """Get mapping dictionary from list of selected pairs.

    Parameters
    ----------
    selected_pairs
        List of selected weathertype pairs

    Returns
    -------
    dict
        Mapping dictionary
    """
    # selected pairs is of form [(wt1, wt2), (wt3, wt4), ...]
    sets_list = [set(i) for i in selected_pairs if i]

    # 3 loops to check all sets against each other until no more merges happen
    # merged_flag indicates if a merge happened in the last full pass -> we need to restart then
    merged_flag = True

    while merged_flag:
        merged_flag = False
        i = 0

        # check each set against each other
        while i < len(sets_list):
            set1 = sets_list[i]
            j = i + 1

            while j < len(sets_list):
                set2 = sets_list[j]

                if set1.isdisjoint(set2) is False:
                    # this means there is an overlap -> merge sets
                    set1.update(set2)
                    sets_list.pop(j)
                    merged_flag = True
                else:
                    j += 1

            i += 1

    return turn_list_to_mapping_dict(sets_list)


def write_mapping_dict(work_dir: str, dataset: str, mapping_dict: dict):
    """Write mapping dictionary to file.

    Parameters
    ----------
    work_dir
        Current working directory
    dataset
        Dataset name
    mapping_dict
        Mapping dictionary in {lwt: slwt, ...} format
    """
    mapping_dict_reformat = convert_dict(mapping_dict)

    with Path(
        f"{work_dir}/wt_mapping_dict_{dataset}.json",
    ).open("w", encoding="utf-8") as file:
        json.dump(mapping_dict_reformat, file)


def convert_dict(dict_: dict) -> dict:
    """Convert dictionary from {lwt: slwt, ...} format to {slwt: [lwt1, lwt2], ...}.

    Parameters
    ----------
    dict_
        Mapping dictionary to be converted

    Returns
    -------
    dict
        Converted dictionary
    """
    new_dict = {}
    for dataset, value in dict_.items():
        if value not in new_dict:
            new_dict[value] = []
        new_dict[value].append(dataset)
    return new_dict


def reverse_convert_dict(originial_dict: dict) -> dict:
    """Convert mapping dictionary.

    From {slwt: [lwt1, lwt2], ...} format to {lwt: slwt, ...}.

    Parameters
    ----------
    original_dict
        Dict in the {slwt: [lwt1, lwt2], ...} format

    Returns
    -------
    dict
        Dict in the  format {lwt: slwt, ...}
    """
    new_dict = {}
    for key, value_list in originial_dict.items():
        for original_key in value_list:
            new_dict[original_key] = key
    return new_dict


def write_corr_rmse_to_csv(
    cfg: dict,
    pattern_correlation_matrix: np.array,
    rmse_matrix: np.array,
    dataset: str,
):
    """Write correlation and rmse matrix to csv files.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    pattern_correlation_matrix
        Correlation matrix
    rmse_matrix
        RMSE matrix
    dataset
        Name of dataset
    """
    logger.info("Writing corr and rsme matrices for %s", dataset)

    work_dir = cfg.get("work_dir")

    df_corr = pd.DataFrame(pattern_correlation_matrix)
    df_corr.index = range(1, len(df_corr) + 1)
    df_corr.columns = range(1, len(df_corr.columns) + 1)
    df_corr.to_csv(
        f"{work_dir}/correlation_matrix_{dataset}.csv",
        index_label="Index",
    )

    df_rmse = pd.DataFrame(rmse_matrix)
    df_rmse.index = range(1, len(df_rmse) + 1)
    df_rmse.columns = range(1, len(df_rmse.columns) + 1)
    df_rmse.to_csv(
        f"{work_dir}/rmse_matrix_{dataset}.csv",
        index_label="Index",
    )


def run_predefined_slwt(
    work_dir: str,
    dataset_name: str,
    lwt: np.array,
    predefined_slwt: dict,
):
    """Run predefined slwt mapping.

    Parameters
    ----------
    work_dir
        Working directory to save mapping dict
    dataset_name
        Name of dataset
    lwt
        lwt array
    predefined_slwt
        Mapping dictionary in {lwt: slwt, ...} format

    Returns
    -------
    np.array
        slwt_era5 array
    np.array
        slwt_eobs array
    """
    predefined_slwt = check_mapping_dict_format(predefined_slwt)
    write_mapping_dict(work_dir, dataset_name, predefined_slwt)
    write_mapping_dict(work_dir, "E-OBS", predefined_slwt)
    slwt_era5 = map_lwt_to_slwt(lwt, predefined_slwt)
    slwt_eobs = map_lwt_to_slwt(lwt, predefined_slwt)

    return slwt_era5, slwt_eobs


def combine_wt_to_file(
    cfg: dict,
    wt_list: list,
    cube: iris.cube.Cube,
    file_name: str,
):
    """Combine lwt and slwt arrays to one file.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    wt_list
        List of weathertype arrays
    cube
        Cube of data to keep time coordinate
    file_name
        Name of output file
    """
    lwt = wt_list[0]
    slwt_era5 = wt_list[1]
    slwt_eobs = wt_list[2]

    logger.info("Writing weathertypes to %s", file_name)

    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg.get("work_dir")

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(lwt, long_name="lwt"))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name="slwt_era5"))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name="slwt_eobs"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f"{write_path}/{file_name}.nc")

    # write to csv file
    d = {
        "date": time_points[:],
        "lwt": np.int8(lwt),
        "slwt_era5": np.int8(slwt_era5),
        "slwt_eobs": np.int8(slwt_eobs),
    }
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f"/{file_name}.csv", index=False)


def write_lwt_to_file(
    cfg: dict,
    lwt: np.array,
    cube: iris.cube.Cube,
    file_name: str,
):
    """Write only lwt to file.

    Parameters
    ----------
    cfg
        Configuration dictionary from recipe
    lwt
        lwt array
    cube
        Cube of data to keep time coordinate
    file_name
        Name of output file
    """
    logger.info("Writing Lamb Weathertype to %s", file_name)

    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg.get("work_dir")

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(np.int8(lwt), long_name="lwt"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    iris.save(wt_cube, f"{write_path}/{file_name}.nc")

    # write to csv file
    d = {"date": time_points[:], "lwt": np.int8(lwt)}
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f"/{file_name}.csv", index=False)


def map_lwt_to_slwt(lwt: np.array, mapping_dict: dict) -> np.array:
    """Map lwt array to slwt array.

    Parameters
    ----------
    lwt
        lwt array
    mapping_dict
        Mapping dictionary in {lwt: slwt, ...} format

    Returns
    -------
    np.array
        array of slwt
    """
    return np.array([np.int8(mapping_dict.get(value, 0)) for value in lwt])


def check_mapping_dict_format(mapping_dict: dict) -> dict:
    """Check format of mapping dict and return in {lwt: slwt, ...} format.

    Parameters
    ----------
    mapping_dict
        mapping dict in any format

    Returns
    -------
    dict
        mapping dict in {lwt: slwt, ...} format
    """
    if isinstance(mapping_dict[next(iter(mapping_dict))], list):
        return reverse_convert_dict(mapping_dict)
    return mapping_dict
