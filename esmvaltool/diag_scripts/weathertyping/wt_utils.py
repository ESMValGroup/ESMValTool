"""Utility functions for weathertyping script."""

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

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import ProvenanceLogger, group_metadata

iris.FUTURE.datum_support = True

logger = logging.getLogger(os.path.basename(__file__))

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")


def get_cfg_vars(cfg: dict):
    """Get list of vars from configuration dict.

    Args:
    ----
        cfg : dict
            Configuration dict from recipe.

    Returns
    -------
        tuple
            cfg vars
    """
    preproc_variables_dict = group_metadata(
        cfg.get("input_data").values(), "dataset"
    )

    return (
        preproc_variables_dict,
        cfg.get("correlation_threshold"),
        cfg.get("rmse_threshold"),
        cfg.get("work_dir"),
        cfg.get("plotting", False),
        cfg.get("automatic_slwt", True),
        cfg.get("predefined_slwt", False),
    )


def load_wt_preprocessors(dataset: str, preproc_variables_dict: dict):
    """Load preprocessor cubes for calculating Lamb weathertypes.

    Args:
    ----
        dataset : str
            dataset name
        preproc_variables_dict : dict
            dictionary of preprocessor variables

    Returns
    -------
        list
            list of preprocessor vars for weathertyping
    """
    wt_preproc = iris.load_cube(
        preproc_variables_dict.get(dataset)[0].get("filename")
    )
    wt_preproc_prcp = iris.load_cube(
        preproc_variables_dict.get(dataset)[1].get("filename")
    )
    wt_preproc_prcp_eobs = iris.load_cube(
        preproc_variables_dict.get("E-OBS")[0].get("filename")
    )

    return wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs


def get_ancestors_era5_eobs(dataset: str, preproc_variables_dict: dict):
    """Get ancestors for ERA5/E-OBS.

    Args:
    ----
       dataset : str
            dataset name
        preproc_variables_dict : dict
            dictionary of preprocessor variables

    Returns
    -------
        tuple(list, list)
            lists of ERA5/E-OBS ancestors
    """
    era5_ancestors = [
        preproc_variables_dict.get(dataset)[0].get("filename"),
        preproc_variables_dict.get(dataset)[1].get("filename"),
    ]
    eobs_ancestors = [
        preproc_variables_dict.get(dataset)[0].get("filename"),
        preproc_variables_dict.get("E-OBS")[0].get("filename"),
    ]
    return era5_ancestors, eobs_ancestors


def get_model_output_filepath(dataset: str, value: list):
    """Generate output filepath for model data.

    Args:
    ----
        dataset : str
            dataset name
        value : dict
            Model variables

    Returns
    -------
        tuple(str, str)
            Output filepath and preprocessor path for
            future referencing.
    """
    timerange = value.get("timerange").replace("/", "-")
    experiment = value.get("exp")
    ensemble = value.get("ensemble")
    driver = value.get("driver", "")

    out_path = f"{dataset}/{driver}/{experiment}/{ensemble}/{timerange}"
    preproc_path = value.get("filename")

    return out_path, preproc_path


def get_preproc_lists(preproc_vars: list):
    """Put preprocessor variables and paths into list for further use.

    Args:
    ----
        preproc_vars : list
            List of variables for specific
            dataset

    Returns
    -------
        tuple(list, list)
            List of preprocessor cubes for mean calculations
            as well as path to files.
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
    """Put preprocessor variables and paths into list for further use.

    Args:
    ----
        preproc_vars : list
            Variable for specific ensemble member.

    Returns
    -------
        tuple(iris.cube.Cube, str)
            Preprocessor cube for mean calculations
            as well as path to files.
    """
    preproc_path = preproc_vars.get("filename")

    preproc = iris.load_cube(preproc_vars.get("filename"))

    # preproc_list = [mean_preproc]
    # preproc_path_list = [preproc_path_psl]

    return preproc, preproc_path


def get_looping_dict(preproc_vars: list):
    """Put cubes into dictionary for further use.

    Args:
    ----
        preproc_vars : list
            Values of dataset dictionary

    Returns
    -------
        Dictionary of the form {'var': [preproc_var, preproc_path]}
    """
    preproc, preproc_path = get_preproc_lists(preproc_vars)
    dict_ = {
        "psl": [preproc[0], preproc_path[0]],
        "pr": [preproc[1], preproc_path[1]],
        "tas": [preproc[2], preproc_path[2]],
    }
    return dict_


def load_wt_files(path: str, only_lwt=False):
    """Load *.nc files of weathertype data.

    If only_lwt is true, only Lamb
    weathertypes will be loaded. (useful for automatic_slwt = False)

    Args:
    ----
        path : str
            Path to weathertype data.
        only_lwt : bool, optional
            If True, only Lamb weathertypes will be loaded.
            Defaults to False.

    Returns
    -------
        list(iris.cube.Cube)
            List of weathertype cubes.
    """
    if not only_lwt:
        lwt_cube = iris.load_cube(path, "lwt")
        slwt_era5_cube = iris.load_cube(path, "slwt_era5")
        slwt_eobs_cube = iris.load_cube(path, "slwt_eobs")
        wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]
    else:
        lwt_cube = iris.load_cube(path, "lwt")
        wt_cubes = [lwt_cube]

    return wt_cubes


def get_provenance_record(
    caption: str,
    ancestors: list,
    long_names: list,
    plot_types: bool | list,
    statistics: bool | list,
) -> dict:
    """Get provenance record.

    Args:
    ----
        caption : str
            Caption for plots
        ancestors : list
            List of ancestors
        long_names : list
            Variable long names
        plot_types : bool | list
            Type of plot. Can be false if output is not a plot.
        statistics : bool | list
            Types of statistics used.

    Returns
    -------
        dict
            Provenance dictionary.
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
    if plot_types is not False:
        record["plot_types"] = plot_types
    if statistics is not False:
        record["statistics"] = statistics
    return record


def log_provenance(filename: str, cfg: dict, provenance_record: dict):
    """Log provenance. Produces xml file provenance info.

    Args:
    ----
        filename : str
            Output name of provenance.
        cfg : dict
            Configuration dictionary provided by recipe.
        provenance_record : dict
            Provenance record dictionary.
    """
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info("Output stored as %s", filename)


def turn_list_to_mapping_dict(list_: list) -> dict:
    """Turn list of combined WT into a dictionary for further processing.

    Args:
    ----
        list_ : list
            List where entries are lists with related WT

    Returns
    -------
        dict
            Mapping dictionary keys are simplified WT, values are Lamb WT
    """
    result_dict = {}

    for i, s in enumerate(list_):
        for elem in s:
            if elem not in result_dict:
                result_dict[elem] = i + 1
            else:
                result_dict[elem].append(i)

    return result_dict


def get_mapping_dict(selected_pairs: list) -> dict:
    """Get mapping dictionary from list of selected pairs.

    Args:
    ----
        selected_pairs : list
            List of selected pairs of related WT based on
            precipitation patterns over specified area and
            correlation and RSME thresholds defined in recipe

    Returns
    -------
        dict
            Mapping dicitonary keys are simplified WT, values are Lamb WT
    """
    mapping_array = []

    for elem in selected_pairs:
        mapping_array.append(elem[0])

    s = [set(i) for i in mapping_array if i]

    def find_intersection(m_list: list) -> list:
        for i, v in enumerate(m_list):
            for j, k in enumerate(m_list[i + 1 :], i + 1):
                if v & k:
                    s[i] = v.union(m_list.pop(j))
                    return find_intersection(m_list)
        return m_list

    merged_tuples = find_intersection(s)
    mapping_dict = turn_list_to_mapping_dict(merged_tuples)

    return mapping_dict


def write_mapping_dict(work_dir: str, dataset: str, mapping_dict: dict):
    """Write mapping dictionary to file.

    Args:
    ----
        work_dir : str
            Working directory to save mapping dict
        dataset : str
            Name of dataset
        mapping_dict : dict
            Mapping dictionary in {lwt: slwt, ...} format
    """
    mapping_dict_reformat = convert_dict(mapping_dict)

    with open(
        f"{work_dir}/wt_mapping_dict_{dataset}.json", "w", encoding="utf-8"
    ) as file:
        json.dump(mapping_dict_reformat, file)


def convert_dict(dict_: dict) -> dict:
    """Convert mapping dictionary.

    From {lwt: slwt, ...} format to {slwt: [lwt1, lwt2], ...}.

    Args:
    ----
        dict_ : dict
            Dict in the {lwt: slwt, ...} format

    Returns
    -------
        dict
            Dict in the {slwt: [lwt1, lwt2], ...} format
    """
    new_dict = {}
    for dataset, value in dict_.items():
        if value not in new_dict:
            new_dict[value] = []
        new_dict[value].append(dataset)
    return new_dict


def reverse_convert_dict(dict_: dict) -> dict:
    """Convert mapping dictionary.

    From {slwt: [lwt1, lwt2], ...} format to {lwt: slwt, ...}.

    Args:
    ----
        original_dict : dict
        Dict in the {slwt: [lwt1, lwt2], ...} format

    Returns
    -------
        dict
            Dict in the  format {lwt: slwt, ...}
    """
    new_dict = {}
    for key, value_list in dict_.items():
        for original_key in value_list:
            new_dict[original_key] = key
    return new_dict


def write_corr_rmse_to_csv(
    cfg: dict,
    pattern_correlation_matrix: np.array,
    rmse_matrix: np.array,
    dataset: str,
):
    """Write correlation and rsme matrix to csv files.

    Args:
    ----
        cfg : dict
            Configuration dictionary from recipe
        pattern_correlation_matrix : np.array
            Correlation matrix
        rmse_matrix : np.array
            RSME matrix
        dataset : str
            Name of dataset
    """
    logger.info("Writing corr and rsme matrices for %s", dataset)

    work_dir = cfg.get("work_dir")

    df_corr = pd.DataFrame(pattern_correlation_matrix)
    df_corr.index = range(1, len(df_corr) + 1)
    df_corr.columns = range(1, len(df_corr.columns) + 1)
    df_corr.to_csv(
        f"{work_dir}/correlation_matrix_{dataset}.csv", index_label="Index"
    )

    df_rmse = pd.DataFrame(rmse_matrix)
    df_rmse.index = range(1, len(df_rmse) + 1)
    df_rmse.columns = range(1, len(df_rmse.columns) + 1)
    df_rmse.to_csv(
        f"{work_dir}/rmse_matrix_{dataset}.csv", index_label="Index"
    )


def run_predefined_slwt(
    work_dir: str, dataset_name: str, lwt: np.array, predefined_slwt: dict
):
    """Run predefined slwt mapping.

    Args:
    ----
        work_dir : str
            Working directory to save mapping dict
        dataset_name : str
            Name of dataset
        lwt : np.array
            lwt array
        predefined_slwt : dict
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
    cfg: dict, wt_list: list, cube: iris.cube.Cube, file_name: str
):
    """Combine lwt and slwt arrays to one file.

    Args:
    ----
        cfg : dict
            Configuration dictionary from recipe
        wt_list : list
            List of weathertype arrays
        cube : iris.cube.Cube
            Cube of data to keep time coordinate
        file_name : str
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
    cfg: dict, lwt: np.array, cube: iris.cube.Cube, file_name: str
):
    """Write only lwt to file.

    Args:
    ----
        cfg : dict
            Configuration dictionary from recipe
        lwt : np.array
            lwt array
        cube : iris.cube.Cube
            Cube of data to keep time coordinate
        file_name : str
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

    Args:
    ----
        lwt : np.array
            lwt array
        mapping_dict : dict
            Mapping dictionary in {lwt: slwt, ...} format

    Returns
    -------
        np.array
            array of slwt
    """
    return np.array([np.int8(mapping_dict.get(value, 0)) for value in lwt])


def check_mapping_dict_format(mapping_dict: dict) -> dict:
    """Check format of mapping dict and return in {lwt: slwt, ...} format.

    Args:
    ----
        mapping_dict : dict
            mapping dict in any format

    Returns
    -------
        dict
            mapping dict in {lwt: slwt, ...} format
    """
    if isinstance(mapping_dict.get(list(mapping_dict.keys())[0]), list):
        dict_ = reverse_convert_dict(mapping_dict)
    else:
        dict_ = mapping_dict

    return dict_
