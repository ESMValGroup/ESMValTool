"""
Arctic-midlatitude teleconnections Diagnostics.

Diagnostic calculates timeseries needed for Causal Model Evaluation
of Arctic-midlatitude teleconnections.

  Description:
    This diagnostics calculates timeseries of the variables that represent
    Arctic-midlatitude teleconenctions that are further used for the
    Causal Model Evaluation of CMIP6 (Galytska et al., 2023). The output of
    this diagnostics is a .nc file per data source. Optionally this diagnostics
    plots the timeseries of the evolution of each selected variable. If the
    user kept "plot_timeseries: True" in recipe_galytska23jgr.yml, then
    "variable_to_plot:" expects the name of the variable to be plotted.
    Possible options for "variable_to_plot:" are:
    Arctic_temperature
    Psl_Ural
    Psl_Sib
    Psl_Aleut
    PV
    heat_flux
    BK_sic
    Ok_sic
  Author: Evgenia Galytska, IUP-UB
          egalytska@iup.physik.uni-bremen.de
  Project: USMILE
"""

import logging
from pathlib import Path

import iris
import numpy as np
import seaborn as sns
from esmvalcore.preprocessor import (
    anomalies,
    area_statistics,
    meridional_statistics,
    zonal_statistics,
)
from matplotlib import pyplot as plt

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
)
from esmvaltool.diag_scripts.shared._base import (
    get_plot_filename,
)

logger = logging.getLogger(Path(__file__).stem)

# Fixed parameters
# list of variables to be ignored per model
ignored_variables = {"HadISST": ["heat_flux"]}

# list of variables per dataset that will be processed
proc_vars = {
    "ERA5": [
        "PV",
        "Arctic_temperature",
        "Psl_Ural",
        "Psl_Sib",
        "Psl_Aleut",
        "heat_flux",
    ],
    "HadISST": ["BK_sic", "Ok_sic"],
    "all_other_datasets": [
        "PV",
        "Arctic_temperature",
        "Psl_Ural",
        "Psl_Sib",
        "Psl_Aleut",
        "heat_flux",
        "BK_sic",
        "Ok_sic",
    ],
}


def get_provenance_record(ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "authors": ["galytska_evgenia"],
        "ancestors": ancestor_files,
        "projects": ["usmile"],
        "references": [
            "galytska23jgr",
        ],
    }
    return record


def calculate_polar_vortex(dict_item):
    """Calculate polar vortex."""
    var = iris.load_cube(dict_item["filename"])
    var = var.collapsed("air_pressure", iris.analysis.MEAN)
    # change the sign of polar vortex so the positive values
    # (negative geopotential height anomalies) stand for
    # the strong polar vortex, similarly to
    # Kretschmer et al., 2016 and Galytska et al., 2023
    var.data *= -1
    var.var_name = "PV"
    return var


def calculate_arctic_tas(dict_item):
    """Read Arctic temperature data."""
    var = iris.load_cube(dict_item["filename"])
    var.var_name = "Arctic_temperature"
    return var


def calculate_slp(dict_item):
    """Get surface pressure."""
    var = iris.load_cube(dict_item["filename"])
    # calculate hPa from Pa.
    var.data /= 100
    return var


def finalize_bk_ice(dict_item):
    """Read sea ice data (Barents-Kara seas)."""
    var = iris.load_cube(dict_item["filename"])
    var.var_name = "BK_sic"
    return var


def finalize_ok_ice(dict_item):
    """Read sea ice data (Sea of Okhotsk)."""
    var = iris.load_cube(dict_item["filename"])
    var.var_name = "Ok_sic"
    return var


def prepare_heat_flux(dict_item):
    """Prepare variables for the heat flux calculations."""
    var = iris.load_cube(dict_item["filename"])
    var_avg = area_statistics(var, operator="mean")
    var_mermean = meridional_statistics(var, operator="mean")
    deviation = var_mermean - var_avg
    return deviation


def calculate_heat_flux(list_va_ta):
    """Calculate eddy poleward heat flux."""
    heat_flux = list_va_ta[0] * list_va_ta[1]
    hf_anom = anomalies(heat_flux, period="monthly")
    hf_anom_zm = zonal_statistics(hf_anom, operator="mean")
    hf_anom_zm.var_name = "heat_flux"
    return hf_anom_zm


def variable_cases(var_name, var):
    """Match preprocessor name and corresponding calculations."""
    if var_name == "pv":
        out_var = calculate_polar_vortex(var)
    elif var_name == "pre_tas":
        out_var = calculate_arctic_tas(var)
    elif var_name == "pressure_ural":
        out_var = calculate_slp(var)
        out_var.var_name = "Psl_Ural"
    elif var_name == "pressure_sib":
        out_var = calculate_slp(var)
        out_var.var_name = "Psl_Sib"
    elif var_name == "pressure_aleut":
        out_var = calculate_slp(var)
        out_var.var_name = "Psl_Aleut"
    elif var_name == "bk_ice":
        out_var = finalize_bk_ice(var)
    elif var_name == "ok_ice":
        out_var = finalize_ok_ice(var)
    elif var_name == "heat_flux":
        out_var = prepare_heat_flux(var)
    else:
        raise NotImplementedError(f"Variable '{var_name}' not yet supported.")
    return out_var


def calculate_variables(dataset_dict):
    """Calculate all necessary variables."""
    logger.debug(
        "Variables are calculated for the following datasources:%s",
        dataset_dict.keys(),
    )
    processed_vars = {}
    for dataset, variables in dataset_dict.items():
        processed_vars[dataset] = {}

        logger.debug(
            "Calculating final variables %s for %s dataset",
            variables,
            dataset,
        )

        if dataset in ignored_variables:
            to_ignore_vars = ignored_variables.get(dataset)
            for var in variables:
                var_name = var["preprocessor"]
                if var_name not in to_ignore_vars:
                    new_var = variable_cases(var_name, var)
                    new_var_name = new_var.var_name
                    processed_vars[dataset][new_var_name] = new_var
        else:
            tmp_list = []
            for var in variables:
                var_name = var["preprocessor"]
                if var_name == "heat_flux":
                    tmp_list.append(variable_cases(var_name, var))
                else:
                    new_var = variable_cases(var_name, var)
                    new_var_name = new_var.var_name
                    processed_vars[dataset][new_var_name] = new_var
            if len(tmp_list) != 2:
                raise IndexError(
                    "The preprocessor heat flux requests two \
                                  variables in the recipe: va and ta",
                )
            heat_flux = calculate_heat_flux(tmp_list)
            processed_vars[dataset][heat_flux.var_name] = heat_flux

    return processed_vars


def plotting_support(cube, key, **kwargs):
    """Help for the pretty plot."""
    if cube.coords("time", dim_coords=True):
        ih.unify_time_coord(cube)
    iris.quickplot.plot(cube, label=key, **kwargs)
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.ylabel("Anomalies, " + str(cube.units))
    plt.title(f"Time series of monthly mean {cube.var_name.upper()} anomalies")
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")


def plot_timeseries(dictionary, var, cfg):
    """Timeseries plot."""
    fig = plt.figure(figsize=(10, 4))
    sns.set_style("whitegrid")
    colors = plt.cm.viridis(np.linspace(0, 1, len(dictionary.keys())))
    baseplotname = f"Timeseries_{var}_anomalies"
    filename = get_plot_filename(baseplotname, cfg)
    for idx, dataset in enumerate(dictionary.keys()):
        if var not in proc_vars["HadISST"]:
            if dataset == "HadISST":
                continue
            if dataset != "ERA5":
                plotting_support(
                    dictionary[dataset][var],
                    dataset,
                    color=colors[idx],
                )
            else:
                plotting_support(
                    dictionary[dataset][var],
                    dataset,
                    color="k",
                    linewidth=2,
                )
        else:
            if dataset == "ERA5":
                continue
            if dataset != "HadISST":
                plotting_support(
                    dictionary[dataset][var],
                    dataset,
                    color=colors[idx],
                )
            else:
                plotting_support(
                    dictionary[dataset][var],
                    dataset,
                    color="blue",
                    linewidth=2,
                )
    fig.savefig(filename, bbox_inches="tight")


def assemble_cube_list(dataset, var, special_datasets):
    """
    Assemble a list of processed vars cubes.

    Depending on what vars are needed per dataset,
    variables list differs per analyzed dataset. Dict holding the
    needed variables per dataset needs updating everytime a new dataset
    or variable gets included.

    Parameters
    ----------
    dataset: str
        dataset name.
    var: dict
        variable dictionary.
    special_datasets: list
        list of datasets to be treated separately,
        with restricted variables.
        type: list of datasets (list of strings).

    Returns
    -------
    iris.cube.CubeList
        list of cubes.
    """
    if dataset not in special_datasets:
        cube_list = iris.cube.CubeList(
            [var[proc_var] for proc_var in proc_vars["all_other_datasets"]],
        )
    else:
        cube_list = iris.cube.CubeList(
            [var[proc_var] for proc_var in proc_vars[dataset]],
        )

    return cube_list


def main(cfg):
    """Calculate and save final variables into .nc files."""
    special_datasets = ["ERA5", "HadISST"]

    my_files_dict = group_metadata(cfg["input_data"].values(), "dataset")
    all_variables = calculate_variables(my_files_dict)

    # Check is timeseries should be plotted
    if cfg["plot_timeseries"]:
        plot_timeseries(all_variables, cfg["variable_to_plot"], cfg)
    for dataset in my_files_dict:
        logger.info("Processing final calculations in dataset %s", dataset)
        prov_record = get_provenance_record([dataset])
        var = all_variables[dataset]
        cube_list = assemble_cube_list(dataset, var, special_datasets)
        save_data(dataset, prov_record, cfg, cube_list)
        logger.info("%s data is saved in .nc", dataset)
    logger.info("Done.")


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
