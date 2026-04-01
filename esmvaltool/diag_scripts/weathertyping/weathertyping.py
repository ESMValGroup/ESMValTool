"""Script for calculating Lamb weathertypes.

It also plots the means and seasonal occurrence of the
weathertypes, and offers the option to calculate simplified
weathertypes based on precipitation patterns.
"""

import iris

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.weathertyping.calc_utils import (
    calc_lwt_model,
    calc_lwt_slwt_model,
    calc_slwt_obs,
    plot_means,
    wt_algorithm,
)
from esmvaltool.diag_scripts.weathertyping.plot_utils import (
    plot_seasonal_occurrence,
)
from esmvaltool.diag_scripts.weathertyping.wt_utils import (
    combine_wt_to_file,
    get_ancestors_era5_eobs,
    get_looping_dict,
    get_model_output_filepath,
    get_preproc_lists_ensemble,
    get_provenance_record,
    load_wt_files,
    load_wt_preprocessors,
    log_provenance,
    run_predefined_slwt,
    write_lwt_to_file,
)


def process_models_automatic_slwt(
    cfg: dict,
    dataset_vars: list,
    data_info: dict,
):
    """Process model data for calculating Lamb and simplified weathertypes.

    Parameters
    ----------
    cfg
        Nested dictionary of metadata.
    dataset_vars
        List of variable dictionaries for a specific dataset.
    data_info
        Dictionary holding dataset information.
    """
    for ensemble_var in dataset_vars:
        if ensemble_var.get("preprocessor") == "weathertype_preproc":
            data_info["ensemble"] = ensemble_var.get("ensemble", "")
            data_info["driver"] = ensemble_var.get("driver", "")

            if data_info["driver"] != "":
                data_info["driver"] = "_" + data_info["driver"]

            wt_preproc = iris.load_cube(ensemble_var.get("filename"))

            output_file_path, preproc_path = get_model_output_filepath(
                data_info["dataset"],
                ensemble_var,
            )
            data_info["output_file_path"] = output_file_path
            data_info["preproc_path"] = preproc_path

            # calculate weathertypes
            calc_lwt_slwt_model(
                cfg,
                wt_preproc,
                data_info,
                cfg.get("predefined_slwt"),
            )
            # plot means
            if cfg.get("plotting", False):
                # load wt files
                wt_cube_path = (
                    f"{cfg.get('work_dir')}/{output_file_path}"
                    f"/{data_info['dataset']}"
                    f"{data_info['driver']}_"
                    f"{data_info['ensemble']}_"
                    f"{data_info['timerange']}.nc"
                )

                wt_cubes = load_wt_files(wt_cube_path)

                var_dict = {
                    f"{ensemble_var.get('short_name')}": get_preproc_lists_ensemble(
                        ensemble_var,
                    ),
                }
                for var_name, var_data in var_dict.items():
                    data_info["var"] = var_name
                    data_info["preproc_path"] = var_data[1]

                    plot_means(cfg, var_data[0], wt_cubes, data_info)
                plot_seasonal_occurrence(
                    cfg,
                    wt_cubes,
                    data_info,
                    wt_cube_path,
                )


def process_era5_automatic_slwt(
    data_info: dict,
    preproc_variables_dict: dict,
    cfg: dict,
    dataset_vars: list,
):
    """Process ERA5 data for calculating Lamb and simplified weathertypes.

    Parameters
    ----------
    data_info
        Dictionary holding dataset information.
    preproc_variables_dict
        Dictionary holding preprocessed variables for all datasets.
    cfg
        Nested dictionary of metadata.
    dataset_vars
        List of variable dictionaries for a specific dataset.
    """
    wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs = load_wt_preprocessors(
        data_info["dataset"],
        preproc_variables_dict,
    )

    # calculate lwt
    lwt = wt_algorithm(wt_preproc, data_info["dataset"])

    era5_ancestors, eobs_ancestors = get_ancestors_era5_eobs(
        data_info["dataset"],
        preproc_variables_dict,
    )

    # calculate simplified lwt based on precipitation
    # patterns or use predefined_slwt
    if not cfg.get("predefined_slwt"):
        slwt_era5 = calc_slwt_obs(
            cfg,
            lwt,
            wt_preproc_prcp,
            data_info["dataset"],
            era5_ancestors,
            data_info["timerange"],
        )
        slwt_eobs = calc_slwt_obs(
            cfg,
            lwt,
            wt_preproc_prcp_eobs,
            "E-OBS",
            eobs_ancestors,
            data_info["timerange"],
        )
    else:
        slwt_era5, slwt_eobs = run_predefined_slwt(
            cfg.get("work_dir"),
            data_info["dataset"],
            lwt,
            cfg.get("predefined_slwt"),
        )

    # write to file
    wt_list = [lwt, slwt_era5, slwt_eobs]
    combine_wt_to_file(cfg, wt_list, wt_preproc, data_info["dataset"])

    # load weathertype files as cubes
    wt_cube_path = f"{cfg.get('work_dir')}/{data_info['dataset']}.nc"
    wt_cubes = load_wt_files(wt_cube_path)

    if cfg.get("plotting", False):
        var_dict = get_looping_dict(
            dataset_vars,
        )  # dataset_vars is list of variables for dataset dataset_name
        # plot means
        for var_name, var_data in var_dict.items():
            data_info["var"] = var_name
            data_info["preproc_path"] = var_data[1]

            plot_means(cfg, var_data[0], wt_cubes, data_info)
        plot_seasonal_occurrence(cfg, wt_cubes, data_info, wt_cube_path)


def run_automatic_slwt(cfg: dict):
    """Run the automated calculation for simplified weathertypes and write to file, and plot the means and seasonal occurrence of the weathertypes.

    Parameters
    ----------
    cfg
        Nested dictionary of metadata
    """
    preproc_variables_dict = group_metadata(
        cfg.get("input_data").values(),
        "dataset",
    )
    for dataset_name, dataset_vars in preproc_variables_dict.items():
        data_info = {
            "timerange": dataset_vars[0].get("timerange").replace("/", "-"),
            "dataset": dataset_name,
        }
        if dataset_name == "ERA5":
            process_era5_automatic_slwt(
                data_info,
                preproc_variables_dict,
                cfg,
                dataset_vars,
            )
        else:
            if data_info["dataset"] == "E-OBS":
                continue
            process_models_automatic_slwt(cfg, dataset_vars, data_info)


def process_era5_lwt(
    preproc_variables_dict: dict,
    cfg: dict,
    dataset_vars: list,
    data_info: dict,
):
    """Process ERA5 data for calculating Lamb weathertypes.

    Parameters
    ----------
    preproc_variables_dict
        Dictionary holding preprocessed variables for all datasets.
    cfg
        Nested dictionary of metadata.
    dataset_vars
        List of variable dictionaries for a specific dataset.
    data_info
        Dictionary holding dataset information.
    """
    wt_preproc, _, _ = load_wt_preprocessors(
        data_info["dataset"],
        preproc_variables_dict,
    )

    # calculate lwt
    lwt = wt_algorithm(wt_preproc, data_info["dataset"])

    era5_ancestors, eobs_ancestors = get_ancestors_era5_eobs(
        data_info["dataset"],
        preproc_variables_dict,
    )

    ancestors = [era5_ancestors, eobs_ancestors]

    provenance_record = get_provenance_record(
        "Lamb Weathertypes",
        ancestors,
        ["Lamb Weathertypes"],
    )

    log_provenance(f"{cfg.get('work_dir')}/lwt_era5", cfg, provenance_record)

    # write only lwt to file
    write_lwt_to_file(cfg, lwt, wt_preproc, data_info["dataset"])

    # load wt files
    wt_cube_path = f"{cfg.get('work_dir')}/{data_info['dataset']}.nc"
    wt_cubes = load_wt_files(
        wt_cube_path,
        mode="lwt",
    )

    if cfg.get("plotting", False):
        var_dict = get_looping_dict(
            dataset_vars,
        )  # dataset_vars is list of variables for dataset dataset_name
        # plot means
        for var_name, var_data in var_dict.items():
            data_info["var"] = var_name
            data_info["preproc_path"] = var_data[1]

            plot_means(cfg, var_data[0], wt_cubes, data_info, "lwt")
        plot_seasonal_occurrence(cfg, wt_cubes, data_info, wt_cube_path)


def process_models_lwt(cfg: dict, dataset_vars: list, data_info: dict):
    """Process model data for calculating Lamb weathertypes.

    Parameters
    ----------
    cfg
        Nested dictionary of metadata.
    dataset_vars
        List of variable dictionaries for a specific dataset.
    data_info
        Dictionary holding dataset information.
    """
    for ensemble_var in dataset_vars:
        if ensemble_var.get("preprocessor") == "weathertype_preproc":
            data_info["ensemble"] = ensemble_var.get("ensemble", "")
            data_info["driver"] = ensemble_var.get("driver", "")

            if data_info["driver"] != "":
                data_info["driver"] = "_" + {data_info["driver"]}

            wt_preproc = iris.load_cube(dataset_vars[0].get("filename"))

            output_file_path, preproc_path = get_model_output_filepath(
                data_info["dataset"],
                dataset_vars,
            )

            data_info["output_file_path"] = output_file_path
            data_info["preproc_path"] = preproc_path

            # calculate weathertypes
            calc_lwt_model(cfg, wt_preproc, data_info)

            # load wt files
            wt_cube_path = (
                f"{cfg.get('work_dir')}/{output_file_path}"
                f"/{data_info['dataset']}"
                f"{data_info['driver']}"
                f"_{data_info['ensemble']}_"
                f"{data_info['timerange']}.nc"
            )

            wt_cubes = load_wt_files(
                wt_cube_path,
                mode="lwt",
            )

            var_dict = get_looping_dict(
                dataset_vars,
            )  # dataset_vars is list of variables for dataset_name

            if cfg.get("plotting", False):
                # plot means
                for var_name, var_data in var_dict.items():
                    data_info["var"] = var_name
                    data_info["preproc_path"] = var_data[1]

                    plot_means(
                        cfg,
                        var_data[0],
                        wt_cubes,
                        data_info,
                        "lwt",
                    )
                plot_seasonal_occurrence(
                    cfg,
                    wt_cubes,
                    data_info,
                    wt_cube_path,
                )


def run_lwt(cfg: dict):
    """Run calculation of weathertypes. Write to file, and plot the means of psl, tas, and pr for each weathertype. Plot seasonal occurrence of the weathertypes.

    Parameters
    ----------
    cfg
        Nested dictionary of metadata.
    """
    preproc_variables_dict = group_metadata(
        cfg.get("input_data").values(),
        "dataset",
    )
    for dataset_name, dataset_vars in preproc_variables_dict.items():
        data_info = {
            "timerange": dataset_vars[0].get("timerange").replace("/", "-"),
            "dataset": dataset_name,
        }

        if dataset_name == "ERA5":
            process_era5_lwt(
                preproc_variables_dict,
                cfg,
                dataset_vars,
                data_info,
            )
        else:
            if data_info["dataset"] == "E-OBS":
                continue
            process_models_lwt(cfg, dataset_vars, data_info)


def run_weathertyping_diagnostic(cfg: dict):
    """Run the weathertyping diagnostic.

    Parameters
    ----------
    cfg
        Nested dictionary of metadata
    """
    automatic_slwt = cfg.get("automatic_slwt")

    # check if user wants to calculate simplified weathertypes automatically
    # for that, automatic_slwt must be true
    # additionally, if a predefined_slwt is given, those will be used
    if automatic_slwt:
        run_automatic_slwt(cfg)
    # if automatic_slwt is false, and predefined_slwt is false,
    # just write selected pairs to file
    else:
        run_lwt(cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        # main function for running the diagnostic
        run_weathertyping_diagnostic(config)
