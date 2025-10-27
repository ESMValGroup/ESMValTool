"""Script for calculating Lamb weathertypes.

It also plots the means and seasonal occurrence of the
weathertypes, and offers the option to calculate simplified
weathertypes based on precipitation patterns.
"""

import iris

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.weathertyping.calc_utils import (
    calc_lwt_model,
    calc_lwt_slwt_model,
    calc_slwt_obs,
    wt_algorithm,
)
from esmvaltool.diag_scripts.weathertyping.plot_utils import (
    plot_means,
    plot_seasonal_occurrence,
)
from esmvaltool.diag_scripts.weathertyping.wt_utils import (
    add_dict_entry,
    combine_wt_to_file,
    get_ancestors_era5_eobs,
    get_cfg_vars,
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


def run_automatic_slwt(cfg: dict):
    """Run the automated calculation for simplified weathertypes \
    and write to file, and plot the means and seasonal occurrence \
    of the weathertypes.

    Args:
        cfg (dict): Nested dictionary of metadata
    """
    preproc_variables_dict, _, _, work_dir, plotting, _, predefined_slwt = (
        get_cfg_vars(cfg)
    )
    for dataset_name, dataset_vars in preproc_variables_dict.items():
        data_info = {
            "timerange": dataset_vars[0].get("timerange").replace("/", "-"),
            "dataset": dataset_name,
        }
        if dataset_name == "ERA5":
            wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs = (
                load_wt_preprocessors(dataset_name, preproc_variables_dict)
            )

            # calculate lwt
            lwt = wt_algorithm(wt_preproc, dataset_name)

            era5_ancestors, eobs_ancestors = get_ancestors_era5_eobs(
                dataset_name, preproc_variables_dict
            )

            # calculate simplified lwt based on precipitation
            # patterns or use predefined_slwt
            if not predefined_slwt:
                slwt_era5 = calc_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp,
                    dataset_name,
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
                    work_dir, dataset_name, lwt, predefined_slwt
                )

            # write to file
            wt_list = [lwt, slwt_era5, slwt_eobs]
            combine_wt_to_file(cfg, wt_list, wt_preproc, dataset_name)

            # load weathertype files as cubes
            wt_cubes = load_wt_files(f"{work_dir}/{dataset_name}.nc")

            var_dict = get_looping_dict(
                dataset_vars
            )  # dataset_vars is list of variables for dataset dataset_name
            if plotting:
                # plot means
                for var_name, var_data in var_dict.items():
                    add_dict_entry(data_info, "var", var_name)
                    add_dict_entry(data_info, "preproc_path", var_data[1])

                    plot_means(cfg, var_data[0], wt_cubes, data_info)
                plot_seasonal_occurrence(cfg, wt_cubes, data_info)
        else:
            if dataset_name == "E-OBS":
                continue
            for ensemble_var in dataset_vars:
                if ensemble_var.get("preprocessor") == "weathertype_preproc":
                    add_dict_entry(
                        data_info, "ensemble", ensemble_var.get("ensemble", "")
                    )
                    add_dict_entry(
                        data_info, "driver", ensemble_var.get("driver", "")
                    )
                    if data_info["driver"] != "":
                        data_info["driver"] = "_" + {data_info["driver"]}

                    wt_preproc = iris.load_cube(ensemble_var.get("filename"))

                    output_file_path, preproc_path = get_model_output_filepath(
                        dataset_name, ensemble_var
                    )

                    data_info["output_file_path"] = output_file_path
                    data_info["preproc_path"] = preproc_path

                    # calculate weathertypes
                    calc_lwt_slwt_model(
                        cfg, wt_preproc, data_info, predefined_slwt
                    )

                    # load wt files
                    wt_cubes = load_wt_files(
                        f"{work_dir}/{output_file_path}"
                        f"/{dataset_name}"
                        f"{data_info['driver']}_"
                        f"{data_info['ensemble']}_"
                        f"{data_info['timerange']}.nc"
                    )

                    var_dict = {
                        f"{ensemble_var.get('short_name')}": get_preproc_lists_ensemble(
                            ensemble_var
                        )
                    }

                    # plot means
                    if plotting:
                        for var_name, var_data in var_dict.items():
                            add_dict_entry(data_info, "var", var_name)
                            add_dict_entry(
                                data_info, "preproc_path", var_data[1]
                            )

                            plot_means(cfg, var_data[0], wt_cubes, data_info)
                        plot_seasonal_occurrence(cfg, wt_cubes, data_info)


def run_lwt(cfg: dict):
    """Run calculation of weathertypes and write to file, and plot the means \
    and seasonal occurrence of the weathertypes.

    Args:
        cfg (dict): Nested dictionary of metadata

    """
    preproc_variables_dict, _, _, work_dir, plotting, _, _ = get_cfg_vars(cfg)
    for dataset_name, dataset_vars in preproc_variables_dict.items():
        data_info = {
            "timerange": dataset_vars[0].get("timerange").replace("/", "-"),
            "dataset": dataset_name,
        }

        if dataset_name == "ERA5":
            wt_preproc, _, _ = load_wt_preprocessors(
                dataset_name, preproc_variables_dict
            )

            # calculate lwt
            lwt = wt_algorithm(wt_preproc, dataset_name)

            era5_ancestors, eobs_ancestors = get_ancestors_era5_eobs(
                dataset_name, preproc_variables_dict
            )

            ancestors = [era5_ancestors, eobs_ancestors]

            provenance_record = get_provenance_record(
                "Lamb Weathertypes",
                ancestors,
                ["Lamb Weathertypes"],
                False,
                False,
            )

            log_provenance(f"{work_dir}/lwt_era5", cfg, provenance_record)

            # write only lwt to file
            write_lwt_to_file(cfg, lwt, wt_preproc, dataset_name)

            # load wt files
            wt_cubes = load_wt_files(
                f"{work_dir}/{dataset_name}.nc", only_lwt=True
            )

            var_dict = get_looping_dict(
                dataset_vars
            )  # dataset_vars is list of variables for dataset dataset_name

            if plotting:
                # plot means
                for var_name, var_data in var_dict.items():
                    add_dict_entry(data_info, "var", var_name)
                    add_dict_entry(data_info, "preproc_path", var_data[1])

                    plot_means(
                        cfg, var_data[0], wt_cubes, data_info, only_lwt=True
                    )
                plot_seasonal_occurrence(cfg, wt_cubes, data_info)
        else:
            if dataset_name == "E-OBS":
                continue
            for ensemble_var in dataset_vars:
                if ensemble_var.get("preprocessor") == "weathertype_preproc":
                    add_dict_entry(
                        data_info, "ensemble", ensemble_var.get("ensemble", "")
                    )
                    add_dict_entry(
                        data_info, "driver", ensemble_var.get("driver", "")
                    )
                    if data_info["driver"] != "":
                        data_info["driver"] = "_" + {data_info["driver"]}

                    wt_preproc = iris.load_cube(
                        dataset_vars[0].get("filename")
                    )

                    output_file_path, preproc_path = get_model_output_filepath(
                        dataset_name, dataset_vars
                    )

                    # calculate weathertypes
                    data_info["output_file_path"] = output_file_path
                    data_info["preproc_path"] = preproc_path

                    calc_lwt_model(cfg, wt_preproc, dataset_name, data_info)

                    # load wt files
                    wt_cubes = load_wt_files(
                        f"{work_dir}/{output_file_path}"
                        f"/{dataset_name}"
                        f"{data_info['driver']}"
                        f"_{data_info['ensemble']}_"
                        f"{data_info['timerange']}.nc",
                        only_lwt=True,
                    )

                    var_dict = get_looping_dict(
                        dataset_vars
                    )  # dataset_vars is list of variables for dataset_name

                    if plotting:
                        # plot means
                        for var_name, var_data in var_dict.items():
                            add_dict_entry(data_info, "var", var_name)
                            add_dict_entry(
                                data_info, "preproc_path", var_data[1]
                            )

                            plot_means(
                                cfg,
                                var_data[0],
                                wt_cubes,
                                data_info,
                                only_lwt=True,
                            )
                        plot_seasonal_occurrence(cfg, wt_cubes, data_info)


def run_my_diagnostic(cfg: dict):
    """Run the weathertyping diagnostic.

    Arguments:
        cfg - nested dictionary of metadata

    Returns
        string; runs the user diagnostic

    """
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable
    automatic_slwt = cfg.get("automatic_slwt")

    if automatic_slwt:
        run_automatic_slwt(cfg)
    # if automatic_slwt is false, and predefined_slwt is false,
    # just write selected pairs to file
    elif not automatic_slwt:
        run_lwt(cfg)


if __name__ == "__main__":
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
