""" This script calculates weathertypes for the input datasets and writes \
    to file. It also plots the means and seasonal occurrence of the
    weathertypes, and offers the option to calculate simplified
    weathertypes based on precipitation patterns."""
import iris
from wt_utils import get_cfg_vars, load_wt_preprocessors, wt_algorithm, \
    get_ancestors_era5_eobs, calc_slwt_obs, run_predefined_slwt, \
    combine_wt_to_file, load_wt_files, get_looping_dict, plot_means, \
    plot_seasonal_occurrence, write_lwt_to_file, get_model_output_filepath, \
    calc_lwt_slwt_model, calc_lwt_model

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic


def run_automatic_slwt(cfg: dict):
    """Run the automated calculation for simplified weathertypes \
    and write to file, and plot the means and seasonal occurrence \
    of the weathertypes.

    Args:
        cfg (dict): Nested dictionary of metadata
    """
    preproc_variables_dict, _, _, \
        work_dir, plotting, _, predefined_slwt = get_cfg_vars(cfg)
    for dataset_name, dataset_vars in preproc_variables_dict.items():
        if dataset_name == 'ERA5':
            wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs = \
                load_wt_preprocessors(dataset_name, preproc_variables_dict)

            # calculate lwt
            lwt = wt_algorithm(wt_preproc, dataset_name)

            era5_ancestors, eobs_ancestors = get_ancestors_era5_eobs(
                dataset_name, preproc_variables_dict)

            # calculate simplified lwt based on precipitation
            # patterns or use predefined_slwt
            if predefined_slwt is False:
                slwt_era5 = calc_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp,
                    dataset_name,
                    era5_ancestors,
                )
                slwt_eobs = calc_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp_eobs,
                    'E-OBS',
                    eobs_ancestors,
                )
            else:
                slwt_era5, slwt_eobs = run_predefined_slwt(work_dir,
                                                           dataset_name,
                                                           lwt,
                                                           predefined_slwt)

            # write to file
            wt_list = [lwt, slwt_era5, slwt_eobs]
            combine_wt_to_file(cfg, wt_list, wt_preproc,
                               dataset_name)

            # load weathertype files as cubes
            wt_cubes = load_wt_files(f'{work_dir}/{dataset_name}.nc')

            var_dict = get_looping_dict(
                dataset_vars
            )  # dataset_vars is list of variables for dataset dataset_name
            if plotting:
                # plot means
                for var_name, var_data in var_dict.items():
                    data_info = {'dataset': dataset_name,
                                 'var': var_name,
                                 'preproc_path': var_data[1],
                                 'ensemble':
                                     dataset_vars[0].get('ensemble', ''),
                                 'timerange': dataset_vars[0].
                                     get('timerange').replace('/', '-')}
                    plot_means(cfg, var_data[0], wt_cubes, data_info)
                plot_seasonal_occurrence(cfg, wt_cubes, dataset_name)
        else:
            if dataset_name == 'E-OBS':
                continue
            wt_preproc = iris.load_cube(dataset_vars[0].get('filename'))

            output_file_path, preproc_path = get_model_output_filepath(
                dataset_name, dataset_vars)

            # calculate weathertypes
            data_info = {'dataset': dataset_name,
                         'preproc_path': preproc_path,
                         'output_file_path': output_file_path,
                         'ensemble': dataset_vars[0].get('ensemble', ''),
                         'timerange': dataset_vars[0].get(
                             'timerange').replace('/', '-')}
            calc_lwt_slwt_model(cfg, wt_preproc, data_info, predefined_slwt)

            # load wt files
            wt_cubes = load_wt_files(
                f'{work_dir}/{output_file_path}/{dataset_name}.nc')

            var_dict = get_looping_dict(
                dataset_vars
            )  # dataset_vars is list of variables for dataset dataset_name

            # plot means
            if plotting:
                for var_name, var_data in var_dict.items():
                    data_info = {'dataset': dataset_name,
                                 'var': var_name,
                                 'preproc_path': var_data[1],
                                 'ensemble': dataset_vars[0].get(
                                     'ensemble', ''),
                                 'timerange':
                                     dataset_vars[0].
                                     get('timerange').replace('/', '-')}
                    plot_means(cfg, var_data[0], wt_cubes, data_info,
                               only_lwt=True)
                plot_seasonal_occurrence(cfg, wt_cubes, dataset_name,
                                         ensemble=dataset_vars[0].
                                         get('ensemble', ''))


def run_lwt(cfg: dict):
    """Run calculation of weathertypes and write to file, and plot the means \
    and seasonal occurrence of the weathertypes.

    Args:
        cfg (dict): Nested dictionary of metadata
    """
    preproc_variables_dict, _, _, \
        work_dir, plotting, _, _ = get_cfg_vars(cfg)
    for dataset_name, dataset_vars in preproc_variables_dict.items():
        if dataset_name == 'ERA5':
            wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs = \
                load_wt_preprocessors(dataset_name, preproc_variables_dict)

            # calculate lwt
            lwt = wt_algorithm(wt_preproc, dataset_name)

            era5_ancestors, eobs_ancestors = get_ancestors_era5_eobs(
                dataset_name, preproc_variables_dict)

            # calculate simplified lwt based on precipitation patterns
            _ = calc_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp,
                    dataset_name,
                    era5_ancestors,
                )
            _ = calc_slwt_obs(
                cfg,
                lwt,
                wt_preproc_prcp_eobs,
                'E-OBS',
                eobs_ancestors,
            )

            # write only lwt to file
            write_lwt_to_file(cfg, lwt, wt_preproc, dataset_name)

            # load wt files
            wt_cubes = load_wt_files(f'{work_dir}/{dataset_name}.nc',
                                     only_lwt=True)

            var_dict = get_looping_dict(
                dataset_vars
            )  # dataset_vars is list of variables for dataset dataset_name

            if plotting:
                # plot means
                for var_name, var_data in var_dict.items():
                    data_info = {'dataset': dataset_name,
                                 'var': var_name,
                                 'preproc_path': var_data[1],
                                 'ensemble': dataset_vars[0].get(
                                     'ensemble', ''),
                                 'timerange': dataset_vars[0].get(
                                     'timerange').replace('/', '-')}
                    plot_means(cfg, var_data[0], wt_cubes, data_info,
                               only_lwt=True)
                plot_seasonal_occurrence(cfg, wt_cubes, dataset_name)
        else:
            if dataset_name == 'E-OBS':
                continue
            wt_preproc = iris.load_cube(dataset_vars[0].get('filename'))

            output_file_path = get_model_output_filepath(
                dataset_name, dataset_vars)[1]

            # calculate weathertypes
            calc_lwt_model(cfg, wt_preproc, dataset_name, output_file_path)

            # load wt files
            wt_cubes = load_wt_files(f'{work_dir}/{dataset_name}.nc',
                                     only_lwt=True)

            var_dict = get_looping_dict(
                dataset_vars
            )  # dataset_vars is list of variables for dataset dataset_name

            if plotting:
                # plot means
                ensemble = dataset_vars[0].get('ensemble')
                for var_name, var_data in var_dict.items():
                    data_info = {'dataset': dataset_name,
                                 'var': var_name,
                                 'preproc_path': var_data[1],
                                 'ensemble': dataset_vars[0].get(
                                     'ensemble', ''),
                                 'timerange': dataset_vars[0].get(
                                     'timerange').replace('/', '-')}
                    plot_means(cfg, var_data[0], wt_cubes, data_info,
                               only_lwt=True)
                plot_seasonal_occurrence(cfg, wt_cubes, dataset_name,
                                         ensemble=ensemble)


def run_my_diagnostic(cfg: dict):
    '''
    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    '''
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable
    automatic_slwt = cfg.get('automatic_slwt')

    if automatic_slwt:
        run_automatic_slwt(cfg)
    # if automatic_slwt is false, and predefined_slwt is false,
    # just write selected pairs to file
    elif not automatic_slwt:
        run_lwt(cfg)


if __name__ == '__main__':

    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
