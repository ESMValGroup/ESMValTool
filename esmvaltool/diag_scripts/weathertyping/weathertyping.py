from utils import *

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    ProvenanceLogger
)


def run_automatic_slwt(cfg: dict, preproc_variables_dict: dict,
                       correlation_threshold: float, 
                       rmse_threshold: float, work_dir: str,
                       plotting: bool, predefined_slwt):
    for key, value in preproc_variables_dict.items():
            if key == 'ERA5':
                wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs = load_wt_preprocessors(
                    value, preproc_variables_dict
                )

                # calculate lwt
                lwt = wt_algorithm(wt_preproc, key)

                era5_ancestors, eobs_ancestors = get_ancestors(value, preproc_variables_dict)

                # calculate simplified lwt based on precipitation patterns or use predefined_slwt
                if predefined_slwt is False:
                    slwt_era5 = calc_slwt_obs(
                        cfg,
                        lwt,
                        wt_preproc_prcp,
                        key,
                        correlation_threshold,
                        rmse_threshold,
                        era5_ancestors,
                    )
                    slwt_eobs = calc_slwt_obs(
                        cfg,
                        lwt,
                        wt_preproc_prcp_eobs,
                        'E-OBS',
                        correlation_threshold,
                        rmse_threshold,
                        eobs_ancestors,
                    )
                else:
                    predefined_slwt = check_mapping_dict_format(
                        predefined_slwt)
                    write_mapping_dict(work_dir, key, predefined_slwt)
                    write_mapping_dict(work_dir, "E-OBS", predefined_slwt)
                    slwt_era5 = map_lwt_to_slwt(lwt, predefined_slwt)
                    slwt_eobs = map_lwt_to_slwt(lwt, predefined_slwt)

                # write to file
                combine_wt_to_file(cfg, lwt, slwt_era5,
                                   slwt_eobs, wt_preproc, key)

                # load weathertype files as cubes
                wt_cubes = load_wt_files(f'{work_dir}/{key}.nc')

                var_dict = get_looping_dict(value) #value is list of variables for dataset key                     

                if plotting:
                    # plot means
                    for var_name, var_data in var_dict.items():
                        print(var_name, var_data[1])
                        print(var_name, var_data[0])
                        plot_means(cfg, var_data[0], wt_cubes, key, var_name, var_data[1])
            else:
                if key == 'E-OBS':
                    continue
                wt_preproc = iris.load_cube(value[0].get('filename'))

                output_file_path, preproc_path = get_model_output_filepath(key, value)

                # calculate weathertypes
                calc_lwt_slwt_model(
                    cfg, wt_preproc, key, preproc_path, output_file_path,
                    predefined_slwt)

                # load wt files
                wt_cubes = load_wt_files(f'{work_dir}/{output_file_path}/{key}.nc')

                var_dict = get_looping_dict(value) #value is list of variables for dataset key                    

                if plotting:
                    # plot means
                    for var_name, var_data in var_dict.items():
                        plot_means(cfg, var_data[0], wt_cubes, key, var_name, var_data[1])


def run_lwt(cfg: dict, preproc_variables_dict: dict,
                       correlation_threshold: float, 
                       rmse_threshold: float, work_dir: str,
                       plotting: bool):
    for key, value in preproc_variables_dict.items():
            if key == 'ERA5':
                wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs = load_wt_preprocessors(
                    value, preproc_variables_dict
                )

                # calculate lwt
                lwt = wt_algorithm(wt_preproc, key)

                era5_ancestors, eobs_ancestors = get_ancestors(value, preproc_variables_dict)

                # calculate simplified lwt based on precipitation patterns
                slwt_era5 = calc_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp,
                    key,
                    correlation_threshold,
                    rmse_threshold,
                    era5_ancestors,
                )
                slwt_eobs = calc_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp_eobs,
                    'E-OBS',
                    correlation_threshold,
                    rmse_threshold,
                    eobs_ancestors,
                )

                # write only lwt to file
                write_lwt_to_file(cfg, lwt, wt_preproc, key)

                # load wt files
                wt_cubes = load_wt_files(f'{work_dir}/{key}.nc', only_lwt=True)

                var_dict = get_looping_dict(value) #value is list of variables for dataset key

                if plotting:
                    # plot means
                    for var_name, var_data in var_dict.items():
                        plot_means(cfg, var_data[0], wt_cubes, key, var_name, var_data[1], only_lwt=True)
            else:
                if key == 'E-OBS':
                    continue
                wt_preproc = iris.load_cube(value[0].get('filename'))

                output_file_path = get_model_output_filepath(key, value)[1]

                # calculate weathertypes
                calc_lwt_model(cfg, wt_preproc, key, output_file_path)

                # load wt files
                wt_cubes = load_wt_files(f'{work_dir}/{key}.nc', only_lwt=True)

                var_dict = get_looping_dict(value) #value is list of variables for dataset key                     

                if plotting:
                    # plot means
                    for var_name, var_data in var_dict.items():
                        plot_means(cfg, var_data[0], wt_cubes, key, var_name, var_data[1], only_lwt=True)



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
    preproc_variables_dict = group_metadata(
        cfg.get('input_data').values(), 'dataset')

    correlation_threshold = cfg.get('correlation_threshold')
    rmse_threshold = cfg.get('rmse_threshold')
    work_dir = cfg.get('work_dir')
    plotting = cfg.get('plotting', False)
    automatic_slwt = cfg.get('automatic_slwt')
    predefined_slwt = cfg.get('predefined_slwt', False)

    # load cubes and run functions
    # key = dataset name, value is dataset
    if automatic_slwt:
        run_automatic_slwt(cfg, preproc_variables_dict,
                       correlation_threshold, 
                       rmse_threshold, work_dir,
                       plotting, predefined_slwt)
    elif not automatic_slwt:  # if automatic_slwt is false, and predefined_slwt is false, just write selected pairs to file
        run_lwt(cfg, preproc_variables_dict,
                       correlation_threshold, 
                       rmse_threshold, work_dir,
                       plotting)

if __name__ == '__main__':

    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
