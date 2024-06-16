from utils import *

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    ProvenanceLogger
)


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
        for key, value in preproc_variables_dict.items():
            if key == 'ERA5':
                wt_preproc = iris.load_cube(value[0].get('filename'))
                wt_preproc_prcp = iris.load_cube(value[1].get('filename'))
                mean_preproc_psl = iris.load_cube(value[2].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[3].get('filename'))
                mean_preproc_tas = iris.load_cube(value[4].get('filename'))
                wt_preproc_prcp_eobs = iris.load_cube(
                    preproc_variables_dict.get('E-OBS')[0].get('filename')
                )

                # calculate lwt
                lwt = wt_algorithm(wt_preproc, key)

                era5_ancestors = [value[0].get('filename'),
                                  value[1].get('filename')]
                eobs_ancestors = [value[0].get('filename'),
                                  preproc_variables_dict.get('E-OBS')[0].get('filename')]

                # calculate simplified lwt based on precipitation patterns or use predefined_slwt
                if predefined_slwt is False:
                    slwt_era5 = calculate_slwt_obs(
                        cfg,
                        lwt,
                        wt_preproc_prcp,
                        key,
                        correlation_threshold,
                        rmse_threshold,
                        era5_ancestors,
                    )
                    slwt_eobs = calculate_slwt_obs(
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
                lwt_cube = iris.load_cube(f'{work_dir}/{key}.nc', 'lwt')
                slwt_era5_cube = iris.load_cube(
                    f'{work_dir}/{key}.nc', 'slwt_era5')
                slwt_eobs_cube = iris.load_cube(
                    f'{work_dir}/{key}.nc', 'slwt_eobs')
                wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]

                preproc_path_psl = value[2].get('filename')
                preproc_path_prcp = value[3].get('filename')
                preproc_path_tas = value[4].get('filename')

                if plotting:
                    # plot means
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='lwt',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='lwt',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='lwt',
                        preproc_path=preproc_path_tas,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='slwt_ERA5',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='slwt_ERA5',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='slwt_ERA5',
                        preproc_path=preproc_path_tas,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='slwt_EOBS',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='slwt_EOBS',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='slwt_EOBS',
                        preproc_path=preproc_path_tas,
                    )
            else:
                if key == 'E-OBS':
                    continue
                wt_preproc = iris.load_cube(value[0].get('filename'))
                mean_preproc_psl = iris.load_cube(value[1].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[2].get('filename'))
                mean_preproc_tas = iris.load_cube(value[3].get('filename'))

                model_name = key
                timerange = value[0].get('timerange').replace('/', '-')
                experiment = value[0].get('exp')
                ensemble = value[0].get('ensemble')

                output_file_path = f'{
                    model_name}/{experiment}/{ensemble}/{timerange}'
                preproc_path = value[0].get('filename')

                # calculate weathertypes
                calculate_lwt_slwt_model(
                    cfg, wt_preproc, key, preproc_path, output_file_path)

                # load wt files
                lwt_cube = iris.load_cube(
                    f'{work_dir}/{output_file_path}/{key}.nc', 'lwt')
                slwt_era5_cube = iris.load_cube(
                    f'{work_dir}/{output_file_path}/{key}.nc', 'slwt_era5')
                slwt_eobs_cube = iris.load_cube(
                    f'{work_dir}/{output_file_path}/{key}.nc', 'slwt_eobs')
                wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]

                preproc_path_psl = value[1].get('filename')
                preproc_path_prcp = value[2].get('filename')
                preproc_path_tas = value[3].get('filename')

                if plotting:
                    # plot means
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='lwt',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='lwt',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='lwt',
                        preproc_path=preproc_path_tas,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='slwt_ERA5',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='slwt_ERA5',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='slwt_ERA5',
                        preproc_path=preproc_path_tas,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='slwt_EOBS',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='slwt_EOBS',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='slwt_EOBS',
                        preproc_path=preproc_path_tas,
                    )
    elif not automatic_slwt:  # if automatic_slwt is false, and predefined_slwt is false, just write selected pairs to file
        for key, value in preproc_variables_dict.items():
            if key == 'ERA5':
                wt_preproc = iris.load_cube(value[0].get('filename'))
                wt_preproc_prcp = iris.load_cube(value[1].get('filename'))
                mean_preproc_psl = iris.load_cube(value[2].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[3].get('filename'))
                mean_preproc_tas = iris.load_cube(value[4].get('filename'))
                wt_preproc_prcp_eobs = iris.load_cube(
                    preproc_variables_dict.get('E-OBS')[0].get('filename')
                )

                # calculate lwt
                lwt = wt_algorithm(wt_preproc, key)

                # calculate simplified lwt based on precipitation patterns
                slwt_era5 = calculate_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp,
                    key,
                    correlation_threshold,
                    rmse_threshold,
                    era5_ancestors,
                )
                slwt_eobs = calculate_slwt_obs(
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

                # load weathertype files as cubes
                lwt_cube = iris.load_cube(f'{work_dir}/{key}.nc', 'lwt')
                wt_cubes = [lwt_cube]

                preproc_path_psl = value[2].get('filename')
                preproc_path_prcp = value[3].get('filename')
                preproc_path_tas = value[4].get('filename')

                if plotting:
                    # plot means
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='lwt',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='lwt',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='lwt',
                        preproc_path=preproc_path_tas,
                    )
            else:
                if key == 'E-OBS':
                    continue
                wt_preproc = iris.load_cube(value[0].get('filename'))
                mean_preproc_psl = iris.load_cube(value[1].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[2].get('filename'))
                mean_preproc_tas = iris.load_cube(value[3].get('filename'))

                model_name = key
                timerange = value[0].get('timerange').replace('/', '-')
                experiment = value[0].get('exp')
                ensemble = value[0].get('ensemble')

                output_file_path = f'{
                    model_name}/{experiment}/{ensemble}/{timerange}'

                # calculate weathertypes
                calculate_lwt_model(cfg, wt_preproc, key, output_file_path)

                # load wt files
                lwt_cube = iris.load_cube(f'{work_dir}/{key}.nc', 'lwt')
                wt_cubes = [lwt_cube]

                # plot means
                preproc_path_psl = value[2].get('filename')
                preproc_path_prcp = value[3].get('filename')
                preproc_path_tas = value[4].get('filename')

                if plotting:
                    # plot means
                    calculate_wt_means(
                        cfg,
                        mean_preproc_psl,
                        wt_cubes,
                        key,
                        var_name='psl',
                        wt_string='lwt',
                        preproc_path=preproc_path_psl,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_prcp,
                        wt_cubes,
                        key,
                        var_name='prcp',
                        wt_string='lwt',
                        preproc_path=preproc_path_prcp,
                    )
                    calculate_wt_means(
                        cfg,
                        mean_preproc_tas,
                        wt_cubes,
                        key,
                        var_name='tas',
                        wt_string='lwt',
                        preproc_path=preproc_path_tas,
                    )


if __name__ == '__main__':

    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
