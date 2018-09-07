"""
Base autoassess area metrics diagnostic

Wrapper that takes two datasets (control_model and exp_model
and observational data (optionally); base for all area diags for
autoassess metrics. Current areas:

version.September.2018
-----------------------
monsoon -- not yet implemented
stratosphere -- implemented
hydrocycle -- not yet implemented
conservation -- not yet implemented
globaltrop -- not yet implemented
land_surface_surfrad -- implemented
land_surface_snow -- implemented
custom -- not yet implemented

Author: Valeriu Predoi, UREAD (valeriu.predoi@ncas.ac.uk)
First version: September 2018
"""
import os
import datetime
import logging
import importlib
import csv
import tempfile
import iris
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(__name__)


def _import_package(area):
    """Import the right area package"""
    root_import = 'esmvaltool.diag_scripts.autoassess.'
    available_areas = [
        'monsoon', 'stratosphere', 'hydrocycle', 'conservation', 'globaltrop',
        'land_surface_surfrad', 'land_surface_snow'
    ]
    if area in available_areas:
        module = root_import + area
        area_package = importlib.import_module(module)
        return area_package
    else:
        raise Exception('Unknown area: ' + area)


def _make_tmp_dir(cfg):
    """Make the tmp and ancil dirs"""
    tmp_dir = os.path.join(cfg['work_dir'], 'tmp')
    ancil_dir = os.path.join(cfg['work_dir'], 'ancil')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not os.path.exists(ancil_dir):
        os.makedirs(ancil_dir)
    return tmp_dir, ancil_dir


def _make_main_dirs(cfg, obs4metrics=False):
    """Create main dirs to hold analysis"""
    locations = []
    suite_loc_m1 = os.path.join(cfg['work_dir'], cfg['control_model'])
    if not os.path.exists(suite_loc_m1):
        os.makedirs(suite_loc_m1)
    locations.append(suite_loc_m1)
    suite_loc_m2 = os.path.join(cfg['work_dir'], cfg['exp_model'])
    if not os.path.exists(suite_loc_m2):
        os.makedirs(suite_loc_m2)
    locations.append(suite_loc_m2)
    if obs4metrics:
        for obs_model in cfg['obs4metrics']:
            suite_loc_obs = os.path.join(cfg['work_dir'], obs_model)
            if not os.path.exists(suite_loc_obs):
                os.makedirs(suite_loc_obs)
        locations.append(suite_loc_obs)
    obs_loc = os.path.join(cfg['work_dir'], 'OBS')
    if not os.path.exists(obs_loc):
        os.makedirs(obs_loc)
    locations.append(obs_loc)

    return locations, obs_loc


def _make_concatenated_data_dirs(suite_locs, area):
    """Create dirs to hold cubeList files"""
    suites_locations = []
    supermeans_locations = []
    for suite_dir in suite_locs:
        suite_data = os.path.join(suite_dir, area)
        if not os.path.exists(suite_data):
            os.makedirs(suite_data)
        suites_locations.append(suite_data)

        # create supermeans directory: [area]_supermeans
        sup_data = os.path.join(suite_dir, area + '_supermeans')
        if not os.path.exists(sup_data):
            os.makedirs(sup_data)
        supermeans_locations.append(sup_data)

    return suites_locations, supermeans_locations


def _get_filelists(cfg, obs4metrics=False, obs_type=None):
    """Put files in lists and return them"""
    files_list_m1 = []
    files_list_m2 = []
    obs4metrics_list = []
    obs_list = []
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        fullpath_file = filename
        if base_file.split('_')[1] == cfg['control_model']:
            files_list_m1.append(fullpath_file)
            if 'fx_files' in attributes:
                files_list_m1.append(attributes['fx_files'][cfg['fx']])
        if base_file.split('_')[1] == cfg['exp_model']:
            files_list_m2.append(fullpath_file)
            if 'fx_files' in attributes:
                files_list_m2.append(attributes['fx_files'][cfg['fx']])
        if obs4metrics and base_file.split('_')[1] in cfg['obs4metrics']:
            obs4metrics_list.append(fullpath_file)
        if obs_type and base_file.split('_')[0] == obs_type:
            obs_list.append(fullpath_file)
    return files_list_m1, files_list_m2, obs4metrics_list, obs_list


def _process_obs(cfg, obs_list, obs_loc):
    """Gather obs files and save them applying specific cases"""
    group_files = [[
        ofile for ofile in obs_list
        if os.path.basename(ofile).split('_')[1] == obs
    ] for obs in cfg['obs_models']]
    for obs_file_group, obs_name in zip(group_files, cfg['obs_models']):
        cubes_list_obs = iris.load(obs_file_group)
        # force add a long_name; supermeans uses extract_strict
        # and for derived vars there is only
        # invalid_standard_name which is an attribute and shit
        for cube in cubes_list_obs:
            if 'invalid_standard_name' in cube.attributes:
                cube.long_name = cube.attributes['invalid_standard_name']
            coord_names = [coord.standard_name for coord in cube.coords()]
            if 'time' in coord_names:
                if not cube.coord('time').has_bounds():
                    cube.coord('time').guess_bounds()
        obs_file_name = obs_name + '_cubeList.nc'
        iris.save(cubes_list_obs, os.path.join(obs_loc, obs_file_name))


def _process_metrics_data(all_files, suites, smeans):
    """Create and save concatenated cubes for ctrl and exp"""
    cubes_lists_paths = []
    for filelist, suite, smean in zip(all_files, suites, smeans):
        if filelist:
            cubelist = iris.load(filelist)

            # save to congragated files; save twice for supermeans as well
            cubes_list_path = os.path.join(suite, 'cubeList.nc')
            cubes_list_smean_path = os.path.join(smean, 'cubeList.nc')
            # force add a long_name; supermeans uses extract_strict
            # and for derived vars there is only
            # invalid_standard_name which is an attribute and shit
            for cube in cubelist:
                if 'invalid_standard_name' in cube.attributes:
                    cube.long_name = cube.attributes['invalid_standard_name']
                coord_names = [coord.standard_name for coord in cube.coords()]
                if 'time' in coord_names:
                    if not cube.coord('time').has_bounds():
                        cube.coord('time').guess_bounds()
            iris.save(cubelist, cubes_list_path)
            iris.save(cubelist, cubes_list_smean_path)
            cubes_lists_paths.append(cubes_list_path)

    return cubes_lists_paths


def create_output_tree(out_dir, ref_suite_id, exp_suite_id, area):
    """
    Create directory tree for area output according to the following scheme:

        `out_dir`/`exp_suite_id`_vs_`ref_suite_id`/`area`

    If the leaf directory `area` exists raises OSError.

    :param str out_dir: Base directory for output.
    :param str suite_id1: Suite Id of reference model run.
    :param str suite_id2: Suite Id of test model run.
    :param str area: Name of asssessment area.
    :returns: Path to area output directory.
    :rtype: str
    :raises: OSError
    """
    assessment_name = exp_suite_id + '_vs_' + ref_suite_id
    # make sure out_dir exists in output folder
    _out_dir = os.path.join(out_dir, assessment_name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # create output folder for area
    area_out_dir = os.path.join(_out_dir, area)
    if not os.path.exists(area_out_dir):
        os.makedirs(area_out_dir)
    return area_out_dir


def create_tmp_dir(tmp_dir, ref_suite_id, exp_suite_id, area):
    """
    Create directory tree for temporary data according to the following scheme:

        `tmp_dir`/`exp_suite_id`_vs_`ref_suite_id`_random/`area`_random

    :param str tmp_dir: Base temporary directory.
    :param str suite_id1: Suite ID of reference model run.
    :param str suite_id2: Suite ID of test model run.
    :param str area: Name of asssessment area.
    :returns: Path to area temporary directory.
    :rtype: str
    """
    assessment_name = exp_suite_id + '_vs_' + ref_suite_id
    # create unique temporary folder in tmp dir
    _tmp_dir = tempfile.mkdtemp(prefix=assessment_name + '_', dir=tmp_dir)

    # create temporary folder for area
    area_tmp_dir = tempfile.mkdtemp(prefix=area + '_', dir=_tmp_dir)
    return area_tmp_dir


def _setup_input(cfg):
    """Assemble all data structures"""
    logger.setLevel(cfg['log_level'].upper())

    # set the main data dirs;
    # check if we need to apply metrics to obs's
    if 'obs4metrics' not in cfg:
        target_locs, obs_loc = _make_main_dirs(cfg)
    else:
        if not cfg['obs4metrics']:
            target_locs, obs_loc = _make_main_dirs(cfg)
        else:
            target_locs, obs_loc = _make_main_dirs(cfg, obs4metrics=True)
    suites, smeans = _make_concatenated_data_dirs(target_locs, cfg['area'])

    # create the ancil and tp dirs
    tmp_dir, ancil_dir = _make_tmp_dir(cfg)

    # get files lists
    if 'obs_type' and 'obs4metrics' in cfg:
        files_list_m1, files_list_m2, obs4metrics_list, obs_list = \
            _get_filelists(cfg, obs4metrics=True, obs_type=cfg['obs_type'])
    else:
        if 'obs_type' in cfg:
            files_list_m1, files_list_m2, obs4metrics_list, obs_list = \
                _get_filelists(cfg, obs4metrics=False,
                               obs_type=cfg['obs_type'])
        elif 'obs4metrics' in cfg:
            files_list_m1, files_list_m2, obs4metrics_list, obs_list = \
                _get_filelists(cfg, obs4metrics=True)
        else:
            files_list_m1, files_list_m2, obs4metrics_list, obs_list = \
                _get_filelists(cfg)

    # spell out the files used
    logger.info("Files for control model for metrics: %s", files_list_m1)
    logger.info("Files for exp model for metrics: %s", files_list_m2)
    logger.info("Files for obs models for metrics: %s", obs4metrics_list)
    logger.info("Files for obs model NOT for metrics: %s", obs_list)
    all_files = [files_list_m1, files_list_m2, obs4metrics_list]

    # load and save control and exp cubelists
    all_cubelists = _process_metrics_data(all_files, suites, smeans)

    # print the paths
    logger.info("Saved control data cubes: %s", str(all_cubelists))

    # separately process the obs's that dont need metrics
    if cfg['obs_models']:
        _process_obs(cfg, obs_list, obs_loc)

    return tmp_dir, obs_loc, ancil_dir


def _create_run_dict(cfg):
    """Create the run dictionary"""
    tmp_dir, obs_loc, ancil_dir = _setup_input(cfg)
    run = {}
    # general parameters (necessary)
    run['suite_id1'] = cfg['control_model']
    run['suite_id2'] = cfg['exp_model']
    run['out_dir'] = cfg['plot_dir']
    run['tmp_dir'] = tmp_dir
    run['_area'] = cfg['area']
    run['_start_date'] = cfg['start']
    run['_end_date'] = cfg['end']
    run['runid'] = cfg['area']
    run['data_root'] = cfg['work_dir']
    run['clim_root'] = obs_loc
    run['ancil_root'] = ancil_dir
    run['start'] = cfg['start']
    run['end'] = cfg['end']

    # optional parameters
    if 'climfiles_root' in cfg:
        run['climfiles_root'] = cfg['climfiles_root']
    if 'obs4metrics' in cfg:
        run['obs4metrics'] = cfg['obs4metrics']

    # specific parameters needed by some areas
    start_year = int(run['start'][0:4])
    end_year = int(run['end'][0:4])
    run['nyear'] = end_year - start_year
    run['period'] = '{:04d}_{:03d}'.format(start_year, run['nyear'])
    year, month, day = [int(s) for s in run['start'].split('/')]
    run['from_instantaneous'] = datetime.datetime(year, month, day)
    run['from_daily'] = datetime.datetime(year, month, day)
    run['from_monthly'] = datetime.datetime(year, month, day)
    run['from_seasonal'] = datetime.datetime(year, month, day)
    run['from_annual'] = datetime.datetime(year, month, day)

    year, month, day = [int(s) for s in run['end'].split('/')]
    assert month == 12 and day == 1  # Climatological year
    run['to_instantaneous'] = datetime.datetime(year, 11, 30)
    run['to_daily'] = datetime.datetime(year, 11, 30)
    run['to_monthly'] = datetime.datetime(year, 11, 1)
    run['to_seasonal'] = datetime.datetime(year, 9, 1)
    run['to_annual'] = datetime.datetime(year - 1, 12, 1)

    return run


def run_area(cfg):
    """Kick start the area diagnostic"""
    run_obj = _create_run_dict(cfg)
    area_out_dir = create_output_tree(run_obj['out_dir'], run_obj['suite_id1'],
                                      run_obj['suite_id2'], run_obj['_area'])

    # the areas write all output to the cwd
    os.chdir(area_out_dir)

    # import area here to allow removal of areas
    area_package = _import_package(run_obj['_area'])

    # assemble the work subjects
    suite_ids = [run_obj['suite_id1'], run_obj['suite_id2']]
    if 'obs4metrics' in cfg:
        if run_obj['obs4metrics']:
            suite_ids.extend(run_obj['obs4metrics'])

    # run the metrics generation
    for suite_id in suite_ids:
        logger.info('Calculating metrics for %s', suite_id)
        all_metrics = {}

        # run each metric function
        for metric_function in area_package.metrics_functions:
            logger.info('# Call: %s', metric_function)
            run_obj['runid'] = suite_id
            metrics = metric_function(run_obj)
            duplicate_metrics = list(
                set(all_metrics.keys()) & set(metrics.keys()))
            if duplicate_metrics:
                raise AssertionError('Duplicate Metrics ' +
                                     str(duplicate_metrics))
            all_metrics.update(metrics)

        # write metrics to file
        csv_dir = os.path.join(area_out_dir, suite_id)
        if not os.path.exists(csv_dir):
            os.makedirs(csv_dir)
        with open(os.path.join(area_out_dir, suite_id, 'metrics.csv'),
                  'w') as file_handle:
            writer = csv.writer(file_handle)
            for metric in all_metrics.items():
                writer.writerow(metric)

    # multimodel functions
    if hasattr(area_package, 'multi_functions'):
        for multi_function in area_package.multi_functions:
            multi_function(run_obj)
    else:
        logger.info('# Area has no multi functions.')


if __name__ == '__main__':

    with run_diagnostic() as config:
        run_area(config)
