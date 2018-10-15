"""Recipe parser"""
import copy
import fnmatch
import inspect
import logging
import os
import subprocess
from collections import OrderedDict

import yamale
import yaml

from . import __version__, preprocessor
from ._config import get_institutes
from ._data_finder import (get_input_filelist, get_input_fx_filelist,
                           get_output_file, get_rootpath, get_start_end_year,
                           get_statistic_output_file)
from ._task import DiagnosticTask, get_independent_tasks, run_tasks, which
from .cmor.table import CMOR_TABLES
from .preprocessor import DEFAULT_ORDER, FINAL_STEPS, INITIAL_STEPS
from .preprocessor._derive import get_required
from .preprocessor._download import synda_search
from .preprocessor._io import DATASET_KEYS, concatenate_callback
from .preprocessor._regrid import get_cmor_levels, get_reference_levels

logger = logging.getLogger(__name__)

TASKSEP = os.sep


class RecipeError(Exception):
    """Recipe contains an error."""


def ordered_safe_load(stream):
    """Load a YAML file using OrderedDict instead of dict"""

    class OrderedSafeLoader(yaml.SafeLoader):
        """Loader class that uses OrderedDict to load a map"""

    def construct_mapping(loader, node):
        """Load a map as an OrderedDict"""
        loader.flatten_mapping(node)
        return OrderedDict(loader.construct_pairs(node))

    OrderedSafeLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)

    return yaml.load(stream, OrderedSafeLoader)


def read_recipe_file(filename, config_user, initialize_tasks=True):
    """Read a recipe from file."""
    raw_recipe = check_recipe(filename)
    return Recipe(
        raw_recipe, config_user, initialize_tasks, recipe_file=filename)


def check_ncl_version():
    """Check the NCL version"""
    ncl = which('ncl')
    if not ncl:
        raise RecipeError("Recipe contains NCL scripts, but cannot find "
                          "an NCL installation.")
    try:
        cmd = [ncl, '-V']
        version = subprocess.check_output(cmd, universal_newlines=True)
    except subprocess.CalledProcessError:
        logger.error("Failed to execute '%s'", ' '.join(' '.join(cmd)))
        raise RecipeError("Recipe contains NCL scripts, but your NCL "
                          "installation appears to be broken.")

    version = version.strip()
    logger.info("Found NCL version %s", version)

    major, minor = (int(i) for i in version.split('.')[:2])
    if major < 6 or (major == 6 and minor < 4):
        raise RecipeError("NCL version 6.4 or higher is required to run "
                          "a recipe containing NCL scripts.")


def check_recipe_with_schema(filename):
    """Check if the recipe content matches schema."""
    schema_file = os.path.join(os.path.dirname(__file__), 'recipe_schema.yml')
    logger.debug("Checking recipe against schema %s", schema_file)
    recipe = yamale.make_data(filename)
    schema = yamale.make_schema(schema_file)
    yamale.validate(schema, recipe)


def check_recipe(filename):
    """Check a recipe file and return it in raw form."""
    # Note that many checks can only be performed after the automatically
    # computed entries have been filled in by creating a Recipe object.
    check_recipe_with_schema(filename)
    with open(filename, 'r') as file:
        contents = file.read()
        raw_recipe = yaml.safe_load(contents)
        raw_recipe['preprocessors'] = ordered_safe_load(contents).get(
            'preprocessors', {})

    check_preprocessors(raw_recipe['preprocessors'])
    check_diagnostics(raw_recipe['diagnostics'])
    return raw_recipe


def check_preprocessors(preprocessors):
    """Check preprocessors in recipe"""
    valid_functions = set(preprocessor.DEFAULT_ORDER)
    for name, profile in preprocessors.items():
        invalid_functions = set(profile) - {'custom_order'} - valid_functions
        if invalid_functions:
            raise RecipeError(
                "Unknown function(s) {} in preprocessor {}, choose from: "
                "{}".format(', '.join(invalid_functions), name,
                            ', '.join(preprocessor.DEFAULT_ORDER)))


def check_diagnostics(diagnostics):
    """Check diagnostics in recipe"""
    for name, diagnostic in diagnostics.items():
        if 'scripts' not in diagnostic:
            raise RecipeError("Missing scripts section in diagnostic {}"
                              .format(name))
        variable_names = tuple(diagnostic.get('variables', {}))
        scripts = diagnostic.get('scripts')
        if scripts is None:
            scripts = {}
        for script_name, script in scripts.items():
            if script_name in variable_names:
                raise RecipeError(
                    "Invalid script name {} encountered in diagnostic {}: "
                    "scripts cannot have the same name as variables.".format(
                        script_name, name))
            if not script.get('script'):
                raise RecipeError(
                    "No script defined for script {} in diagnostic {}".format(
                        script_name, name))


def check_preprocessor_settings(settings):
    """Check preprocessor settings."""
    # The inspect functions getargspec and getcallargs are deprecated
    # in Python 3, but their replacements are not available in Python 2.
    # TODO: Use the new Python 3 inspect API
    for step in settings:
        if step not in preprocessor.DEFAULT_ORDER:
            raise RecipeError(
                "Unknown preprocessor function '{}', choose from: {}".format(
                    step, ', '.join(preprocessor.DEFAULT_ORDER)))

        function = getattr(preprocessor, step)
        argspec = inspect.getargspec(function)
        args = argspec.args[1:]
        # Check for invalid arguments
        invalid_args = set(settings[step]) - set(args)
        if invalid_args:
            raise RecipeError(
                "Invalid argument(s): {} encountered for preprocessor "
                "function {}. \nValid arguments are: [{}]".format(
                    ', '.join(invalid_args), step, ', '.join(args)))

        # Check for missing arguments
        defaults = argspec.defaults
        end = None if defaults is None else -len(defaults)
        missing_args = set(args[:end]) - set(settings[step])
        if missing_args:
            raise RecipeError(
                "Missing required argument(s) {} for preprocessor "
                "function {}".format(missing_args, step))
        # Final sanity check in case the above fails to catch a mistake
        try:
            inspect.getcallargs(function, None, **settings[step])
        except TypeError:
            logger.error(
                "Wrong preprocessor function arguments in "
                "function '%s'", step)
            raise


def check_duplicate_datasets(datasets):
    """Check for duplicate datasets."""
    checked_datasets_ = []
    for dataset in datasets:
        if dataset in checked_datasets_:
            raise RecipeError(
                "Duplicate dataset {} in datasets section".format(dataset))
        checked_datasets_.append(dataset)


def check_variable(variable, required_keys):
    """Check variables as derived from recipe"""
    required = set(required_keys)
    missing = required - set(variable)
    if missing:
        raise RecipeError(
            "Missing keys {} from variable {} in diagnostic {}".format(
                missing, variable.get('short_name'),
                variable.get('diagnostic')))


def check_data_availability(input_files, variable):
    """Check if the required input data is available"""
    if not input_files:
        raise RecipeError("No input files found for variable {}"
                          .format(variable))

    required_years = set(
        range(variable['start_year'], variable['end_year'] + 1))
    available_years = set()
    for filename in input_files:
        start, end = get_start_end_year(filename)
        available_years.update(range(start, end + 1))

    missing_years = required_years - available_years
    if missing_years:
        raise RecipeError(
            "No input data available for years {} in files {}".format(
                ", ".join(str(year) for year in missing_years), input_files))


def _get_value(key, datasets):
    """Get a value for key by looking at the other datasets."""
    values = {dataset[key] for dataset in datasets if key in dataset}

    if len(values) > 1:
        raise RecipeError("Ambigous values {} for property {}".format(
            values, key))

    value = None
    if len(values) == 1:
        value = values.pop()

    return value


def _update_from_others(variable, keys, datasets):
    """Get values for keys by copying from the other datasets."""
    for key in keys:
        if key not in variable:
            value = _get_value(key, datasets)
            if value is not None:
                variable[key] = value


def _update_cmor_table(table, mip, short_name):
    """Try to add an ESMValTool custom CMOR table file."""
    cmor_table = CMOR_TABLES[table]
    var_info = cmor_table.get_variable(mip, short_name)

    if var_info is None and hasattr(cmor_table, 'add_custom_table_file'):
        table_file = os.path.join(
            os.path.dirname(__file__), 'cmor', 'tables', 'custom',
            'CMOR_' + short_name + '.dat')
        if os.path.exists(table_file):
            logger.debug("Loading custom CMOR table from %s", table_file)
            cmor_table.add_custom_table_file(table_file, mip)
            var_info = cmor_table.get_variable(mip, short_name)

    if var_info is None:
        raise RecipeError(
            "Unable to load CMOR table '{}' for variable '{}' with mip '{}'"
            .format(table, short_name, mip))


def _add_cmor_info(variable, override=False):
    """Add information from CMOR tables to variable."""
    logger.debug("If not present: adding keys from CMOR table to %s", variable)

    if 'cmor_table' not in variable or 'mip' not in variable:
        logger.debug("Skipping because cmor_table or mip not specified")
        return

    if variable['cmor_table'] not in CMOR_TABLES:
        logger.warning("Unknown CMOR table %s", variable['cmor_table'])

    # Copy the following keys from CMOR table
    cmor_keys = [
        'standard_name', 'long_name', 'units', 'modeling_realm', 'frequency'
    ]
    table_entry = CMOR_TABLES[variable['cmor_table']].get_variable(
        variable['mip'], variable['short_name'])

    for key in cmor_keys:
        if key not in variable or override:
            value = getattr(table_entry, key, None)
            if value is not None:
                variable[key] = value
            else:
                logger.debug(
                    "Failed to add key %s to variable %s from CMOR table", key,
                    variable)

    # Check that keys are available
    check_variable(variable, required_keys=cmor_keys)


def _special_name_to_dataset(variable, special_name):
    """Convert special names to dataset names."""
    if special_name in ('reference_dataset', 'alternative_dataset'):
        if special_name not in variable:
            raise RecipeError(
                "Preprocessor {} uses {}, but {} is not defined for "
                "variable {} of diagnostic {}".format(
                    variable['preprocessor'], special_name, special_name,
                    variable['short_name'], variable['diagnostic']))
        special_name = variable[special_name]

    return special_name


def _update_target_levels(variable, variables, settings, config_user):
    """Replace the target levels dataset name with a filename if needed."""
    levels = settings.get('extract_levels', {}).get('levels')
    if not levels:
        return

    levels = _special_name_to_dataset(variable, levels)

    # If levels is a dataset name, replace it by a dict with a 'dataset' entry
    if any(levels == v['dataset'] for v in variables):
        settings['extract_levels']['levels'] = {'dataset': levels}
        levels = settings['extract_levels']['levels']

    if not isinstance(levels, dict):
        return

    if 'cmor_table' in levels and 'coordinate' in levels:
        settings['extract_levels']['levels'] = get_cmor_levels(
            levels['cmor_table'], levels['coordinate'])
    elif 'dataset' in levels:
        dataset = levels['dataset']
        if variable['dataset'] == dataset:
            del settings['extract_levels']
        else:
            variable_data = _get_dataset_info(dataset, variables)
            filename = \
                _dataset_to_file(variable_data, config_user)
            coordinate = levels.get('coordinate', 'air_pressure')
            settings['extract_levels']['levels'] = get_reference_levels(
                filename,
                variable_data['project'], dataset, variable_data['short_name'],
                os.path.splitext(variable_data['filename'])[0] + '_fixed',
                coordinate)


def _update_target_grid(variable, variables, settings, config_user):
    """Replace the target grid dataset name with a filename if needed."""
    grid = settings.get('regrid', {}).get('target_grid')
    if not grid:
        return

    grid = _special_name_to_dataset(variable, grid)

    if variable['dataset'] == grid:
        del settings['regrid']
    elif any(grid == v['dataset'] for v in variables):
        settings['regrid']['target_grid'] = _dataset_to_file(
            _get_dataset_info(grid, variables), config_user)


def _get_dataset_info(dataset, variables):
    for var in variables:
        if var['dataset'] == dataset:
            return var
    raise RecipeError(
        "Unable to find matching file for dataset"
        "{}".format(dataset)
    )


def _dataset_to_file(variable, config_user):
    """Find the first file belonging to dataset from variable info."""
    files = get_input_filelist(
        variable=variable,
        rootpath=config_user['rootpath'],
        drs=config_user['drs'])
    if not files and variable.get('derive'):
        variable = copy.deepcopy(variable)
        variable['short_name'], variable['field'] = get_required(
            variable['short_name'], variable['field'])[0]
        files = get_input_filelist(
            variable=variable,
            rootpath=config_user['rootpath'],
            drs=config_user['drs'])
    check_data_availability(files, variable)
    return files[0]


def _limit_datasets(variables, profile, max_datasets=None):
    """Try to limit the number of datasets to max_datasets."""
    if not max_datasets:
        return variables

    logger.info("Limiting the number of datasets to %s", max_datasets)

    required_datasets = (
        profile.get('extract_levels', {}).get('levels'),
        profile.get('regrid', {}).get('target_grid'),
        variables[0].get('reference_dataset'),
        variables[0].get('alternative_dataset'),
    )

    limited = []

    for variable in variables:
        if variable['dataset'] in required_datasets:
            limited.append(variable)

    for variable in variables[::-1]:
        if len(limited) >= max_datasets:
            break
        if variable not in limited:
            limited.append(variable)

    logger.info("Only considering %s",
                ', '.join(v['dataset'] for v in limited))

    return limited


def _get_default_settings(variable, config_user, derive=False):
    """Get default preprocessor settings."""
    settings = {}

    # Set up downloading using synda if requested.
    if config_user['synda_download']:
        # TODO: make this respect drs or download to preproc dir?
        local_dir = get_rootpath(config_user['rootpath'], variable['project'])
        settings['download'] = {
            'dest_folder': local_dir,
        }

    # Configure loading
    settings['load_cubes'] = {
        'callback': concatenate_callback,
        'filename': variable['filename'],
        'metadata': variable,
    }
    if not derive:
        settings['load_cubes']['constraints'] = variable['standard_name']
    # Configure merge
    settings['concatenate'] = {}

    # Configure fixes
    fix = {
        'project': variable['project'],
        'dataset': variable['dataset'],
        'short_name': variable['short_name'],
    }
    # File fixes
    fix_dir = os.path.splitext(variable['filename'])[0] + '_fixed'
    if not derive:
        settings['fix_file'] = dict(fix)
        settings['fix_file']['output_dir'] = fix_dir
    # Cube fixes
    # Only supply mip if the CMOR check fixes are implemented.
    if variable.get('cmor_table'):
        fix['cmor_table'] = variable['cmor_table']
        fix['mip'] = variable['mip']
    settings['fix_data'] = dict(fix)
    settings['fix_metadata'] = dict(fix)

    # Configure time extraction
    settings['extract_time'] = {
        'start_year': variable['start_year'],
        'end_year': variable['end_year'] + 1,
        'start_month': 1,
        'end_month': 1,
        'start_day': 1,
        'end_day': 1,
    }

    if derive:
        settings['derive'] = {'variable': variable}

    # Configure CMOR metadata check
    if variable.get('cmor_table'):
        settings['cmor_check_metadata'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
        }
    # Configure final CMOR data check
    if variable.get('cmor_table'):
        settings['cmor_check_data'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
        }

    # Clean up fixed files
    if not config_user['save_intermediary_cubes']:
        settings['cleanup'] = {
            'remove': [fix_dir],
        }

    # Configure saving cubes to file
    settings['save'] = {'compress': config_user['compress_netcdf']}

    return settings


def _update_fx_settings(settings, variable, config_user):
    """Find and set the FX mask settings"""
    # update for landsea
    if 'mask_landsea' in settings.keys():
        # Configure ingestion of land/sea masks
        logger.debug('Getting fx mask settings now...')

        # settings[mask_landsea][fx_file] is a list to store ALL
        # available masks
        settings['mask_landsea']['fx_files'] = []

        # fx_files already in variable
        variable = dict(variable)
        variable['fx_files'] = ['sftlf', 'sftof']
        fx_files_dict = get_input_fx_filelist(
            variable=variable,
            rootpath=config_user['rootpath'],
            drs=config_user['drs'])

        # allow both sftlf and sftof
        if fx_files_dict['sftlf']:
            settings['mask_landsea']['fx_files'].append(fx_files_dict['sftlf'])
        if fx_files_dict['sftof']:
            settings['mask_landsea']['fx_files'].append(fx_files_dict['sftof'])
    # update for landseaice
    if 'mask_landseaice' in settings.keys():
        # Configure ingestion of land/sea masks
        logger.debug('Getting fx mask settings now...')

        # settings[mask_landseaice][fx_file] is a list to store ALL
        # available masks
        settings['mask_landseaice']['fx_files'] = []

        # fx_files already in variable
        variable = dict(variable)
        variable['fx_files'] = ['sftgif']
        fx_files_dict = get_input_fx_filelist(
            variable=variable,
            rootpath=config_user['rootpath'],
            drs=config_user['drs'])

        # allow sftgif (only, for now)
        if fx_files_dict['sftgif']:
            settings['mask_landseaice']['fx_files'].append(
                fx_files_dict['sftgif'])


def _get_input_files(variable, config_user):
    """Get the input files for a single dataset"""
    # Find input files locally.
    input_files = get_input_filelist(
        variable=variable,
        rootpath=config_user['rootpath'],
        drs=config_user['drs'])

    # Set up downloading using synda if requested.
    # Do not download if files are already available locally.
    if config_user['synda_download'] and not input_files:
        input_files = synda_search(variable)

    logger.info("Using input files for variable %s of dataset %s:\n%s",
                variable['short_name'], variable['dataset'],
                '\n'.join(input_files))
    check_data_availability(input_files, variable)

    return input_files


def _apply_preprocessor_settings(settings, profile_settings):
    """Apply settings from preprocessor profile."""
    for step, args in profile_settings.items():
        # Remove disabled preprocessor functions
        if args is False:
            if step in settings:
                del settings[step]
            continue
        # Enable/update functions without keywords
        if step not in settings:
            settings[step] = {}
        if isinstance(args, dict):
            settings[step].update(args)


def _update_multi_model_statistics(variables, settings, preproc_dir):
    """Configure multi dataset statistics."""
    if settings.get('multi_model_statistics', False):
        if settings['multi_model_statistics'] is True:
            settings['multi_model_statistics'] = {}
        stat_settings = settings['multi_model_statistics']

        variable = variables[0]

        # Define output files
        stat_settings['filenames'] = {}
        for statistic in stat_settings['statistics']:
            stat_settings['filenames'][statistic] = get_statistic_output_file(
                variable, statistic, preproc_dir)

        # Define datasets to exclude
        exclude_datasets = set(stat_settings.get('exclude', {}))
        for key in 'reference_dataset', 'alternative_dataset':
            if key in exclude_datasets and key in variable:
                exclude_datasets.remove(key)
                exclude_datasets.add(variable[key])
        exclude_files = {
            v['filename']
            for v in variables if v['dataset'] in exclude_datasets
        }
        logger.debug('Multidataset excludes files %s', exclude_files)
        stat_settings['exclude'] = {'_filename': exclude_files}


def _get_preprocessor_settings(variables, profile, config_user):
    """Get preprocessor settings for a set of datasets."""
    all_settings = {}
    profile = copy.deepcopy(profile)
    _update_multi_model_statistics(variables, profile,
                                   config_user['preproc_dir'])

    for variable in variables:
        derive = 'derive' in profile
        settings = _get_default_settings(variable, config_user, derive=derive)
        _apply_preprocessor_settings(settings, profile)
        # if the target grid is a dataset name, replace it with a file name
        # TODO: call _update_target_grid only once per variable?
        _update_target_levels(
            variable=variable,
            variables=variables,
            settings=settings,
            config_user=config_user)
        _update_fx_settings(
            settings=settings, variable=variable, config_user=config_user)
        _update_target_grid(
            variable=variable,
            variables=variables,
            settings=settings,
            config_user=config_user)
        check_preprocessor_settings(settings)
        all_settings[variable['filename']] = settings

    _check_multi_model_settings(all_settings)
    return all_settings


def _check_multi_model_settings(all_settings):
    """Check that multi dataset settings are identical for all datasets."""
    multi_model_steps = (step for step in preprocessor.MULTI_MODEL_FUNCTIONS
                         if any(step in settings
                                for settings in all_settings.values()))
    for step in multi_model_steps:
        result = None
        for filename, settings in all_settings.items():
            if result is None:
                result = settings[step]
            elif result != settings[step]:
                raise RecipeError(
                    "Unable to combine differing multi-dataset settings "
                    "{} and {} for output file {}".format(
                        result, settings[step], filename))


def _extract_preprocessor_order(profile):
    """Extract the order of the preprocessing steps from the profile."""
    custom_order = profile.pop('custom_order', False)
    if not custom_order:
        return DEFAULT_ORDER
    order = tuple(p for p in profile if p not in INITIAL_STEPS + FINAL_STEPS)
    return INITIAL_STEPS + order + FINAL_STEPS


def _split_derive_profile(profile):
    """Split the derive preprocessor profile"""
    order = _extract_preprocessor_order(profile)
    before, after = preprocessor.split_settings(profile, 'derive', order)
    after['derive'] = {}
    if order != DEFAULT_ORDER:
        before['custom_order'] = True
        after['custom_order'] = True
    return before, after


def _get_single_preprocessor_task(variables,
                                  profile,
                                  config_user,
                                  ancestors=None):
    """Create preprocessor tasks for a set of datasets."""
    # Configure preprocessor
    order = _extract_preprocessor_order(profile)
    all_settings = _get_preprocessor_settings(
        variables=variables, profile=profile, config_user=config_user)

    # Input files, used by tasks without ancestors
    input_files = None
    if not ancestors:
        input_files = [
            filename for variable in variables
            for filename in _get_input_files(variable, config_user)
        ]

    output_dir = os.path.dirname(variables[0]['filename'])

    task = preprocessor.PreprocessingTask(
        settings=all_settings,
        output_dir=output_dir,
        ancestors=ancestors,
        input_files=input_files,
        order=order,
        debug=config_user['save_intermediary_cubes'])

    return task


def _get_preprocessor_task(variables,
                           profiles,
                           config_user,
                           write_ncl_interface=False):
    """Create preprocessor task(s) for a set of datasets."""
    # First set up the preprocessor profile
    variable = variables[0]
    preproc_name = variable.get('preprocessor')
    if preproc_name not in profiles:
        raise RecipeError(
            "Unknown preprocessor {} in variable {} of diagnostic {}".format(
                preproc_name, variable['short_name'], variable['diagnostic']))
    profile = copy.deepcopy(profiles[variable['preprocessor']])
    logger.info("Creating preprocessor '%s' task for variable '%s'",
                variable['preprocessor'], variable['short_name'])
    variables = _limit_datasets(variables, profile,
                                config_user.get('max_datasets'))

    # Create preprocessor task(s)
    derive_tasks = []
    if variable.get('derive'):
        # Create tasks to prepare the input data for the derive step
        derive_profile, profile = _split_derive_profile(profile)

        derive_input = {}
        for variable in variables:
            _update_cmor_table(
                table=variable['cmor_table'],
                mip=variable['mip'],
                short_name=variable['short_name'])
            _add_cmor_info(variable)
            if not variable.get('force_derivation') and get_input_filelist(
                    variable=variable,
                    rootpath=config_user['rootpath'],
                    drs=config_user['drs']):
                # No need to derive, just process normally up to derive step
                short_name = variable['short_name']
                if short_name not in derive_input:
                    derive_input[short_name] = []
                derive_input[short_name].append(variable)
            else:
                # Process input data needed to derive variable
                for short_name, field in get_required(variable['short_name'],
                                                      variable['field']):
                    if short_name not in derive_input:
                        derive_input[short_name] = []
                    variable = copy.deepcopy(variable)
                    variable['short_name'] = short_name
                    variable['field'] = field
                    variable['filename'] = get_output_file(
                        variable, config_user['preproc_dir'])
                    _add_cmor_info(variable, override=True)
                    derive_input[short_name].append(variable)

        for derive_variables in derive_input.values():
            task = _get_single_preprocessor_task(derive_variables,
                                                 derive_profile, config_user)
            derive_tasks.append(task)

    # Add CMOR info
    for variable in variables:
        _add_cmor_info(variable)

    # Create (final) preprocessor task
    profile['extract_metadata'] = {'write_ncl': write_ncl_interface}
    task = _get_single_preprocessor_task(
        variables, profile, config_user, ancestors=derive_tasks)

    return task


class Recipe(object):
    """Recipe object"""

    def __init__(self,
                 raw_recipe,
                 config_user,
                 initialize_tasks=True,
                 recipe_file=None):
        """Parse a recipe file into an object."""
        self._cfg = config_user
        self._recipe_file = os.path.basename(recipe_file)
        self._preprocessors = raw_recipe.get('preprocessors', {})
        if 'default' not in self._preprocessors:
            self._preprocessors['default'] = {}
        self._support_ncl = self._need_ncl(raw_recipe['diagnostics'])
        self.diagnostics = self._initialize_diagnostics(
            raw_recipe['diagnostics'], raw_recipe.get('datasets', []))
        self.tasks = self.initialize_tasks() if initialize_tasks else None

    @staticmethod
    def _need_ncl(raw_diagnostics):
        if not raw_diagnostics:
            return False
        for diagnostic in raw_diagnostics.values():
            if not diagnostic.get('scripts'):
                continue
            for script in diagnostic['scripts'].values():
                if script.get('script', '').lower().endswith('.ncl'):
                    logger.info("NCL script detected, checking NCL version")
                    check_ncl_version()
                    return True
        return False

    def _initialize_diagnostics(self, raw_diagnostics, raw_datasets):
        """Define diagnostics in recipe"""
        logger.debug("Retrieving diagnostics from recipe")

        diagnostics = {}

        for name, raw_diagnostic in raw_diagnostics.items():
            diagnostic = {}
            diagnostic['name'] = name
            diagnostic['preprocessor_output'] = \
                self._initialize_preprocessor_output(
                    name,
                    raw_diagnostic.get('variables', {}),
                    raw_datasets +
                    raw_diagnostic.get('additional_datasets', []))
            variable_names = tuple(raw_diagnostic.get('variables', {}))
            diagnostic['scripts'] = self._initialize_scripts(
                name, raw_diagnostic.get('scripts'), variable_names)
            diagnostics[name] = diagnostic

        return diagnostics

    @staticmethod
    def _initialize_datasets(raw_datasets):
        """Define datasets used by variable"""
        datasets = copy.deepcopy(raw_datasets)

        for dataset in datasets:
            for key in dataset:
                DATASET_KEYS.add(key)

        check_duplicate_datasets(datasets)
        return datasets

    def _initialize_variables(self, raw_variable, raw_datasets):
        """Define variables for all datasets."""
        variables = []

        datasets = self._initialize_datasets(
            raw_datasets + raw_variable.pop('additional_datasets', []))

        for dataset in datasets:
            variable = dict(raw_variable)
            variable.update(dataset)
            if ('cmor_table' not in variable
                    and variable.get('project') in CMOR_TABLES):
                variable['cmor_table'] = variable['project']
            if 'end_year' in variable and 'max_years' in self._cfg:
                variable['end_year'] = min(
                    variable['end_year'],
                    variable['start_year'] + self._cfg['max_years'] - 1)
            variables.append(variable)

        required_keys = {
            'short_name', 'field', 'dataset', 'project', 'start_year',
            'end_year', 'preprocessor', 'diagnostic'
        }

        for variable in variables:
            _update_from_others(variable, ['cmor_table', 'mip'], datasets)
            institute = get_institutes(variable['dataset'])
            if institute:
                variable['institute'] = institute
            check_variable(variable, required_keys)
            variable['filename'] = get_output_file(variable,
                                                   self._cfg['preproc_dir'])
            if 'fx_files' in variable:
                for fx_file in variable['fx_files']:
                    DATASET_KEYS.add(fx_file)
                # Get the fx files
                variable['fx_files'] = get_input_fx_filelist(
                    variable=variable,
                    rootpath=self._cfg['rootpath'],
                    drs=self._cfg['drs'])
                logger.info("Using fx files for var %s of dataset %s:\n%s",
                            variable['short_name'], variable['dataset'],
                            variable['fx_files'])

        return variables

    def _initialize_preprocessor_output(self, diagnostic_name, raw_variables,
                                        raw_datasets):
        """Define variables in diagnostic"""
        logger.debug("Populating list of variables for diagnostic %s",
                     diagnostic_name)

        preprocessor_output = {}

        for variable_name, raw_variable in raw_variables.items():
            if 'short_name' not in raw_variable:
                raw_variable['short_name'] = variable_name
            raw_variable['diagnostic'] = diagnostic_name
            raw_variable['preprocessor'] = str(
                raw_variable.get('preprocessor', 'default'))
            preprocessor_output[variable_name] = \
                self._initialize_variables(raw_variable, raw_datasets)

        return preprocessor_output

    def _initialize_scripts(self, diagnostic_name, raw_scripts,
                            variable_names):
        """Define script in diagnostic"""
        if not raw_scripts:
            return {}

        logger.debug("Setting script for diagnostic %s", diagnostic_name)

        scripts = {}

        for script_name, raw_settings in raw_scripts.items():
            raw_script = raw_settings.pop('script')
            ancestors = []
            for id_glob in raw_settings.pop('ancestors', variable_names):
                if TASKSEP not in id_glob:
                    id_glob = diagnostic_name + TASKSEP + id_glob
                ancestors.append(id_glob)
            settings = dict(copy.deepcopy(raw_settings))
            settings['recipe'] = self._recipe_file
            settings['version'] = __version__
            settings['script'] = script_name
            # Add output dirs to settings
            for dir_name in ('run_dir', 'plot_dir', 'work_dir'):
                settings[dir_name] = os.path.join(self._cfg[dir_name],
                                                  diagnostic_name, script_name)
            # Copy other settings
            if self._support_ncl:
                settings['exit_on_ncl_warning'] = self._cfg['exit_on_warning']
            for key in ('max_data_filesize', 'output_file_type', 'log_level',
                        'write_plots', 'write_netcdf', 'profile_diagnostic'):
                settings[key] = self._cfg[key]

            scripts[script_name] = {
                'script': raw_script,
                'output_dir': settings['work_dir'],
                'settings': settings,
                'ancestors': ancestors,
            }

        return scripts

    def _resolve_diagnostic_ancestors(self, tasks):
        """Resolve diagnostic ancestors"""
        for diagnostic_name, diagnostic in self.diagnostics.items():
            for script_name, script_cfg in diagnostic['scripts'].items():
                task_id = diagnostic_name + TASKSEP + script_name
                if isinstance(tasks[task_id], DiagnosticTask):
                    logger.debug("Linking tasks for diagnostic %s script %s",
                                 diagnostic_name, script_name)
                    ancestors = []
                    for id_glob in script_cfg['ancestors']:
                        ancestor_ids = fnmatch.filter(tasks, id_glob)
                        if not ancestor_ids:
                            raise RecipeError(
                                "Could not find any ancestors matching {}"
                                .format(id_glob))
                        logger.debug("Pattern %s matches %s", id_glob,
                                     ancestor_ids)
                        ancestors.extend(tasks[a] for a in ancestor_ids)
                    tasks[task_id].ancestors = ancestors

    def initialize_tasks(self):
        """Define tasks in recipe"""
        logger.info("Creating tasks from recipe")
        tasks = {}

        for diagnostic_name, diagnostic in self.diagnostics.items():
            logger.info("Creating tasks for diagnostic %s", diagnostic_name)

            # Create preprocessor tasks
            for variable_name in diagnostic['preprocessor_output']:
                task_id = diagnostic_name + TASKSEP + variable_name
                logger.info("Creating preprocessor task %s", task_id)
                task = _get_preprocessor_task(
                    variables=diagnostic['preprocessor_output'][variable_name],
                    profiles=self._preprocessors,
                    config_user=self._cfg,
                    write_ncl_interface=self._support_ncl)
                tasks[task_id] = task

            if not self._cfg['run_diagnostic']:
                continue

            # Create diagnostic tasks
            for script_name, script_cfg in diagnostic['scripts'].items():
                task_id = diagnostic_name + TASKSEP + script_name
                logger.info("Creating diagnostic task %s", task_id)
                task = DiagnosticTask(
                    script=script_cfg['script'],
                    output_dir=script_cfg['output_dir'],
                    settings=script_cfg['settings'])
                tasks[task_id] = task

        # Resolve diagnostic ancestors
        self._resolve_diagnostic_ancestors(tasks)

        # TODO: check that no loops are created (will throw RecursionError)

        # Return smallest possible set of tasks
        return get_independent_tasks(tasks.values())

    def __str__(self):
        """Get human readable summary."""
        return '\n\n'.join(str(task) for task in self.tasks)

    def run(self):
        """Run all tasks in the recipe."""
        run_tasks(
            self.tasks, max_parallel_tasks=self._cfg['max_parallel_tasks'])
