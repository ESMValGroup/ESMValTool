"""Recipe parser."""
import fnmatch
import logging
import os
from collections import OrderedDict
from copy import deepcopy

import yaml
from netCDF4 import Dataset

from . import __version__
from . import _recipe_checks as check
from ._config import TAGS, get_institutes, replace_tags
from ._data_finder import (get_input_filelist, get_input_fx_filelist,
                           get_output_file, get_statistic_output_file)
from ._provenance import TrackedFile, get_recipe_provenance
from ._recipe_checks import RecipeError
from ._task import (DiagnosticTask, get_flattened_tasks, get_independent_tasks,
                    run_tasks)
from .cmor.table import CMOR_TABLES
from .preprocessor import (DEFAULT_ORDER, FINAL_STEPS, INITIAL_STEPS,
                           MULTI_MODEL_FUNCTIONS, PreprocessingTask,
                           PreprocessorFile)
from .preprocessor._derive import get_required
from .preprocessor._download import synda_search
from .preprocessor._io import DATASET_KEYS, concatenate_callback
from .preprocessor._regrid import (get_cmor_levels, get_reference_levels,
                                   parse_cell_spec)

logger = logging.getLogger(__name__)

TASKSEP = os.sep


def ordered_safe_load(stream):
    """Load a YAML file using OrderedDict instead of dict."""
    class OrderedSafeLoader(yaml.SafeLoader):
        """Loader class that uses OrderedDict to load a map."""

    def construct_mapping(loader, node):
        """Load a map as an OrderedDict."""
        loader.flatten_mapping(node)
        return OrderedDict(loader.construct_pairs(node))

    OrderedSafeLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)

    return yaml.load(stream, OrderedSafeLoader)


def load_raw_recipe(filename):
    """Check a recipe file and return it in raw form."""
    # Note that many checks can only be performed after the automatically
    # computed entries have been filled in by creating a Recipe object.
    check.recipe_with_schema(filename)
    with open(filename, 'r') as file:
        contents = file.read()
        raw_recipe = yaml.safe_load(contents)
        raw_recipe['preprocessors'] = ordered_safe_load(contents).get(
            'preprocessors', {})

    check.diagnostics(raw_recipe['diagnostics'])
    return raw_recipe


def read_recipe_file(filename, config_user, initialize_tasks=True):
    """Read a recipe from file."""
    raw_recipe = load_raw_recipe(filename)
    return Recipe(
        raw_recipe, config_user, initialize_tasks, recipe_file=filename)


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


def _add_cmor_info(variable, override=False):
    """Add information from CMOR tables to variable."""
    logger.debug("If not present: adding keys from CMOR table to %s", variable)

    if 'cmor_table' not in variable or 'mip' not in variable:
        logger.debug("Skipping because cmor_table or mip not specified")
        return

    if variable['cmor_table'] not in CMOR_TABLES:
        logger.warning("Unknown CMOR table %s", variable['cmor_table'])

    derive = variable.get('derive', False)
    # Copy the following keys from CMOR table
    cmor_keys = [
        'standard_name', 'long_name', 'units', 'modeling_realm', 'frequency'
    ]
    cmor_table = variable['cmor_table']
    mip = variable['mip']
    short_name = variable['short_name']
    table_entry = CMOR_TABLES[cmor_table].get_variable(mip, short_name)

    if derive and table_entry is None:
        custom_table = CMOR_TABLES['custom']
        table_entry = custom_table.get_variable(mip, short_name)

    if table_entry is None:
        raise RecipeError(
            "Unable to load CMOR table '{}' for variable '{}' with mip '{}'".
            format(cmor_table, short_name, mip))

    mip_info = CMOR_TABLES[cmor_table].get_table(mip)
    if mip_info:
        table_entry.frequency = mip_info.frequency

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
    check.variable(variable, required_keys=cmor_keys)


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
            settings['extract_levels']['levels'] = get_reference_levels(
                filename, variable_data['project'], dataset,
                variable_data['short_name'],
                os.path.splitext(variable_data['filename'])[0] + '_fixed')


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
    else:
        # Check that MxN grid spec is correct
        parse_cell_spec(settings['regrid']['target_grid'])


def _update_regrid_time(variable, settings):
    """Input data frequency automatically for regrid_time preprocessor."""
    regrid_time = settings.get('regrid_time')
    if regrid_time is None:
        return
    frequency = settings.get('regrid_time', {}).get('frequency')
    if not frequency:
        settings['regrid_time']['frequency'] = variable['frequency']


def _get_dataset_info(dataset, variables):
    for var in variables:
        if var['dataset'] == dataset:
            return var
    raise RecipeError("Unable to find matching file for dataset"
                      "{}".format(dataset))


def _augment(base, update):
    """Update dict base with values from dict update."""
    for key in update:
        if key not in base:
            base[key] = update[key]


def _dataset_to_file(variable, config_user):
    """Find the first file belonging to dataset from variable info."""
    files = get_input_filelist(
        variable=variable,
        rootpath=config_user['rootpath'],
        drs=config_user['drs'])
    if not files and variable.get('derive'):
        first_required = get_required(variable['short_name'])[0]
        _augment(first_required, variable)
        files = get_input_filelist(
            variable=first_required,
            rootpath=config_user['rootpath'],
            drs=config_user['drs'])
    check.data_availability(files, variable)
    return files[0]


def _limit_datasets(variables, profile, max_datasets=0):
    """Try to limit the number of datasets to max_datasets."""
    if not max_datasets:
        return variables

    logger.info("Limiting the number of datasets to %s", max_datasets)

    required_datasets = [
        (profile.get('extract_levels') or {}).get('levels'),
        (profile.get('regrid') or {}).get('target_grid'),
        variables[0].get('reference_dataset'),
        variables[0].get('alternative_dataset'),
    ]

    limited = [v for v in variables if v['dataset'] in required_datasets]
    for variable in variables:
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
        download_folder = os.path.join(config_user['preproc_dir'], 'downloads')
        settings['download'] = {
            'dest_folder': download_folder,
        }

    # Configure loading
    settings['load'] = {
        'callback': concatenate_callback,
    }
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
    settings['fix_file'] = dict(fix)
    settings['fix_file']['output_dir'] = fix_dir
    # Cube fixes
    # Only supply mip if the CMOR check fixes are implemented.
    if variable.get('cmor_table'):
        fix['cmor_table'] = variable['cmor_table']
        fix['mip'] = variable['mip']
        fix['frequency'] = variable['frequency']
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
        settings['derive'] = {
            'short_name': variable['short_name'],
            'standard_name': variable['standard_name'],
            'long_name': variable['long_name'],
            'units': variable['units'],
        }

    # Configure CMOR metadata check
    if variable.get('cmor_table'):
        settings['cmor_check_metadata'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
            'frequency': variable['frequency'],
        }
    # Configure final CMOR data check
    if variable.get('cmor_table'):
        settings['cmor_check_data'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
            'frequency': variable['frequency'],
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
    """Find and set the FX derive/mask settings."""
    # update for derive
    if 'derive' in settings:
        fx_files = {}
        for var in get_required(variable['short_name']):
            if 'fx_files' in var:
                _augment(var, variable)
                fx_files.update(
                    get_input_fx_filelist(
                        variable=var,
                        rootpath=config_user['rootpath'],
                        drs=config_user['drs']))
        settings['derive']['fx_files'] = fx_files

    # update for landsea
    if 'mask_landsea' in settings:
        # Configure ingestion of land/sea masks
        logger.debug('Getting fx mask settings now...')

        settings['mask_landsea']['fx_files'] = []

        var = dict(variable)
        var['fx_files'] = ['sftlf', 'sftof']
        fx_files_dict = get_input_fx_filelist(
            variable=var,
            rootpath=config_user['rootpath'],
            drs=config_user['drs'])

        # allow both sftlf and sftof
        if fx_files_dict['sftlf']:
            settings['mask_landsea']['fx_files'].append(fx_files_dict['sftlf'])
        if fx_files_dict['sftof']:
            settings['mask_landsea']['fx_files'].append(fx_files_dict['sftof'])

    if 'mask_landseaice' in settings:
        logger.debug('Getting fx mask settings now...')

        settings['mask_landseaice']['fx_files'] = []

        var = dict(variable)
        var['fx_files'] = ['sftgif']
        fx_files_dict = get_input_fx_filelist(
            variable=var,
            rootpath=config_user['rootpath'],
            drs=config_user['drs'])

        # allow sftgif (only, for now)
        if fx_files_dict['sftgif']:
            settings['mask_landseaice']['fx_files'].append(
                fx_files_dict['sftgif'])

    for step in ('average_region', 'average_volume'):
        if settings.get(step, {}).get('fx_files'):
            settings[step]['fx_files'] = get_input_fx_filelist(
                variable=variable,
                rootpath=config_user['rootpath'],
                drs=config_user['drs'],
            )


def _read_attributes(filename):
    """Read the attributes from a netcdf file."""
    attributes = {}
    if not (os.path.exists(filename)
            and os.path.splitext(filename)[1].lower() == '.nc'):
        return attributes

    with Dataset(filename, 'r') as dataset:
        for attr in dataset.ncattrs():
            attributes[attr] = getattr(dataset, attr)
    return attributes


def _get_input_files(variable, config_user):
    """Get the input files for a single dataset."""
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
    if (not config_user.get('skip-nonexistent')
            or variable['dataset'] == variable.get('reference_dataset')):
        check.data_availability(input_files, variable)

    # Set up provenance tracking
    for i, filename in enumerate(input_files):
        attributes = _read_attributes(filename)
        input_files[i] = TrackedFile(filename, attributes)

    return input_files


def _apply_preprocessor_profile(settings, profile_settings):
    """Apply settings from preprocessor profile."""
    profile_settings = deepcopy(profile_settings)
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


def _get_statistic_attributes(products):
    """Get attributes for the statistic output products."""
    attributes = {}
    some_product = next(iter(products))
    for key, value in some_product.attributes.items():
        if all(p.attributes.get(key, object()) == value for p in products):
            attributes[key] = value

    # Ensure start_year and end_year attributes are available
    for product in products:
        start = product.attributes['start_year']
        if 'start_year' not in attributes or start < attributes['start_year']:
            attributes['start_year'] = start
        end = product.attributes['end_year']
        if 'end_year' not in attributes or end > attributes['end_year']:
            attributes['end_year'] = end

    return attributes


def _get_remaining_common_settings(step, order, products):
    """Get preprocessor settings that are shared between products."""
    settings = {}
    remaining_steps = order[order.index(step) + 1:]
    some_product = next(iter(products))
    for key, value in some_product.settings.items():
        if key in remaining_steps:
            if all(p.settings.get(key, object()) == value for p in products):
                settings[key] = value
    return settings


def _update_multi_dataset_settings(variable, settings):
    """Configure multi dataset statistics."""
    for step in MULTI_MODEL_FUNCTIONS:
        if not settings.get(step):
            continue
        # Exclude dataset if requested
        exclude = {
            _special_name_to_dataset(variable, dataset)
            for dataset in settings[step].pop('exclude', [])
        }
        if variable['dataset'] in exclude:
            settings.pop(step)


def _update_statistic_settings(products, order, preproc_dir):
    """Define statistic output products."""
    # TODO: move this to multi model statistics function?
    # But how to check, with a dry-run option?
    step = 'multi_model_statistics'

    products = {p for p in products if step in p.settings}
    if not products:
        return

    some_product = next(iter(products))
    for statistic in some_product.settings[step]['statistics']:
        attributes = _get_statistic_attributes(products)
        attributes['dataset'] = 'MultiModel{}'.format(statistic.title())
        attributes['filename'] = get_statistic_output_file(
            attributes, preproc_dir)
        common_settings = _get_remaining_common_settings(step, order, products)
        statistic_product = PreprocessorFile(attributes, common_settings)
        for product in products:
            settings = product.settings[step]
            if 'output_products' not in settings:
                settings['output_products'] = {}
            settings['output_products'][statistic] = statistic_product


def _match_products(products, variables):
    """Match a list of input products to output product attributes."""
    grouped_products = {}

    def get_matching(attributes):
        """Find the output filename which matches input attributes best."""
        score = 0
        filenames = []
        for variable in variables:
            filename = variable['filename']
            tmp = sum(v == variable.get(k) for k, v in attributes.items())
            if tmp > score:
                score = tmp
                filenames = [filename]
            elif tmp == score:
                filenames.append(filename)
        if not filenames:
            logger.warning(
                "Unable to find matching output file for input file %s",
                filename)
        return filenames

    # Group input files by output file
    for product in products:
        for filename in get_matching(product.attributes):
            if filename not in grouped_products:
                grouped_products[filename] = []
            grouped_products[filename].append(product)

    return grouped_products


def _get_preprocessor_products(variables, profile, order, ancestor_products,
                               config_user):
    """Get preprocessor product definitions for a set of datasets."""
    products = set()

    for variable in variables:
        variable['filename'] = get_output_file(variable,
                                               config_user['preproc_dir'])

    if ancestor_products:
        grouped_ancestors = _match_products(ancestor_products, variables)
    else:
        grouped_ancestors = {}

    for variable in variables:
        settings = _get_default_settings(
            variable, config_user, derive='derive' in profile)
        _apply_preprocessor_profile(settings, profile)
        _update_multi_dataset_settings(variable, settings)
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
        _update_regrid_time(variable, settings)
        ancestors = grouped_ancestors.get(variable['filename'])
        if not ancestors:
            ancestors = _get_input_files(variable, config_user)
            if config_user.get('skip-nonexistent') and not ancestors:
                logger.info("Skipping: no data found for %s", variable)
                continue
        product = PreprocessorFile(
            attributes=variable, settings=settings, ancestors=ancestors)
        products.add(product)

    _update_statistic_settings(products, order, config_user['preproc_dir'])

    for product in products:
        product.check()

    return products


def _get_single_preprocessor_task(variables,
                                  profile,
                                  config_user,
                                  name,
                                  ancestor_tasks=None):
    """Create preprocessor tasks for a set of datasets."""
    if ancestor_tasks is None:
        ancestor_tasks = []
    order = _extract_preprocessor_order(profile)
    ancestor_products = [p for task in ancestor_tasks for p in task.products]
    products = _get_preprocessor_products(
        variables=variables,
        profile=profile,
        order=order,
        ancestor_products=ancestor_products,
        config_user=config_user,
    )

    if not products:
        raise RecipeError(
            "Did not find any input data for task {}".format(name))

    task = PreprocessingTask(
        products=products,
        ancestors=ancestor_tasks,
        name=name,
        order=order,
        debug=config_user['save_intermediary_cubes'],
        write_ncl_interface=config_user['write_ncl_interface'],
    )

    logger.info("PreprocessingTask %s created. It will create the files:\n%s",
                task.name, '\n'.join(p.filename for p in task.products))

    return task


def _extract_preprocessor_order(profile):
    """Extract the order of the preprocessing steps from the profile."""
    custom_order = profile.pop('custom_order', False)
    if not custom_order:
        return DEFAULT_ORDER
    order = tuple(p for p in profile if p not in INITIAL_STEPS + FINAL_STEPS)
    return INITIAL_STEPS + order + FINAL_STEPS


def _split_settings(settings, step, order=DEFAULT_ORDER):
    """Split settings, using step as a separator."""
    before = {}
    for _step in order:
        if _step == step:
            break
        if _step in settings:
            before[_step] = settings[_step]
    after = {
        k: v
        for k, v in settings.items() if not (k == step or k in before)
    }
    return before, after


def _split_derive_profile(profile):
    """Split the derive preprocessor profile."""
    order = _extract_preprocessor_order(profile)
    before, after = _split_settings(profile, 'derive', order)
    after['derive'] = True
    after['fix_file'] = False
    after['fix_metadata'] = False
    after['fix_data'] = False
    if order != DEFAULT_ORDER:
        before['custom_order'] = True
        after['custom_order'] = True
    return before, after


def _get_derive_input_variables(variables, config_user):
    """Determine the input sets of `variables` needed for deriving."""
    derive_input = {}

    def append(group_prefix, var):
        """Append variable `var` to a derive input group."""
        group = group_prefix + var['short_name']
        var['variable_group'] = group
        if group not in derive_input:
            derive_input[group] = []
        derive_input[group].append(var)

    for variable in variables:
        group_prefix = variable['variable_group'] + '_derive_input_'
        if not variable.get('force_derivation') and get_input_filelist(
                variable=variable,
                rootpath=config_user['rootpath'],
                drs=config_user['drs']):
            # No need to derive, just process normally up to derive step
            var = deepcopy(variable)
            append(group_prefix, var)
        else:
            # Process input data needed to derive variable
            for var in get_required(variable['short_name']):
                _augment(var, variable)
                append(group_prefix, var)

    return derive_input


def _get_preprocessor_task(variables, profiles, config_user, task_name):
    """Create preprocessor task(s) for a set of datasets."""
    # First set up the preprocessor profile
    variable = variables[0]
    preproc_name = variable.get('preprocessor')
    if preproc_name not in profiles:
        raise RecipeError(
            "Unknown preprocessor {} in variable {} of diagnostic {}".format(
                preproc_name, variable['short_name'], variable['diagnostic']))
    profile = deepcopy(profiles[variable['preprocessor']])
    logger.info("Creating preprocessor '%s' task for variable '%s'",
                variable['preprocessor'], variable['short_name'])
    variables = _limit_datasets(variables, profile,
                                config_user.get('max_datasets'))
    for variable in variables:
        _add_cmor_info(variable)
    # Create preprocessor task(s)
    derive_tasks = []
    if variable.get('derive'):
        # Create tasks to prepare the input data for the derive step
        derive_profile, profile = _split_derive_profile(profile)
        derive_input = _get_derive_input_variables(variables, config_user)

        for derive_variables in derive_input.values():
            for derive_variable in derive_variables:
                _add_cmor_info(derive_variable, override=True)
            derive_name = task_name.split(
                TASKSEP)[0] + TASKSEP + derive_variables[0]['variable_group']
            task = _get_single_preprocessor_task(
                derive_variables,
                derive_profile,
                config_user,
                name=derive_name)
            derive_tasks.append(task)

    # Create (final) preprocessor task
    task = _get_single_preprocessor_task(
        variables,
        profile,
        config_user,
        ancestor_tasks=derive_tasks,
        name=task_name)

    return task


class Recipe:
    """Recipe object."""

    def __init__(self,
                 raw_recipe,
                 config_user,
                 initialize_tasks=True,
                 recipe_file=None):
        """Parse a recipe file into an object."""
        self._cfg = deepcopy(config_user)
        self._cfg['write_ncl_interface'] = self._need_ncl(
            raw_recipe['diagnostics'])
        self._filename = os.path.basename(recipe_file)
        self._preprocessors = raw_recipe.get('preprocessors', {})
        if 'default' not in self._preprocessors:
            self._preprocessors['default'] = {}
        self.diagnostics = self._initialize_diagnostics(
            raw_recipe['diagnostics'], raw_recipe.get('datasets', []))
        self.entity = self._initalize_provenance(
            raw_recipe.get('documentation', {}))
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
                    check.ncl_version()
                    return True
        return False

    def _initalize_provenance(self, raw_documentation):
        """Initialize the recipe provenance."""
        doc = deepcopy(raw_documentation)
        for key in doc:
            if key in TAGS:
                doc[key] = replace_tags(key, doc[key])

        return get_recipe_provenance(doc, self._filename)

    def _initialize_diagnostics(self, raw_diagnostics, raw_datasets):
        """Define diagnostics in recipe."""
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
            for key in ('themes', 'realms'):
                if key in raw_diagnostic:
                    for script in diagnostic['scripts'].values():
                        script['settings'][key] = raw_diagnostic[key]
            diagnostics[name] = diagnostic

        return diagnostics

    @staticmethod
    def _initialize_datasets(raw_datasets):
        """Define datasets used by variable."""
        datasets = deepcopy(raw_datasets)

        for dataset in datasets:
            for key in dataset:
                DATASET_KEYS.add(key)

        check.duplicate_datasets(datasets)
        return datasets

    def _initialize_variables(self, raw_variable, raw_datasets):
        """Define variables for all datasets."""
        variables = []

        raw_variable = deepcopy(raw_variable)
        datasets = self._initialize_datasets(
            raw_datasets + raw_variable.pop('additional_datasets', []))

        for index, dataset in enumerate(datasets):
            variable = deepcopy(raw_variable)
            variable.update(dataset)
            variable['recipe_dataset_index'] = index
            if ('cmor_table' not in variable
                    and variable.get('project') in CMOR_TABLES):
                variable['cmor_table'] = variable['project']
            if 'end_year' in variable and 'max_years' in self._cfg:
                variable['end_year'] = min(
                    variable['end_year'],
                    variable['start_year'] + self._cfg['max_years'] - 1)
            variables.append(variable)

        required_keys = {
            'short_name',
            'mip',
            'dataset',
            'project',
            'start_year',
            'end_year',
            'preprocessor',
            'diagnostic',
        }

        for variable in variables:
            _update_from_others(variable, ['cmor_table', 'mip'], datasets)
            institute = get_institutes(variable)
            if institute:
                variable['institute'] = institute
            check.variable(variable, required_keys)
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
        """Define variables in diagnostic."""
        logger.debug("Populating list of variables for diagnostic %s",
                     diagnostic_name)

        preprocessor_output = {}

        for variable_group, raw_variable in raw_variables.items():
            if raw_variable is None:
                raw_variable = {}
            else:
                raw_variable = deepcopy(raw_variable)
            raw_variable['variable_group'] = variable_group
            if 'short_name' not in raw_variable:
                raw_variable['short_name'] = variable_group
            raw_variable['diagnostic'] = diagnostic_name
            raw_variable['preprocessor'] = str(
                raw_variable.get('preprocessor', 'default'))
            preprocessor_output[variable_group] = \
                self._initialize_variables(raw_variable, raw_datasets)

        return preprocessor_output

    def _initialize_scripts(self, diagnostic_name, raw_scripts,
                            variable_names):
        """Define script in diagnostic."""
        if not raw_scripts:
            return {}

        logger.debug("Setting script for diagnostic %s", diagnostic_name)

        scripts = {}

        for script_name, raw_settings in raw_scripts.items():
            settings = deepcopy(raw_settings)
            script = settings.pop('script')
            ancestors = []
            for id_glob in settings.pop('ancestors', variable_names):
                if TASKSEP not in id_glob:
                    id_glob = diagnostic_name + TASKSEP + id_glob
                ancestors.append(id_glob)
            settings['recipe'] = self._filename
            settings['version'] = __version__
            settings['script'] = script_name
            # Add output dirs to settings
            for dir_name in ('run_dir', 'plot_dir', 'work_dir'):
                settings[dir_name] = os.path.join(self._cfg[dir_name],
                                                  diagnostic_name, script_name)
            # Copy other settings
            if self._cfg['write_ncl_interface']:
                settings['exit_on_ncl_warning'] = self._cfg['exit_on_warning']
            for key in (
                    'max_data_filesize',
                    'output_file_type',
                    'log_level',
                    'write_plots',
                    'write_netcdf',
                    'profile_diagnostic',
                    'auxiliary_data_dir',
            ):
                settings[key] = self._cfg[key]

            scripts[script_name] = {
                'script': script,
                'output_dir': settings['work_dir'],
                'settings': settings,
                'ancestors': ancestors,
            }

        return scripts

    def _resolve_diagnostic_ancestors(self, tasks):
        """Resolve diagnostic ancestors."""
        tasks = {t.name: t for t in tasks}
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
                                "Could not find any ancestors matching {}".
                                format(id_glob))
                        logger.debug("Pattern %s matches %s", id_glob,
                                     ancestor_ids)
                        ancestors.extend(tasks[a] for a in ancestor_ids)
                    tasks[task_id].ancestors = ancestors

    def initialize_tasks(self):
        """Define tasks in recipe."""
        logger.info("Creating tasks from recipe")
        tasks = set()

        for diagnostic_name, diagnostic in self.diagnostics.items():
            logger.info("Creating tasks for diagnostic %s", diagnostic_name)

            # Create preprocessor tasks
            for variable_group in diagnostic['preprocessor_output']:
                task_name = diagnostic_name + TASKSEP + variable_group
                logger.info("Creating preprocessor task %s", task_name)
                task = _get_preprocessor_task(
                    variables=diagnostic['preprocessor_output']
                    [variable_group],
                    profiles=self._preprocessors,
                    config_user=self._cfg,
                    task_name=task_name)
                tasks.add(task)

            # Create diagnostic tasks
            for script_name, script_cfg in diagnostic['scripts'].items():
                task_name = diagnostic_name + TASKSEP + script_name
                logger.info("Creating diagnostic task %s", task_name)
                task = DiagnosticTask(
                    script=script_cfg['script'],
                    output_dir=script_cfg['output_dir'],
                    settings=script_cfg['settings'],
                    name=task_name)
                tasks.add(task)

        check.tasks_valid(tasks)

        # Resolve diagnostic ancestors
        self._resolve_diagnostic_ancestors(tasks)

        # Select only requested tasks
        tasks = get_flattened_tasks(tasks)
        if not self._cfg.get('run_diagnostic'):
            tasks = {t for t in tasks if isinstance(t, PreprocessingTask)}
        if self._cfg.get('diagnostics'):
            names = {t.name for t in tasks}
            selection = set()
            for pattern in self._cfg.get('diagnostics'):
                selection |= set(fnmatch.filter(names, pattern))
            tasks = {t for t in tasks if t.name in selection}

        tasks = get_flattened_tasks(tasks)
        logger.info("These tasks will be executed: %s",
                    ', '.join(t.name for t in tasks))

        # Initialize task provenance
        for task in tasks:
            task.initialize_provenance(self.entity)

        # TODO: check that no loops are created (will throw RecursionError)

        # Return smallest possible set of tasks
        return get_independent_tasks(tasks)

    def __str__(self):
        """Get human readable summary."""
        return '\n\n'.join(str(task) for task in self.tasks)

    def run(self):
        """Run all tasks in the recipe."""
        run_tasks(
            self.tasks, max_parallel_tasks=self._cfg['max_parallel_tasks'])
