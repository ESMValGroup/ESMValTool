import inspect
import logging
import os
import subprocess

import yamale

from . import preprocessor
from ._data_finder import get_start_end_year
from ._task import which

logger = logging.getLogger(__name__)


class RecipeError(Exception):
    """Recipe contains an error."""


def ncl_version():
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


def recipe_with_schema(filename):
    """Check if the recipe content matches schema."""
    schema_file = os.path.join(os.path.dirname(__file__), 'recipe_schema.yml')
    logger.debug("Checking recipe against schema %s", schema_file)
    recipe = yamale.make_data(filename)
    schema = yamale.make_schema(schema_file)
    yamale.validate(schema, recipe)


def preprocessors(profiles):
    """Check preprocessors in recipe"""
    valid_functions = set(preprocessor.DEFAULT_ORDER)
    for name, profile in profiles.items():
        invalid_functions = set(profile) - {'custom_order'} - valid_functions
        if invalid_functions:
            raise RecipeError(
                "Unknown function(s) {} in preprocessor {}, choose from: "
                "{}".format(', '.join(invalid_functions), name,
                            ', '.join(preprocessor.DEFAULT_ORDER)))


def diagnostics(diags):
    """Check diagnostics in recipe"""
    for name, diagnostic in diags.items():
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


def preprocessor_settings(settings):
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


def duplicate_datasets(datasets):
    """Check for duplicate datasets."""
    checked_datasets_ = []
    for dataset in datasets:
        if dataset in checked_datasets_:
            raise RecipeError(
                "Duplicate dataset {} in datasets section".format(dataset))
        checked_datasets_.append(dataset)


def variable(var, required_keys):
    """Check variables as derived from recipe"""
    required = set(required_keys)
    missing = required - set(var)
    if missing:
        raise RecipeError(
            "Missing keys {} from variable {} in diagnostic {}".format(
                missing, var.get('short_name'), var.get('diagnostic')))


def data_availability(input_files, var):
    """Check if the required input data is available"""
    if not input_files:
        raise RecipeError("No input files found for variable {}".format(var))

    required_years = set(range(var['start_year'], var['end_year'] + 1))
    available_years = set()
    for filename in input_files:
        start, end = get_start_end_year(filename)
        available_years.update(range(start, end + 1))

    missing_years = required_years - available_years
    if missing_years:
        raise RecipeError(
            "No input data available for years {} in files {}".format(
                ", ".join(str(year) for year in missing_years), input_files))


def multi_model_settings(all_settings):
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
