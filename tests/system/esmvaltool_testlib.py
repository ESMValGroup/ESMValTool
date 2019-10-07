"""Provide a class for testing esmvaltool."""

import glob
import os
import shutil
import sys
from unittest import SkipTest

import numpy as np
import yaml
from easytest import EasyTest

import esmvaltool


def _load_config(filename=None):
    """Load test configuration"""
    if filename is None:
        # look in default locations for config-test.yml
        config_file = 'config-test.yml'
        default_locations = [
            '.',
            '~/.esmvaltool',
            os.path.dirname(__file__),
        ]
        for path in default_locations:
            filepath = os.path.join(os.path.expanduser(path), config_file)
            if os.path.exists(filepath):
                filename = os.path.abspath(filepath)
                break

    with open(filename, 'r') as file:
        cfg = yaml.safe_load(file)

    cfg['configfile'] = filename
    cfg['reference']['output'] = os.path.abspath(
        os.path.expanduser(cfg['reference']['output']))

    if cfg['test'].get('recipes', []) == []:
        script_root = esmvaltool.get_script_root()
        recipe_glob = os.path.join(script_root, 'nml', 'recipe_*.yml')
        cfg['test']['recipes'] = glob.glob(recipe_glob)

    return cfg


_CFG = _load_config()

RECIPES = _CFG['test']['recipes']


def _create_config_user_file(output_directory):
    """Write a config-user.yml file.

    Write a configuration file for running ESMValTool
    such that it writes all output to `output_directory`.
    """
    cfg = _CFG['user']

    cfg['output_dir'] = output_directory

    # write to file
    filename = os.path.join(output_directory, 'config-user.yml')
    with open(filename, 'w') as file:
        yaml.safe_dump(cfg, file)

    return filename


class ESMValToolTest(EasyTest):
    """Main class for ESMValTool test runs."""

    def __init__(self, recipe, output_directory, ignore='', **kwargs):
        """
        Create ESMValToolTest instance

        recipe: str
            The filename of the recipe that should be tested.
        output_directory : str
            The name of a directory where results can be stored.
        ignore: str or iterable of str
            Glob patterns of files to be ignored when testing.
        """
        if not _CFG['test']['run']:
            raise SkipTest("System tests disabled in {}".format(
                _CFG['configfile']))

        self.ignore = (ignore, ) if isinstance(ignore, str) else ignore

        script_root = esmvaltool.get_script_root()

        # Set recipe path
        if not os.path.exists(recipe):
            recipe = os.path.join(
                os.path.dirname(script_root), 'recipes', recipe)
        self.recipe_file = os.path.abspath(recipe)

        # Simulate input data?
        self.simulate_input = _CFG['test']['simulate_input']

        # Create reference output?
        self.create_reference_output = _CFG['reference']['generate']

        # Define reference output path
        reference_dir = os.path.join(
            _CFG['reference']['output'],
            os.path.splitext(os.path.basename(self.recipe_file))[0])

        # If reference data is neither available nor should be generated, skip
        if not (os.path.exists(reference_dir) or self.create_reference_output):
            raise SkipTest(
                "No reference data available for recipe {} in {}".format(
                    recipe, _CFG['reference']['output']))

        # Write ESMValTool configuration file
        self.config_user_file = _create_config_user_file(output_directory)

        super(ESMValToolTest, self).__init__(
            exe='esmvaltool',
            args=['-n', self.recipe_file, '-c', self.config_user_file],
            output_directory=output_directory,
            refdirectory=reference_dir,
            **kwargs)

    def run(self, **kwargs):
        """Run tests or generate reference data."""
        if self.simulate_input:
            from .data_simulator import simulate_input_data
            simulate_input_data(
                recipe_file=self.recipe_file,
                config_user_file=self.config_user_file)

        if self.create_reference_output:
            self.generate_reference_output()
            raise SkipTest("Generated reference data instead of running test")
        else:
            super(ESMValToolTest, self).run_tests(**kwargs)

    def generate_reference_output(self):
        """Generate reference output.

        Generate reference data by executing the recipe and then moving
        results to the output directory.
        """
        if not os.path.exists(self.refdirectory):
            self._execute()
            shutil.move(self.output_directory,
                        os.path.dirname(self.refdirectory))
        else:
            print("Warning: not generating reference data, reference "
                  "directory {} already exists.".format(self.refdirectory))

    def _execute(self):
        """Execute ESMValTool

        Override the _execute method because we want to run in our own
        Python instance to get coverage reporting and we want to update
        the location of `self.output_directory` afterwards.
        """
        # run ESMValTool
        sys.argv[1:] = self.args
        esmvaltool.main.run()

        # Update the output directory to point to the output of the run
        output_directory = self.output_directory  # noqa

        output = []
        for path in os.listdir(output_directory):
            path = os.path.join(output_directory, path)
            if os.path.isdir(path):
                output.append(path)

        if not output:
            raise OSError(
                "Output directory not found in location {}. "
                "Probably ESMValTool failed to create any output.".format(
                    output_directory))

        if len(output) > 1:
            print("Warning: found multiple output directories:\n{}\nin output "
                  "location {}\nusing the first one.".format(
                      output, output_directory))

        self.output_directory = output[0] + os.sep  # noqa

    def _get_files_from_refdir(self):
        """Get a list of files from reference directory.

        Ignore files that match patterns in self.ignore.

        Override this method of easytest.EasyTest to be able to ignore certain
        files.
        """
        from fnmatch import fnmatchcase

        matches = []
        for root, _, filenames in os.walk(self.refdirectory):
            for filename in filenames:
                path = os.path.join(root, filename)
                relpath = os.path.relpath(path, start=self.refdirectory)
                for pattern in self.ignore:
                    if fnmatchcase(relpath, pattern):
                        break
                else:
                    matches.append(path)

        return matches

    def _compare_netcdf_values(self, f1, f2, allow_subset=False):
        """Compare two netCDF4 Dataset instances.

        Check if dataset2 contains the same variable values as dataset1.

        Override this method of easytest.EasyTest because it is broken
        for the case where value1 and value2 have no length.
        """
        if allow_subset:  # allow that only a subset of data is compared
            raise NotImplementedError

        for key in f1.variables:
            values1 = f1.variables[key][:]
            values2 = f2.variables[key][:]

            if not np.array_equal(values1, values2):
                return False

        return True
