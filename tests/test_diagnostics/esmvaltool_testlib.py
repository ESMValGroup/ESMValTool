"""
This modules provides generic functions needed for ESMValTool
performance metrics testing
"""

import glob
import os
import shutil
import tempfile
from unittest import SkipTest

import numpy as np
import yaml
from easytest import EasyTest

import esmvaltool


def _load_config(filename=None):
    """ Load test configuration
    """
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

    cfg['reference']['output'] = os.path.abspath(
        os.path.expanduser(cfg['reference']['output']))

    return cfg


_CFG = _load_config()


def _create_config_user_file(output_directory):
    """ Write a config-user.yml file for running ESMValTool
        such that it writes all output to `output_directory`.
    """

    cfg = _CFG['user']

    # update/insert output directory definitions
    cfg['run_dir'] = output_directory
    for task in ('preproc', 'work', 'plot'):
        cfg[task + '_dir'] = os.path.join(output_directory, task)

    # write to file
    fd, filename = tempfile.mkstemp(text=True)
    with os.fdopen(fd, 'w') as file:
        yaml.safe_dump(cfg, file)
    # use code below once moving run_dir is fixed
    #     filename = os.path.join(output_directory, 'config-user.yml')
    #     with open(filename, 'w') as file:
    #         yaml.safe_dump(cfg, file)

    return filename


def _get_ref_directory(namelist):
    """ Get the name of the directory containing the reference data"""

    namelist_base = os.path.splitext(os.path.basename(namelist))[0]
    namelist_glob = os.path.join(_CFG['reference']['output'],
                                 '*' + namelist_base + '*')
    refs = glob.glob(namelist_glob)
    refdirectory = refs[0] if refs else None

    if len(refs) > 1:
        print("Warning: found multiple possible reference datasets:\n{}\n"
              "for namelist {}\n"
              "Using only the first one, {}".format("\n".join(refs), namelist,
                                                    refdirectory))
    return refdirectory


class ESMValToolTest(EasyTest):
    """
    main class for all ESMValTool tests
    """

    def __init__(self, namelist, output_directory, ignore='', **kwargs):
        """
        namelist: str
            The filename of the namelist that should be tested.
        output_directory : str
            The name of a temporary directory where results can be stored.
        ignore: str or iterable of str
            Glob patterns of files to be ignored when testing.
        """
        if not _CFG['test']['run']:
            raise SkipTest("System tests disabled")

        self.ignore = (ignore, ) if isinstance(ignore, str) else ignore

        script_root = esmvaltool.get_script_root()

        # Normalize input paths
        if not os.path.exists(namelist):
            namelist = os.path.join(script_root, 'nml', namelist)
        namelist = os.path.abspath(namelist)
        output_directory = os.path.abspath(output_directory)

        # Are we asked to generate reference data?
        self.generate_ref = _CFG['reference'].get('generate', False)

        # Find the location of the reference data
        refdirectory = _get_ref_directory(namelist)

        # If no reference data is available
        if refdirectory is None:

            # Should reference data be generated?
            if not self.generate_ref:
                raise SkipTest("No reference data available")

            # Should it be generated for this namelist?
            generate_nml = _CFG['reference'].get('generate_namelists', [])
            namelist_base = os.path.basename(namelist)
            if not (generate_nml == [] or namelist_base in generate_nml):
                raise SkipTest("No reference data available")

            # Define a new reference data directory
            refdirectory = os.path.join(_CFG['reference']['output'],
                                        os.path.splitext(namelist_base)[0])

        # the command to call ESMValTool
        exe = os.path.join(script_root, 'main.py')

        # required ESMValTool configuration file
        config = _create_config_user_file(output_directory)

        super(ESMValToolTest, self).__init__(
            exe=exe,
            args=['-n', namelist, '-c', config],
            output_directory=output_directory,
            refdirectory=refdirectory,
            **kwargs)

    def run(self, **kwargs):
        """ Run tests, unless we are asked to generate the reference data instead.
        """
        if self.generate_ref:
            self.generate_reference_data()
            raise SkipTest("Generated reference data instead of running test")
        else:
            super(ESMValToolTest, self).run_tests(**kwargs)

    def generate_reference_data(self):
        """ Generate reference data by executing the namelist and then copying
            results to the output directory.
        """
        if not os.path.exists(self.refdirectory):
            self._execute()
            shutil.copytree(self.output_directory, self.refdirectory)
        else:
            print("Warning: not generating reference data, reference "
                  "directory {} already exists.".format(self.refdirectory))

    # Overwrite this method of easytest.EasyTest to be able to ignore certain files
    def _get_files_from_refdir(self):
        """ Get a list of files from reference directory, while ignoring files
            that match patterns in self.ignore.
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

    # Overwrite this method of easytest.EasyTest because it is broken
    # for the case where x1 and x2 have no length
    def _compare_netcdf_values(self, F1, F2, allow_subset=False):
        """
        compare if two netCDF files have the same values
        """
        for k in F1.variables.keys():
            #             print('Comparing variable', k)
            x1 = F1.variables[k][:]
            x2 = F2.variables[k][:]

            if allow_subset:  # allow that only a subset of data is compared
                raise NotImplementedError
            else:
                if not np.array_equal(x1, x2):
                    return False

        return True
