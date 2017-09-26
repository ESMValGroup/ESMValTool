"""
This modules provides generic functions needed for ESMValTool
performance metrics testing
"""

import glob
import os
import shutil

import yaml
from easytest import EasyTest

import esmvaltool


def _load_config(filename=None):
    """ Load test configuration
    """
    if filename is None:
        config_file = 'config-test.yml'
        if os.path.isfile(config_file):
            filename = config_file
        else:
            filename = os.path.join(os.path.dirname(__file__), config_file)

    with open(filename, 'r') as file:
        cfg = yaml.safe_load(file)

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
    filename = os.path.join(output_directory, 'config-user.yml')
    with open(filename, 'w') as file:
        yaml.safe_dump(cfg, file)

    return filename


class ESMValToolTest(EasyTest):
    """
    main class for all ESMValTool tests
    """

    def __init__(self, namelist, output_directory, **kwargs):
        """
        namelist: str
            The filename of the namelist that should be tested.
        output_directory : str
            The name of a temporary directory where results can be stored.
        """
        script_root = esmvaltool.get_script_root()

        # the command to call ESMValTool
        exe = os.path.join(script_root, 'main.py')

        # required ESMValTool configuration files
        config = _create_config_user_file(output_directory)
        if not os.path.isfile(namelist):
            namelist = os.path.join(script_root, 'nml', namelist)

        # the location of the reference output
        namelist_base = os.path.splitext(os.path.basename(namelist))[0]
        refs = glob.glob(
            os.path.join(_CFG['reference_output'], '*' + namelist_base + '*'))
        refdirectory = None if not refs else refs[0]
        if len(refs) > 1:
            print("Warning: found multiple possible reference datasets:\n{}"
                  "for namelist {}\n"
                  "Using only {}".format("\n".join(refs), namelist,
                                         refdirectory))

        super(ESMValToolTest, self).__init__(
            exe=exe,
            args=['-n', namelist, '-c', config],
            output_directory=output_directory,
            refdirectory=refdirectory,
            checksum_exclude=['ps', 'epsi', 'eps'],
            **kwargs)

    def generate_reference_data(self):
        """
        generate reference data by executing the namelist once and then copy
        results to the output directory
        """
        if not os.path.exists(self.refdirectory):
            self._execute()
            shutil.copytree(self.output_directory, self.refdirectory)
