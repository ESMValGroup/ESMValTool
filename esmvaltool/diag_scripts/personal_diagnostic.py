"""
Run your own diagnostic.

This is a wrapper for your own personal diagnostic. This has to be
a Python module that you can build and run without needing to be an ESMValTool
developer. An example of such a module can be found in:
esmvaltool/diag_scripts/examples/my_personal_diagnostic_example
The initial version of this script is barebones and passes just a preprocessed
files dictionary to your personal diagnostic module.
"""
import logging
import importlib
import os
import sys
from esmvaltool.diag_scripts.shared import run_diagnostic


logger = logging.getLogger(__name__)


def _get_my_files(cfg):
    """Put files in dicts of datasets and return them."""
    files_dict = {}
    for filename, attributes in cfg['input_data'].items():
        base_file = os.path.basename(filename)
        dataset = base_file.split('_')[1]
        files_dict[dataset] = {}
        files_dict[dataset]['file'] = filename
        if 'fx_files' in attributes:
            for fx_var in attributes['fx_files']:
                files_dict[dataset][fx_var] = attributes['fx_files'][fx_var]

    return files_dict


def _import_my_diag(my_diag_place, my_diag):
    """Do an in-place import of your script of choice."""
    sys.path.insert(0, my_diag_place)
    area_package = importlib.import_module(my_diag)
    return area_package


def main(cfg):
    """Run your preferred diagnostic using preprocessed files."""
    my_files_dict = _get_my_files(cfg)
    to_run = _import_my_diag(cfg['myDiagPlace'], cfg['myDiag'])
    for analysis_function in to_run.analyses:
        analysis_function(my_files_dict)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
