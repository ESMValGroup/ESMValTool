"""
Run your own diagnostic.

This is a wrapper for your own personal diagnostic. This has to be
a Python module that you can build and run without needing to be an ESMValTool
developer. An example of such a module can be found in:
esmvaltool/diag_scripts/examples/my_personal_diagnostic_example

This script is barebones and passes the config file
to your personal diagnostic module. This wrapper should not do
more than this since it only acts as a bridge between the preprocessing
suite and the personal diagnostic that, in turn, may use the full
capacity of preprocessor.
"""
import logging
import importlib
import sys
from esmvaltool.diag_scripts.shared import run_diagnostic


logger = logging.getLogger(__name__)


def _import_my_diag(my_diag_place, my_diag):
    """Do an in-place import of your script of choice."""
    sys.path.insert(0, my_diag_place)
    area_package = importlib.import_module(my_diag)
    return area_package


def main(cfg):
    """Run your preferred diagnostic using preprocessed files."""
    logger.info('-----------------------------------------------')
    logger.info('Running ESMValTool in personal diagnostic mode')
    logger.info('-----------------------------------------------')
    to_run = _import_my_diag(cfg['myDiagPlace'], cfg['myDiag'])
    for analysis_function in to_run.analyses:
        logger.info('Running function: %s', str(analysis_function))
        analysis_function(cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
