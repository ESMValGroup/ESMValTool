"""Barebone example of executing an interactive diagnostic."""
import os

import yaml
from esmvaltool.diag_scripts.shared import run_diagnostic_interactive
from IPython import embed


def main(cfg):
    """
    Function to jump in ipython.

    This basic example loads the runtime environment
    in the iPython console and reads in the preprcessor
    configuration settings (in a human-readable manner)
    in a nested dictionary ALL_SETTINGS.
    """
    parsed_config = os.path.join(cfg['run_dir'], 'parsed_settings.yml')
    with open(parsed_config, 'r') as file:
        ALL_SETTINGS = yaml.safe_load(file)
    embed()


if __name__ == '__main__':

    with run_diagnostic_interactive() as config:
        main(config)
