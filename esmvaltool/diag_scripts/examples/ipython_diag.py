"""Barebone example of executing an interactive diagnostic."""
import logging
import os

# import subprocess
from IPython import embed
import yaml
from esmvaltool.diag_scripts.shared import run_diagnostic_interactive


logger = logging.getLogger(os.path.basename(__file__))


def parse_settings(settings_file):
    """
    Read and output the contents of the settings file.

    Returns a nested dictionary that has the shape:
    DIAG_NAMES: VARS: DATASETS: DATSET_PROPERTIES eg
    {'validation_basic': {'tas': {'MPI-ESM-LR':{}}}}
    """
    parser_meta = {}
    for metadata_file in settings_file['input_files']:
        diag = metadata_file.split(os.sep)[-3]
        parser_meta[diag] = {}
        var = metadata_file.split(os.sep)[-2]
        parser_meta[diag][var] = {}
        logger.debug("Diagnostic %s for variable %s", diag, var)
        with open(metadata_file, 'r') as file:
            meta_dict = yaml.safe_load(file)
        for key in meta_dict:
            dataset = meta_dict[key]['dataset']
            parser_meta[diag][var][dataset] = {}
            for sub_key in meta_dict[key]:
                parser_meta[diag][var][dataset][sub_key] = \
                    meta_dict[key][sub_key]
    logger.debug("****************************************")
    logger.debug("Global configuration dictionary: %s ", parser_meta)
    logger.debug("****************************************")
    return parser_meta


def main(cfg):
    """Function to create a parsed_settings file."""
    parsed_settings = parse_settings(cfg)
    parsed_settings_file = os.path.join(cfg['run_dir'],
                                        'parsed_settings.yml')
    with open(parsed_settings_file, 'w') as file:
        yaml.safe_dump(parsed_settings, file)
    parsed_config = os.path.join(cfg['run_dir'], 'parsed_settings.yml')
    return parsed_config


if __name__ == '__main__':

    with run_diagnostic_interactive() as config:
        # launching ipython
        env = dict(os.environ)
        env['parsed_settings'] = main(config)
        logger.info("Running interactive ipython")

    # This is waiting for process to communicate to task
    # since task.run(stdout=subprocess.PIPE) but ipython
    # desnt work with PIPE
    # process = subprocess.Popen(['ipython'], env=env)
    # process.communicate()

# task.run still waiting and not launching ipython console
embed()
