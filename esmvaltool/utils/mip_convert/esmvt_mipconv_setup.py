"""
Run the first communication between esmvaltool's recipe and mip_convert.

This script sets up the correct rose suite directories to run mip_convert
on different UM suite data.
"""
import argparse
import os
import shutil
import configparser as ConfigParser
import yaml


# print the header
HEADER = r"""
______________________________________________________________________
          _____ ____  __  ____     __    _ _____           _
         | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
         |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
         | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
         |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
______________________________________________________________________

""" + __doc__


def get_args():
    """Define the `esmvaltool` command line."""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=HEADER,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-c',
        '--config-file',
        default=os.path.join(os.path.dirname(__file__), 'config-user.yml'),
        help='Configuration file')
    parser.add_argument(
        '-r',
        '--recipe-files',
        type=str,
        help='Recipe files (list or single file)')
    args = parser.parse_args()
    return args


def read_yaml_file(yaml_file):
    """Read recipe into a dictionary."""
    with open(yaml_file, 'r') as file:
        loaded_file = yaml.safe_load(file)
    return loaded_file


def map_var_to_stream(diagnostics, stream_map):
    """Map variable standard name to stream string."""
    stream_list = []
    for _, diag in diagnostics.items():
        for var in diag['variables']:
            stream = stream_map[var]
            stream_list.append(stream)
    stream_list = list(set(stream_list))
    return stream_list


def write_rose_conf(rose_config_template, recipe_file, config_file):
    """Write the new rose conf file per suite."""
    # Build the ConfigParser object
    Config = ConfigParser.ConfigParser()
    Config.optionxform = str
    Config.read(rose_config_template)
    recipe_object = read_yaml_file(recipe_file)
    conf_file = read_yaml_file(config_file)
    datasets = recipe_object['datasets']
    diagnostics = recipe_object['diagnostics']
    active_streams = map_var_to_stream(diagnostics,
                                       conf_file['STREAM_MAP'])

    # set stream overrides to None and set components
    # also set CYCLING_FREQUENCIES to P1Y overall
    stream_overrides = {}
    stream_components = {}
    cycling_frequencies = {}
    for stream in active_streams:
        stream_overrides[stream] = 'None'
        stream_components[stream] = conf_file['STREAM_COMPONENTS'][stream]
        cycling_frequencies[stream] = 'P1Y'

    # loop through datasets (different suites for different datasets)
    for dataset in datasets:

        # set correct paths
        rose_suite = os.path.join(
            conf_file['ROSES_ROOT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']]
        )
        rose_output = os.path.join(
            conf_file['ROSES_OUTPUT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']]
        )
        if os.path.exists(rose_suite):
            shutil.rmtree(rose_suite)
        shutil.copytree(conf_file['AK283'], rose_suite)
        if not os.path.exists(rose_output):
            os.makedirs(rose_output)

        # write the file
        rose_config = open(os.path.join(rose_suite, 'rose-suite.conf'), 'w')
        Config.set('jinja2:suite.rc', 'INPUT_DIR',
                   '"' + conf_file['INPUT_DIR'] + '"')
        Config.set('jinja2:suite.rc', 'OUTPUT_DIR', '"' + rose_output + '"')
        Config.set('jinja2:suite.rc', 'ACTIVE_STREAMS', str(active_streams))
        Config.set('jinja2:suite.rc', 'STREAM_TIME_OVERRIDES',
                   str(stream_overrides))
        Config.set('jinja2:suite.rc', 'FIRST_YEAR', str(dataset['start_year']))
        Config.set('jinja2:suite.rc', 'REF_YEAR', str(dataset['start_year']))
        Config.set('jinja2:suite.rc', 'FINAL_YEAR', str(dataset['end_year']))
        Config.set('jinja2:suite.rc', 'STREAM_COMPONENTS',
                   str(stream_components))
        Config.set('jinja2:suite.rc', 'CYCLING_FREQUENCIES',
                   str(cycling_frequencies))
        Config.set(
            'jinja2:suite.rc', 'TARGET_SUITE_NAME',
            '"' + conf_file['DATASET_TO_SUITE'][dataset['dataset']] + '"'
        )
        Config.write(rose_config)


def main():
    """Run the the meat of the code."""
    args = get_args()
    rose_config_template = "rose-suite-template.conf"
    recipe_files = args.recipe_files
    config_file = args.config_file
    for recipe_file in recipe_files.split(','):
        write_rose_conf(rose_config_template, recipe_file, config_file)


if __name__ == '__main__':
    main()
