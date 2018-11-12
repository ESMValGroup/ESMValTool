"""
Run the first communication between esmvaltool's recipe and mip_convert.

This script sets up the correct rose suite directories to run mip_convert
on different UM suite data.
"""
import argparse
import logging
import os
import sys
import shutil
import subprocess
from distutils.version import LooseVersion
# configparser has changed names in python 3.x
if LooseVersion(sys.version) < LooseVersion("3.0"):
    import ConfigParser
else:
    import configparser as ConfigParser
import yaml  # noqa


# set up logging
logger = logging.getLogger(__name__)


# print the header
HEADER = r"""
______________________________________________________________________

     ESMValTool + mip_convert: linking mip_convert to ESMValTool
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
        nargs='+',
        help='Recipe files (list or single file)')
    parser.add_argument(
        '-m',
        '--mode',
        choices=['setup-only', 'setup-run-suites', 'postproc'],
        help='How to run: setup: sets up mipconvert suites only;\n' +
        'or setup-run-suites: sets up suites and runs them as well;\n' +
        'or postproc: grab the output from mip_convert and use it.')
    parser.add_argument(
        '-l',
        '--log-level',
        default='info',
        choices=['debug', 'info', 'warning', 'error'])
    args = parser.parse_args()
    return args


def _set_logger(logging, out_dir, log_file, log_level):
    # set logging for screen and file output
    root_logger = logging.getLogger()
    out_fmt = "%(asctime)s %(levelname)-8s %(name)s,%(lineno)s\t%(message)s"
    logging.basicConfig(
        filename=os.path.join(out_dir, log_file),
        filemode='a',
        format=out_fmt,
        datefmt='%H:%M:%S',
        level=logging.DEBUG)
    root_logger.setLevel(log_level.upper())
    logfmt = logging.Formatter(out_fmt)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logfmt)
    root_logger.addHandler(console_handler)


def read_yaml_file(yaml_file):
    """Read recipe into a dictionary."""
    with open(yaml_file, 'r') as yfile:
        loaded_file = yaml.safe_load(yfile)
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


def write_rose_conf(rose_config_template, recipe_file,
                    config_file, log_level):
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

    # set the logger to start outputting
    _set_logger(logging, conf_file['ROSES_OUTPUT'],
                'rose_suites_setup.log', log_level)
    logger.info(HEADER)

    # store the rose suite locations
    rose_suite_locations = []

    # loop through datasets (different suites for different datasets)
    for dataset in datasets:

        # set correct paths
        rose_suite = os.path.join(
            conf_file['ROSES_ROOT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']]
        )
        rose_suite_locations.append(rose_suite)
        rose_output = os.path.join(
            conf_file['ROSES_OUTPUT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']]
        )
        if os.path.exists(rose_suite):
            shutil.rmtree(rose_suite)
        shutil.copytree(conf_file['AK283'], rose_suite)
        if not os.path.exists(rose_output):
            os.makedirs(rose_output)

        # start logging
        logger.info("Working on dataset: %s", dataset)
        logger.info("Mapping dataset to suite: %s", rose_suite)
        logger.info("Output and logs written to: %s", rose_output)
        logger.info("Creating rose suite directories...")
        logger.info("Use rose-suite.conf template %s", rose_config_template)
        logger.info("Use user config file %s", config_file)

        # write the file
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
        with open(os.path.join(rose_suite, 'rose-suite.conf'), 'w') as r_c:
            logger.info("Writing rose-suite.conf file %s",
                        os.path.join(rose_suite, 'rose-suite.conf'))
            Config.write(r_c)

    return rose_suite_locations


def _put_in_env(env_script):
    """Put new system vars in environment."""
    logger.info("Setting environment for suite submission...")

    # First make it executable.
    chmod_command = ["chmod", "+x", env_script]
    proc = subprocess.Popen(chmod_command, stdout=subprocess.PIPE)
    proc.communicate()
    logger.info("Script %s is now executable.", env_script)

    # set the environment
    for line in open(env_script, 'r'):
        if line.split("=")[0] == 'export PATH':
            logger.info("Appending %s to path...",
                        line.split("=")[1].strip("\n"))
            add_path = line.split("=")[1].strip("\n").strip(":$PATH")
            os.environ["PATH"] += os.pathsep + add_path
        elif line.split("=")[0] == 'export PYTHONPATH':
            logger.info("Exporting %s as PYTHONPATH...",
                        line.split("=")[1].strip("\n"))
            os.environ["PYTHONPATH"] = line.split("=")[1].strip("\n")

    # print and check
    logger.info("New path: %s", str(os.environ["PATH"]))
    logger.info("mip_convert PYTHONPATH: %s", str(os.environ["PYTHONPATH"]))
    proc = subprocess.Popen(["which", "rose"], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    logger.info("rose: %s %s", out, err)
    proc = subprocess.Popen(["which", "mip_convert"], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    logger.info("mip_convert: %s %s", out, err)


def _source_envs(suite):
    """Source relevant environments."""
    # source the Met Office rose/cylc environment
    # and the suite specific environment
    suite_env = os.path.join(suite,
                             'env_setup_command_line.sh')  # suite env
    env_file_mo = os.path.join(suite, 'sourcepaths.sh')   # metomi env
    _put_in_env(suite_env)
    _put_in_env(env_file_mo)


def _run_suite(suite):
    """Run the mip_convert suite."""
    os.chdir(suite)
    logger.info("Submitting suite from %s", suite)
    proc = subprocess.Popen(["rose", "suite-run"], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    logger.info("Rose communications: %s %s", str(out), str(err))


def symlink_data(recipe_file, config_file, log_level):
    """Grab the mip_converted output and manage it for ESMValTool."""
    # get configuration and recipe
    recipe_object = read_yaml_file(recipe_file)
    conf_file = read_yaml_file(config_file)
    datasets = recipe_object['datasets']

    # create directory that stores all the output netCDF files
    sym_output_dir = os.path.join(conf_file['ROSES_OUTPUT'],
                                  'symlinks', recipe_file.strip('.yml'))
    if os.path.exists(sym_output_dir):
        shutil.rmtree(sym_output_dir)
    os.makedirs(sym_output_dir)

    # set the logger to start outputting
    _set_logger(logging, conf_file['ROSES_OUTPUT'],
                'file_simlink.log', log_level)
    logger.info(HEADER)

    # loop through all datasets to symlink output
    for dataset in datasets:
        rose_output = os.path.join(
            conf_file['ROSES_OUTPUT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']]
        )
        logger.info("Working on dataset: %s", dataset)
        logger.info("Output and logs written to: %s", rose_output)

        # loop through files
        for root, _, files in os.walk(rose_output):
            for xfile in files:
                real_file = os.path.join(root, xfile)
                imag_file = os.path.join(sym_output_dir, xfile)

                # symlink it if nc file
                if real_file.endswith('nc') and \
                        real_file.split('_')[2] == dataset['dataset']:
                    if not os.path.islink(imag_file):
                        logger.info("File to symlink: %s", real_file)
                        logger.info("Symlinked file: %s", imag_file)
                        os.symlink(real_file, imag_file)
                    else:
                        logger.info("Symlinked file exists...")
                        logger.info("Original file: %s", real_file)
                        logger.info("Symlinked file: %s", imag_file)


def main():
    """Run the the meat of the code."""
    args = get_args()
    rose_config_template = os.path.join(os.path.dirname(__file__),
                                        "rose-suite-template.conf")
    recipe_files = args.recipe_files
    config_file = args.config_file
    log_level = args.log_level
    for recipe_file in recipe_files:
        if args.mode == 'setup-only':
            # set up the rose suites
            write_rose_conf(rose_config_template, recipe_file,
                            config_file, log_level)
        elif args.mode == 'setup-run-suites':
            # setup roses
            roses = write_rose_conf(rose_config_template, recipe_file,
                                    config_file, log_level)
            # set up the environment and submit
            for rose in roses:
                _source_envs(rose)
                _run_suite(rose)
        elif args.mode == 'postproc':
            symlink_data(recipe_file, config_file, log_level)


if __name__ == '__main__':
    main()
