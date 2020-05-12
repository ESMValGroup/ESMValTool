"""
Run the first communication between esmvaltool's recipe and mip_convert.

Description:
------------

This script sets up the correct rose suite directories to run mip_convert
on different UM suite data. You can run this tool in three different ways:
 - (with -m --mode option) setup-only: will set up the mip convert rose
   directories only; it will use the -c configuration file for user options;
 - (with -m --mode option) setup-run-suites: will set up the mip convert rose
   suites and will go ahead and submit them to cylc via rose suite-run;
 - (with -m --mode option) postproc: will symlink newly created netCDF data
   into a directory per esmvaltool recipe; note that for now, there is no
   DRS-like path set up in that directory;

Usage:
------
-c --config-file:  [REQUIRED]  user specific configuration file;
-r --recipe-file:  [REQUIRED]  single or multiple (space-sep) recipe files;
-m --mode:         [OPTIONAL]  running mode (setup-only, setup-run-suites,
                               postproc), default=setup-only
-l --log-level:    [OPTIONAL]  log level, default=info

Environment
-----------
current JASMIN rose/cyclc need python2.7; esmvaltool needs python3.x
So it is impossible at the moment to run this script as executable from an
esmvaltool environment. Instead, you can run it as a stand-alone tool in a
python 2.7 environment, intwo stages:

[set up mip_convert suites and run them]
python esmvt_mipconv_setup.py -c config.yml -r recipe.yml -m setup-run-suites
[check succesful completion of mip_convert suites]
[run the symlinking]
python esmvt_mipconv_setup.py -c config.yml -r recipe.yml -m postproc

A practical example of running the tool can be found on JASMIN:
/home/users/valeriu/esmvaltool_mip_convert
There you will find the two component shells: run_conversion
and run_symlink, as well as an example how to set the configuration file.

The suite used is now on MOSRS (as of 3 December 2018): u-bd681
You can use the default location on Jasmin:
DEFAULT_SUITE_LOCATION = "/home/users/valeriu/roses/u-bd681"
alternatively this can be turned off, should you want to check out the suite
off MOSRS and use it locally.

Contact:
--------
author: Valeriu Predoi (UREAD, valeriu.predoi@ncas.ac.uk)
"""
import argparse
import datetime
import logging
import os
import sys
import shutil
import subprocess
import socket
from distutils.version import LooseVersion
# configparser has changed names in python 3.x
if LooseVersion(sys.version) < LooseVersion("3.0"):
    import ConfigParser
else:
    import configparser as ConfigParser
import yaml  # noqa

####################
# global variables #
####################

# the tool uses a specially tailored mip_convert Rose suite
# locations of the suite depends on the host
host_name = socket.gethostname().split('.')
if len(host_name) > 1:
    if host_name[1] == 'ceda':
        # default location for mip_convert suite on JASMIN:
        # previous suite: u-ak283_esmvt; new one u-bd681
        # DEFAULT_SUITE_LOCATION = "/home/users/valeriu/roses/u-ak283_esmvt"
        DEFAULT_SUITE_LOCATION = "/home/users/valeriu/roses/u-bd681"
        # note that you can fcm checkout it straight from the MOSRS

# stream mapping; taken from hadsdk.streams
# these are used to set defaults if not overrides
STREAM_MAP = {
    'CMIP5': {
        '3hr': 'apk',
        '6hrPlev': 'apc',
        '6hrlev': 'apg',
        'Amon': 'apm',
        'Lmon': 'apm',
        'LImon': 'apm',
        'Oday': 'opa',
        'Omon': 'opm',
        'Oyr': 'opy',
        'CF3hr': 'apk',
        'CFday': 'apa',
        'CFmon': 'apm',
        'CFsubhr': 'ape',
        'day': 'apa'
    },
    'CMIP6': {
        '3hr': 'ap8',
        '6hrLev': 'ap7',
        '6hrPlev': 'ap7',
        '6hrPlevPt': 'ap7',
        'AERday': 'ap6',
        'AERhr': 'ap9',
        'AERmon': 'ap4',
        'AERmonZ': 'ap4',
        'Amon': 'ap5',
        'CF3hr': 'ap8',
        'CFday': 'ap6',
        'CFmon': 'ap5',
        'E1hr': 'ap9',
        'E1hrClimMon': 'ap9',
        'E3hr': 'ap8',
        'E3hrPt': 'ap8',
        'E6hrZ': 'ap7',
        'Eday': 'ap6',
        'EdayZ': 'ap6',
        'Efx': 'ancil',
        'Emon': 'ap5',
        'EmonZ': 'ap5',
        'Esubhr': 'ap8',
        'Eyr': 'ap5',
        'LImon': 'ap5',
        'Lmon': 'ap5',
        'Oday': 'ond',
        'Ofx': 'ancil',
        'Omon': 'onm',
        'SIday': 'ind',
        'SImon': 'inm',
        'day': 'ap6',
        'fx': 'ancil',
        'prim1hrpt': 'ap9',
        'prim3hr': 'ap8',
        'prim3hrpt': 'ap8',
        'prim6hr': 'ap7',
        'prim6hrpt': 'ap7',
        'primDay': 'ap6',
        'primMon': 'ap5',
        'primSIday': 'ap6'
    }
}

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
        default='setup-only',
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


def write_rose_conf(rose_config_template, recipe_file, config_file, log_level):
    """Write the new rose conf file per suite."""
    # Build the ConfigParser object
    config = ConfigParser.ConfigParser()
    config.optionxform = str
    config.read(rose_config_template)
    recipe_object = read_yaml_file(recipe_file)
    conf_file = read_yaml_file(config_file)
    datasets = recipe_object['datasets']

    # check if dataset needs analysis
    datasets_to_analyze = []
    for dataset in datasets:
        if dataset['dataset'] not in conf_file['DATASET_TO_SUITE']:
            logger.warning("Dataset %s has no mapping to suite",
                           dataset['dataset'])
            logger.warning("Assuming data retrival from elsewhere.")
        else:
            datasets_to_analyze.append(dataset)
    diagnostics = recipe_object['diagnostics']
    active_streams = map_var_to_stream(diagnostics, conf_file['STREAM_MAP'])

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
    if not os.path.exists(conf_file['ROSES_OUTPUT']):
        os.makedirs(conf_file['ROSES_OUTPUT'])
    _set_logger(logging, conf_file['ROSES_OUTPUT'], 'rose_suites_setup.log',
                log_level)
    logger.info(HEADER)

    # store the rose suite locations
    rose_suite_locations = []

    # loop through datasets (different suites for different datasets)
    for dataset in datasets_to_analyze:

        # set correct paths
        rose_suite = os.path.join(
            conf_file['ROSES_ROOT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']])
        rose_suite_locations.append(rose_suite)
        rose_output = os.path.join(
            conf_file['ROSES_OUTPUT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']])
        if os.path.exists(rose_suite):
            shutil.rmtree(rose_suite)
        if os.path.exists(DEFAULT_SUITE_LOCATION):
            shutil.copytree(DEFAULT_SUITE_LOCATION, rose_suite)
        else:
            logger.error("Default Suite Location not found: %s",
                         DEFAULT_SUITE_LOCATION)
            break
        if not os.path.exists(rose_output):
            os.makedirs(rose_output)
        new_mipconv_config = os.path.join(rose_suite, 'mip_convert_config')

        # start logging
        logger.info("Working on dataset: %s", dataset)
        logger.info("Mapping dataset to suite: %s", rose_suite)
        logger.info("Output and logs written to: %s", rose_output)
        logger.info("Creating rose suite directories...")
        logger.info("Use rose-suite.conf template %s", rose_config_template)
        logger.info("Use user config file %s", config_file)

        # write the file
        config.set('jinja2:suite.rc', 'INPUT_DIR',
                   '"' + conf_file['INPUT_DIR'] + '"')
        config.set('jinja2:suite.rc', 'OUTPUT_DIR', '"' + rose_output + '"')
        config.set('jinja2:suite.rc', 'CDDS_DIR',
                   '"' + DEFAULT_SUITE_LOCATION + '"')
        config.set('jinja2:suite.rc', 'MIP_CONVERT_CONFIG_DIR',
                   '"' + new_mipconv_config + '"')
        config.set('jinja2:suite.rc', 'ACTIVE_STREAMS', str(active_streams))
        config.set('jinja2:suite.rc', 'STREAM_TIME_OVERRIDES',
                   str(stream_overrides))
        config.set('jinja2:suite.rc', 'FIRST_YEAR', str(dataset['start_year']))
        config.set('jinja2:suite.rc', 'REF_YEAR', str(dataset['start_year']))
        config.set('jinja2:suite.rc', 'FINAL_YEAR', str(dataset['end_year']))
        config.set('jinja2:suite.rc', 'STREAM_COMPONENTS',
                   str(stream_components))
        config.set('jinja2:suite.rc', 'CYCLING_FREQUENCIES',
                   str(cycling_frequencies))
        config.set(
            'jinja2:suite.rc', 'TARGET_SUITE_NAME',
            '"' + conf_file['DATASET_TO_SUITE'][dataset['dataset']] + '"')
        with open(os.path.join(rose_suite, 'rose-suite.conf'), 'w') as r_c:
            logger.info("Writing rose-suite.conf file %s",
                        os.path.join(rose_suite, 'rose-suite.conf'))
            config.write(r_c)

        # now that we have to conf file set up we need to
        # edit the mip_convert configuration file with the correct data
        for key, values in conf_file['STREAM_COMPONENTS'].items():
            for comp in values:
                mipconv_config = os.path.join(new_mipconv_config,
                                              'mip_convert.cfg.' + comp)
                _edit_mip_convert_config(mipconv_config, conf_file, dataset,
                                         key)

    return rose_suite_locations


def _edit_mip_convert_config(mipconv_config, conf_file, dataset, stream):
    """Edit the mip_convert file for correct runs."""
    # set the correct variables
    base_date = str(dataset['start_year']) + '-01-01-00-00-00'
    suite_id = conf_file['DATASET_TO_SUITE'][dataset['dataset']]
    cdds_dir = os.path.join(DEFAULT_SUITE_LOCATION, 'mip_convert_aux')

    # Build the ConfigParser object
    config = ConfigParser.ConfigParser()
    config.optionxform = str
    config.read(mipconv_config)

    # set the correct fields
    config.set('COMMON', 'cdds_dir', cdds_dir)
    config.set('request', 'base_date', base_date)
    config.set('request', 'suite_id', suite_id)
    stream_section = '_'.join(['stream', stream])
    # add the section if not there already
    if not config.has_section(stream_section):
        config.add_section(stream_section)
    if 'mip' not in dataset:
        # can work without any mip in dataset
        # will not take it from diagnostic (will assemble
        # all possible mappings instead)
        logger.warning("No mip in the recipe dataset section.")
        logger.warning("Assigning mapping from default dictionary.")
        stream_map_default = STREAM_MAP[dataset['project']]
        variables = []
        cmip_types = []
        for key, val in conf_file['STREAM_MAP'].items():
            for key_def, val_def in stream_map_default.items():
                if val == val_def:
                    cmip_types.append('_'.join([dataset['project'], key_def]))
                    variables.append(key)
        str_variables = ' '.join(list(set([v for v in variables])))
        if variables:
            for cmip_type in cmip_types:
                config.set(stream_section, cmip_type, str_variables)
    else:
        cmip_type = '_'.join([dataset['project'], dataset['mip']])
        all_vars = conf_file['STREAM_MAP'].keys()
        str_variables = ' '.join(
            [v for v in all_vars if conf_file['STREAM_MAP'][v] == stream])
        config.set(stream_section, cmip_type, str_variables)

    # write to file
    with open(mipconv_config, 'w') as r_c:
        logger.info("Writing mip_convert config file %s", mipconv_config)
        config.write(r_c)


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
    suite_env = os.path.join(suite, 'env_setup_command_line.sh')  # suite env
    env_file_mo = os.path.join(suite, 'sourcepaths.sh')  # metomi env
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
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join((recipe_file.strip('.yml'), now))
    sym_output_dir = os.path.join(conf_file['ROSES_OUTPUT'],
                                  'mip_convert_symlinks', new_subdir)
    if not os.path.exists(sym_output_dir):
        os.makedirs(sym_output_dir)

    # set the logger to start outputting
    _set_logger(logging, conf_file['ROSES_OUTPUT'], 'file_simlink.log',
                log_level)
    logger.info(HEADER)

    # loop through all datasets to symlink output
    for dataset in datasets:
        rose_output = os.path.join(
            conf_file['ROSES_OUTPUT'],
            conf_file['DATASET_TO_SUITE'][dataset['dataset']])
        logger.info("Working on dataset: %s", dataset)
        logger.info("Output and logs written to: %s", rose_output)

        # create the dataset dir
        dataset_output = os.path.join(sym_output_dir, dataset['dataset'])
        if os.path.exists(dataset_output):
            shutil.rmtree(dataset_output)
        os.makedirs(dataset_output)

        # loop through files
        for root, _, files in os.walk(rose_output):
            for xfile in files:
                real_file = os.path.join(root, xfile)
                imag_file = os.path.join(dataset_output, xfile)

                # symlink it if nc file
                if real_file.endswith('.nc') and \
                        xfile.split('_')[2] == dataset['dataset']:
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
    logger.info("Running main function...")
    args = get_args()
    rose_config_template = os.path.join(
        os.path.dirname(__file__), "rose-suite-template.conf")

    # make sure the file is retrieved nonetheless
    if not os.path.isfile(rose_config_template):
        logger.info("Fetching rose template config from suite %s",
                    DEFAULT_SUITE_LOCATION)
        rose_config_template = os.path.join(DEFAULT_SUITE_LOCATION,
                                            "rose-suite-template.conf")

    recipe_files = args.recipe_files
    config_file = args.config_file
    log_level = args.log_level
    for recipe_file in recipe_files:
        if args.mode == 'setup-only':
            # set up the rose suites
            write_rose_conf(rose_config_template, recipe_file, config_file,
                            log_level)
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
