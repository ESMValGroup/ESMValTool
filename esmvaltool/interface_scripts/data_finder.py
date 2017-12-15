"""Data finder module for the ESMValTool."""
# Authors:
# Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
# Mattia Righi (DLR, Germany - mattia.righi@dlr.de)

import logging
import os
import re
import subprocess

import yaml

logger = logging.getLogger(__name__)


def _read_config_file(cfg_file=None):
    """Parse the developer's configuration file."""
    if cfg_file is None:
        cfg_file = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'config-developer.yml',
        )

    with open(cfg_file, 'r') as file:
        cfg = yaml.safe_load(file)

    return cfg


_CFG = _read_config_file()


def cmip5_model2inst(model):
    """Return the institute given the model name in CMIP5."""
    logger.debug("Retrieving institute for CMIP5 model %s", model)
    return _CFG['CMIP5']['institute'][model]


def cmip5_mip2realm_freq(mip):
    """Return realm and frequency given the mip in CMIP5."""
    logger.debug("Retrieving realm and frequency for CMIP5 mip %s", mip)
    return _CFG['CMIP5']['realm_frequency'][mip]


def replace_tags(path, variable):
    """Replace tags in the config-developer's file with actual values."""
    path = path.strip('/')

    tlist = re.findall(r'\[([^]]*)\]', path)

    for tag in tlist:

        if tag == 'var':
            replacewith = variable['short_name']
        elif tag == 'field':
            replacewith = variable['field']
        elif tag in ('institute', 'freq', 'realm'):
            if tag in variable:
                replacewith = str(variable[tag])
            else:
                if tag == 'institute':
                    replacewith = cmip5_model2inst(variable['model'])
                elif tag == 'freq':
                    replacewith = cmip5_mip2realm_freq(variable['mip'])[1]
                elif tag == 'realm':
                    replacewith = cmip5_mip2realm_freq(variable['mip'])[0]
        elif tag == 'latestversion':  # handled separately later
            pass
        elif tag == 'tier':
            replacewith = ''.join(('Tier', str(variable['tier'])))
        elif tag == 'model':
            replacewith = variable['model']
        else:  # all other cases use the corresponding model dictionary key
            if tag in variable:
                replacewith = str(variable[tag])
            else:
                raise KeyError(
                    "Model key {} must be specified for project {}, check "
                    "your namelist entry".format(tag, variable['project']))

        path = path.replace('[' + tag + ']', replacewith)

    return path


def read_config_file(project, cfg_file=None):
    """Get developer-configuration for project."""
    logger.debug("Reading specifications for project %s from "
                 "config-developer file", project)

    if cfg_file is None:
        cfg = _CFG
    else:
        cfg = read_config_file(project, cfg_file)

    return cfg[project]


def get_input_dirname_template(variable, rootpath, drs):
    """Return a template of the full path to input directory."""
    project = variable['project']

    cfg = read_config_file(project)

    # Set the rootpath
    if project in rootpath:
        dir1 = rootpath[project]
    elif 'default' in rootpath:
        dir1 = rootpath['default']
    else:
        raise KeyError(
            'default rootpath must be specified in config-user file')

    # Set the drs
    _drs = drs.get(project, 'default')

    if _drs in cfg['input_dir']:
        dir2 = replace_tags(cfg['input_dir'][_drs], variable)
    else:
        raise KeyError(
            'drs {} for {} project not specified in config-developer file'
            .format(_drs, project))

    dirname_template = os.path.join(dir1, dir2)

    return dirname_template


def get_input_filename(variable, rootpath, drs):
    """Simulate a path to input file.

    This function should match the function get_input_filelist below.
    """
    dirname_template = get_input_dirname_template(variable, rootpath, drs)
    # Simulate a latest version if required
    if '[latestversion]' in dirname_template:
        part1, part2 = dirname_template.split('[latestversion]')
        dirname = os.path.join(part1, 'latestversion', part2)
    else:
        dirname = dirname_template

    # Set the filename
    cfg = read_config_file(variable['project'])
    filename = replace_tags(cfg['input_file'], variable)
    if filename.endswith('*'):
        filename = filename.rstrip(
            '*') + "{start_year}01-{end_year}12.nc".format(**variable)

    # Full path to files
    return os.path.join(dirname, filename)


def get_input_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
    dirname_template = get_input_dirname_template(variable, rootpath, drs)

    def check_isdir(path):
        """Raise an OSError if path is not a directory."""
        path = os.path.abspath(path)
        if not os.path.isdir(path):
            raise OSError('Directory not found: {}'.format(path))

    check_isdir(os.path.dirname(dirname_template))

    # Find latest version if required
    if '[latestversion]' in dirname_template:
        part1, part2 = dirname_template.split('[latestversion]')
        list_versions = os.listdir(part1)
        list_versions.sort()
        latest = os.path.basename(list_versions[-1])
        dirname = os.path.join(part1, latest, part2)
    else:
        dirname = dirname_template

    check_isdir(dirname)

    # Set the filename glob
    cfg = read_config_file(variable['project'])
    filename_glob = replace_tags(cfg['input_file'], variable)

    # Find files
    files = find_files(dirname, filename_glob)

    # Select files within the required time interval
    files = select_files(files, variable['start_year'], variable['end_year'])

    return files


def get_output_file(variable, preproc_dir):
    """Return the full path to the output (preprocessed) file"""
    cfg = read_config_file(variable['project'])

    outfile = os.path.join(preproc_dir, variable['preprocessor'],
                           replace_tags(cfg['output_file'], variable))
    outfile = ''.join((outfile, '.nc'))

    return outfile


def find_files(dirname, filename):
    """Find files

    Function that performs local search for files using `find'
    The depth is as high as possible so that find is fast.
    """
    flist = []

    # work only with existing dirs or allowed permission dirs
    strfindic = 'find {dirname} -follow -type f -iname "*{filename}*"'.format(
        dirname=dirname, filename=filename)
    logger.debug("Running %s", strfindic)
    proc = subprocess.Popen(
        strfindic, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    out, err = proc.communicate()
    if err:
        logger.warning("'%s' says:\n%s", strfindic, err)
    logger.debug("Result:\n%s", out)
    for line in out.split('\n')[0:-1]:
        flist.append(line)
    return flist


def get_start_end_year(filename):
    """Get the start and end year from a file name.

    This works for filenames matching *_YYYY*-YYYY*.*
    """
    name = os.path.splitext(filename)[0]
    start, end = name.split('_')[-1].split('-')
    start_year, end_year = int(start[:4]), int(end[:4])
    return start_year, end_year


def select_files(filenames, start_year, end_year):
    """Select files containing data between start_year and end_year.

    This works for filenames matching *_YYYY*-YYYY*.*
    """
    selection = []
    for filename in filenames:
        start, end = get_start_end_year(filename)
        if start <= end_year and end >= start_year:
            selection.append(filename)
    return selection
