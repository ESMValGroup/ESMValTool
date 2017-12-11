"""Data finder module for the ESMValTool."""
# Authors:
# Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
# Mattia Righi (DLR, Germany - mattia.righi@dlr.de)

import logging
import os
import re
import subprocess
from datetime import datetime

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


def get_input_filename(variable, rootpath, drs):
    """Return the expected path to input file.

    This function should match the function get_input_filelist below.
    """
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

    dirname = os.path.join(dir1, dir2)

    # Simulate a latest version if required
    if '[latestversion]' in dirname:
        part1, part2 = dirname.split('[latestversion]')
        dirname = os.path.join(part1, 'dummy', part2)

    # Set the filename
    filename = replace_tags(cfg['input_file'], variable)
    if filename.endswith('*'):
        filename = filename.rstrip(
            '*') + "{start_year}01-{end_year}12.nc".format(**variable)

    # Full path to files
    return os.path.join(dirname, filename)


def get_input_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
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

    if not os.path.isdir(dir1):
        raise OSError('directory not found', dir1)

    # Set the drs
    _drs = drs.get(project, 'default')

    if _drs in cfg['input_dir']:
        dir2 = replace_tags(cfg['input_dir'][_drs], variable)
    else:
        raise KeyError(
            'drs %s for %s project not specified in config-developer file' %
            (_drs, project))

    dirname = os.path.join(dir1, dir2)

    # Find latest version if required
    if '[latestversion]' in dirname:
        part1 = dirname.split('[latestversion]')[0]
        part2 = dirname.split('[latestversion]')[1]
        list_versions = os.listdir(part1)
        list_versions.sort()
        latest = os.path.basename(list_versions[-1])
        dirname = os.path.join(part1, latest, part2)

    if not os.path.isdir(dirname):
        raise OSError('directory not found', dirname)

    # Set the filename
    filename = replace_tags(cfg['input_file'], variable)

    # Full path to files
    files = veto_files(variable, dirname, filename)

    return files


def get_output_file(variable, preproc_dir):
    """Return the full path to the output (preprocessed) file"""
    cfg = read_config_file(variable['project'])

    outfile = os.path.join(preproc_dir,
                           variable['project'],
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


def veto_files(variable, dirname, filename):
    """Find files between start_year and end_year.

    Function that does direct parsing of available datasource files
    and establishes if files are the needed ones or not.
    """
    files = find_files(dirname, filename)

    yr1 = int(variable['start_year'])
    yr2 = int(variable['end_year'])
    files = [s for s in files if time_check(s, yr1, yr2)]

    return files


# TODO: simplify: the functions below are way to complicated for what they do
def time_handling(year1, year1_model, year2, year2_model):
    """

    This function is responsible for finding the correct
    files for the needed timespan:

    year1 - the start year in files
    year1_model - the needed start year of data
    year2 - the last year in files
    year2_model - the needed last year of data
    WARNINGS:
    we reduce our analysis only to years

    """
    # model interval < data interval / file
    # model requirements completely within data stretch
    if year1 <= int(year1_model) and year2 >= int(year2_model):
        return True
    # model interval > data interval / file
    # data stretch completely within model requirements
    elif year1 >= int(year1_model) and year2 <= int(year2_model):
        return True
    # left/right overlaps and complete misses
    elif year1 <= int(year1_model) and year2 <= int(year2_model):
        # data is entirely before model
        if year2 < int(year1_model):
            return False
        # edge on
        elif year2 == int(year1_model):
            return True
        # data overlaps to the left
        elif year2 > int(year1_model):
            return True
    elif year1 >= int(year1_model) and year2 >= int(year2_model):
        # data is entirely after model
        if year1 > int(year2_model):
            return False
        # edge on
        elif year1 == int(year2_model):
            return True
        # data overlaps to the right
        elif year1 < int(year2_model):
            return True


# ---- function to handle various date formats
def date_handling(time1, time2):
    """
    This function deals with different input date formats e.g.
    time1 = 198204 or
    time1 = 19820422 or
    time1 = 198204220511 etc
    More formats can be coded in at this stage.
    Returns year 1 and year 2
    """
    # yyyymm
    if len(list(time1)) == 6 and len(list(time2)) == 6:
        y1 = datetime.strptime(time1, '%Y%m')
        year1 = y1.year
        y2 = datetime.strptime(time2, '%Y%m')
        year2 = y2.year
    else:
        # yyyymmdd
        if len(list(time1)) == 8 and len(list(time2)) == 8:
            y1 = datetime.strptime(time1, '%Y%m%d')
            year1 = y1.year
            y2 = datetime.strptime(time2, '%Y%m%d')
            year2 = y2.year
        # yyyymmddHHMM
        if len(list(time1)) == 12 and len(list(time2)) == 12:
            y1 = datetime.strptime(time1, '%Y%m%d%H%M')
            year1 = y1.year
            y2 = datetime.strptime(time2, '%Y%m%d%H%M')
            year2 = y2.year
    return year1, year2


# ---- function that does time checking on a file
def time_check(fpath, yr1, yr2):
    """
    fpath: full path to file
    yr1, yr2: model['start_year'], model['end_year']
    """
    ssp = fpath.split('/')
    av = ssp[-1]
    time_range = av.split('_')[-1].strip('.nc')
    time1 = time_range.split('-')[0]
    time2 = time_range.split('-')[1]
    year1 = date_handling(time1, time2)[0]
    year2 = date_handling(time1, time2)[1]
    return time_handling(year1, yr1, year2, yr2)
