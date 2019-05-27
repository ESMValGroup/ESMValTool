"""Data finder module for the ESMValTool."""
# Authors:
# Bouwe Andela (eScience, NL - b.andela@esciencecenter.nl)
# Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
# Mattia Righi (DLR, Germany - mattia.righi@dlr.de)

import fnmatch
import logging
import os
import re

import six

from ._config import get_project_config, replace_mip_fx
from .cmor.table import CMOR_TABLES

logger = logging.getLogger(__name__)


def find_files(dirnames, filenames):
    """Find files matching filenames in dirnames."""
    logger.debug("Looking for files matching %s in %s", filenames, dirnames)

    result = []
    for dirname in dirnames:
        for path, _, files in os.walk(dirname, followlinks=True):
            for filename in filenames:
                matches = fnmatch.filter(files, filename)
                result.extend(os.path.join(path, f) for f in matches)

    return result


def get_start_end_year(filename):
    """Get the start and end year from a file name.

    This works for filenames matching

    *[-,_]YYYY*[-,_]YYYY*.*
      or
    *[-,_]YYYY*.*
      or
    YYYY*[-,_]*.*
      or
    YYYY*[-,_]YYYY*[-,_]*.*
      or
    YYYY*[-,_]*[-,_]YYYY*.* (Does this make sense? Is this worth catching?)
    """
    name = os.path.splitext(filename)[0]

    filename = name.split(os.sep)[-1]
    filename_list = [elem.split('-') for elem in filename.split('_')]
    filename_list = [elem for sublist in filename_list for elem in sublist]

    pos_ydates = [elem.isdigit() and len(elem) >= 4 for elem in filename_list]
    pos_ydates_l = list(pos_ydates)
    pos_ydates_r = list(pos_ydates)

    for ind, _ in enumerate(pos_ydates_l):
        if ind != 0:
            pos_ydates_l[ind] = (pos_ydates_l[ind - 1] and pos_ydates_l[ind])

    for ind, _ in enumerate(pos_ydates_r):
        if ind != 0:
            pos_ydates_r[-ind - 1] = (pos_ydates_r[-ind]
                                      and pos_ydates_r[-ind - 1])

    dates = [
        filename_list[ind] for ind, _ in enumerate(pos_ydates)
        if pos_ydates_r[ind] or pos_ydates_l[ind]
    ]

    if len(dates) == 1:
        start_year = int(dates[0][:4])
        end_year = start_year
    elif len(dates) == 2:
        start_year, end_year = int(dates[0][:4]), int(dates[1][:4])
    else:
        raise ValueError('Name {0} dates do not match a recognized '
                         'pattern'.format(name))

    return start_year, end_year


def select_files(filenames, start_year, end_year):
    """Select files containing data between start_year and end_year.

    This works for filenames matching *_YYYY*-YYYY*.* or *_YYYY*.*
    """
    selection = []
    for filename in filenames:
        start, end = get_start_end_year(filename)
        if start <= end_year and end >= start_year:
            selection.append(filename)
    return selection


def _replace_tags(path, variable, fx_var=None):
    """Replace tags in the config-developer's file with actual values."""
    path = path.strip('/')

    tlist = re.findall(r'\[([^]]*)\]', path)

    paths = [path]
    for tag in tlist:
        original_tag = tag
        tag, _, _ = _get_caps_options(tag)

        if tag == 'fx_var':
            replacewith = fx_var
        elif tag == 'latestversion':  # handled separately later
            continue
        elif tag in variable:
            replacewith = variable[tag]
        else:
            raise KeyError("Dataset key {} must be specified for {}, check "
                           "your recipe entry".format(tag, variable))

        paths = _replace_tag(paths, original_tag, replacewith)

    return paths


def _replace_tag(paths, tag, replacewith):
    """Replace tag by replacewith in paths."""
    _, lower, upper = _get_caps_options(tag)
    result = []
    if isinstance(replacewith, (list, tuple)):
        for item in replacewith:
            result.extend(_replace_tag(paths, tag, item))
    else:
        text = _apply_caps(str(replacewith), lower, upper)
        result.extend(p.replace('[' + tag + ']', text) for p in paths)
    return result


def _get_caps_options(tag):
    lower = False
    upper = False
    if tag.endswith('.lower'):
        lower = True
        tag = tag[0:-6]
    elif tag.endswith('.upper'):
        upper = True
        tag = tag[0:-6]
    return tag, lower, upper


def _apply_caps(original, lower, upper):
    if lower:
        return original.lower()
    if upper:
        return original.upper()
    return original


def _resolve_latestversion(dirname_template):
    """Resolve the 'latestversion' tag."""
    if '[latestversion]' not in dirname_template:
        return dirname_template

    # Find latest version
    part1, part2 = dirname_template.split('[latestversion]')
    part2 = part2.lstrip(os.sep)
    if os.path.exists(part1):
        versions = os.listdir(part1)
        versions.sort(reverse=True)
        for version in ['latest'] + versions:
            dirname = os.path.join(part1, version, part2)
            if os.path.isdir(dirname):
                return dirname

    return dirname_template


def _select_drs(input_type, drs, project):
    """Select the directory structure of input path."""
    cfg = get_project_config(project)
    input_path = cfg[input_type]
    if isinstance(input_path, six.string_types):
        return input_path

    structure = drs.get(project, 'default')
    if structure in input_path:
        return input_path[structure]

    raise KeyError(
        'drs {} for {} project not specified in config-developer file'.format(
            structure, project))


def get_rootpath(rootpath, project):
    """Select the rootpath."""
    if project in rootpath:
        return rootpath[project]
    if 'default' in rootpath:
        return rootpath['default']
    raise KeyError('default rootpath must be specified in config-user file')


def _find_input_dirs(variable, rootpath, drs, fx_var=None):
    """Return a the full paths to input directories."""
    project = variable['project']

    root = get_rootpath(rootpath, project)
    input_type = 'input_{}dir'.format('fx_' if fx_var else '')
    path_template = _select_drs(input_type, drs, project)

    dirnames = []
    for dirname_template in _replace_tags(path_template, variable, fx_var):
        for base_path in root:
            dirname = os.path.join(base_path, dirname_template)
            dirname = _resolve_latestversion(dirname)
            if os.path.exists(dirname):
                logger.debug("Found %s", dirname)
                dirnames.append(dirname)
            else:
                logger.debug("Skipping non-existent %s", dirname)

    return dirnames


def _get_filenames_glob(variable, drs, fx_var=None):
    """Return patterns that can be used to look for input files."""
    input_type = 'input_{}file'.format('fx_' if fx_var else '')
    path_template = _select_drs(input_type, drs, variable['project'])
    filenames_glob = _replace_tags(path_template, variable, fx_var)
    return filenames_glob


def _find_input_files(variable, rootpath, drs, fx_var=None):
    logger.debug("Looking for input %sfiles for variable %s of dataset %s",
                 fx_var + ' fx ' if fx_var else '', variable['short_name'],
                 variable['dataset'])

    input_dirs = _find_input_dirs(variable, rootpath, drs, fx_var)
    filenames_glob = _get_filenames_glob(variable, drs, fx_var)
    files = find_files(input_dirs, filenames_glob)

    return files


def get_input_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
    files = _find_input_files(variable, rootpath, drs)
    files = select_files(files, variable['start_year'], variable['end_year'])
    return files


def get_input_fx_filelist(variable, rootpath, drs):
    """Return a dict with the full path to fx input files."""
    fx_files = {}
    for fx_var in variable['fx_files']:
        var = dict(variable)
        var['mip'] = replace_mip_fx(fx_var)
        table = CMOR_TABLES[var['cmor_table']].get_table(var['mip'])
        var['frequency'] = table.frequency
        realm = getattr(table.get(var['short_name']), 'modeling_realm', None)
        var['modeling_realm'] = realm if realm else table.realm

        files = _find_input_files(var, rootpath, drs, fx_var)
        fx_files[fx_var] = files[0] if files else None

    return fx_files


def get_output_file(variable, preproc_dir):
    """Return the full path to the output (preprocessed) file."""
    cfg = get_project_config(variable['project'])

    # Join different experiment names
    if isinstance(variable.get('exp'), (list, tuple)):
        variable = dict(variable)
        variable['exp'] = '-'.join(variable['exp'])

    outfile = os.path.join(
        preproc_dir,
        variable['diagnostic'],
        variable['variable_group'],
        _replace_tags(cfg['output_file'], variable)[0] + '.nc',
    )

    return outfile


def get_statistic_output_file(variable, preproc_dir):
    """Get multi model statistic filename depending on settings."""
    template = os.path.join(
        preproc_dir,
        '{diagnostic}',
        '{variable_group}',
        '{dataset}_{mip}_{short_name}_{start_year}-{end_year}.nc',
    )

    outfile = template.format(**variable)

    return outfile
