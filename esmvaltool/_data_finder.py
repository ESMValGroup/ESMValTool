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

from ._config import (cmip5_dataset2inst, cmip5_mip2realm_freq,
                      get_project_config)

logger = logging.getLogger(__name__)


def find_files(dirname, filename):
    """Find files matching filename."""
    logger.debug("Looking for files matching %s in %s", filename, dirname)

    result = []
    for path, _, files in os.walk(dirname, followlinks=True):
        files = fnmatch.filter(files, filename)
        if files:
            result.extend(os.path.join(path, f) for f in files)

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


def replace_tags(path, variable, j=None, i=None):
    """Replace tags in the config-developer's file with actual values."""
    path = path.strip('/')

    tlist = re.findall(r'\[([^]]*)\]', path)

    for tag in tlist:
        original_tag = tag
        tag, lower, upper = _get_caps_options(tag)

        if tag == 'var':
            replacewith = variable['short_name']
        elif tag == 'fx_var':
            replacewith = variable['fx_files'][i]
        elif tag == 'field':
            replacewith = variable['field']
        elif tag in ('institute', 'freq', 'realm'):
            if tag in variable:
                replacewith = str(variable[tag])
            else:
                if tag == 'institute':
                    replacewith = cmip5_dataset2inst(variable['dataset'])
                elif tag == 'freq':
                    replacewith = cmip5_mip2realm_freq(variable['mip'])[1]
                elif tag == 'realm':
                    replacewith = cmip5_mip2realm_freq(variable['mip'])[0]
        elif tag == 'latestversion':  # handled separately later
            continue
        elif tag == 'tier':
            replacewith = ''.join(('Tier', str(variable['tier'])))
        elif tag == 'dataset':
            replacewith = variable['dataset']
        else:  # all other cases use the corresponding dataset dictionary key
            if tag in variable:
                replacewith = str(variable[tag])
            else:
                raise KeyError(
                    "Dataset key {} must be specified for project {}, check "
                    "your recipe entry".format(tag, variable['project']))

        if not isinstance(replacewith, list):
            path = path.replace('[' + original_tag + ']',
                                _apply_caps(replacewith, lower, upper))
        else:
            path = [
                path.replace('[' + original_tag + ']',
                             _apply_caps(dkrz_place, lower, upper))
                for dkrz_place in replacewith
            ][j]
    return path


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
    elif upper:
        return original.upper()
    return original


def get_input_dirname_template(variable, rootpath, drs):
    """Return a template of the full path to input directory."""
    project = variable['project']

    cfg = get_project_config(project)

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
    input_dir = cfg['input_dir']
    if isinstance(input_dir, six.string_types):
        dir2 = replace_tags(input_dir, variable)
    elif _drs in input_dir:
        try:
            insts = cmip5_dataset2inst(variable['dataset'])
        except KeyError as msg:
            logger.debug('CMIP5 dataset2inst: %s', msg)
            insts = 0
        dirs2 = []
        if isinstance(insts, list):
            for j in range(len(insts)):
                dir2 = replace_tags(input_dir[_drs], variable, j)
                dirs2.append(dir2)
        else:
            dir2 = replace_tags(input_dir[_drs], variable)
            dirs2.append(dir2)
    else:
        raise KeyError(
            'drs {} for {} project not specified in config-developer file'
            .format(_drs, project))

    dirname_template = [os.path.join(dir1, dir_2) for dir_2 in dirs2]

    return dirname_template


def get_input_fx_dirname_template(variable, rootpath, drs):
    """Return a template of the full path to input directory."""
    project = variable['project']

    cfg = get_project_config(project)

    dirs = []
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
    input_dir = cfg['fx_dir']

    # Loop through fx_files types
    for fx_ind in range(len(variable['fx_files'])):

        # Need to reassign the mip so we can find sftlf/of
        # make a copy of variable -> new_variable for this
        new_variable = dict(variable)
        if new_variable['fx_files'][fx_ind] == 'sftlf' or \
                new_variable['fx_files'][fx_ind] == 'areacella':
            new_variable['mip'] = 'Amon'
        elif new_variable['fx_files'][fx_ind] == 'sftof' or \
                new_variable['fx_files'][fx_ind] == 'areacello':
            new_variable['mip'] = 'Omon'

        if isinstance(input_dir, six.string_types):
            dir2 = replace_tags(input_dir, new_variable, i=fx_ind)
        elif _drs in input_dir:
            dir2 = replace_tags(input_dir[_drs], new_variable, i=fx_ind)
        else:
            raise KeyError(
                'drs {} for {} project not specified in config-developer file'
                .format(_drs, project))

        dirname_template = os.path.join(dir1, dir2)
        dirs.append(dirname_template)

    return dirs


def get_input_filename(variable, rootpath, drs):
    """Simulate a path to input file.

    This function should match the function get_input_filelist below.
    """
    dirname_templates = get_input_dirname_template(variable, rootpath, drs)
    for dirname_template in dirname_templates:
        # Simulate a latest version if required
        if '[latestversion]' in dirname_template:
            part1, part2 = dirname_template.split('[latestversion]')
            dirname = os.path.join(part1, 'latestversion', part2)
        else:
            dirname = dirname_template

        # Set the filename
        filename = _get_filename(variable, drs)
        if filename.endswith('*'):
            filename = filename.rstrip(
                '*') + "{start_year}01-{end_year}12.nc".format(**variable)

        # Full path to files
        return os.path.join(dirname, filename)


def get_input_fx_filename(variable, rootpath, drs):
    """Simulate a path to input file.

    This function should match the function get_input_filelist below.
    """
    files = []
    dirname_templates = get_input_fx_dirname_template(variable, rootpath, drs)
    for j, dirname_template in zip(
            range(len(dirname_templates)), dirname_templates):
        # Simulate a latest version if required
        if '[latestversion]' in dirname_template:
            part1, part2 = dirname_template.split('[latestversion]')
            dirname = os.path.join(part1, 'latestversion', part2)
        else:
            dirname = dirname_template

        # Set the filename
        filename = _get_fx_filename(variable, drs, j)

        # Full path to files
        files.append(os.path.join(dirname, filename))

    return files


def _get_filename(variable, drs):
    project = variable['project']
    cfg = get_project_config(project)

    input_file = cfg['input_file']
    _drs = drs.get(project, 'default')
    if not isinstance(input_file, six.string_types):
        if _drs in input_file:
            input_file = input_file[_drs]
        else:
            raise KeyError(
                'drs {} for {} project not specified for input_file '
                'in config-developer file'.format(_drs, project))
    filename = replace_tags(input_file, variable)
    return filename


def _get_fx_filename(variable, drs, j):
    project = variable['project']
    cfg = get_project_config(project)

    input_file = cfg['fx_file']
    _drs = drs.get(project, 'default')
    if not isinstance(input_file, six.string_types):
        if _drs in input_file:
            input_file = input_file[_drs]
        else:
            raise KeyError(
                'drs {} for {} project not specified for input_file '
                'in config-developer file'.format(_drs, project))
    filename = replace_tags(input_file, variable, i=j)
    return filename


def get_input_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
    all_files = []
    dirname_templates = get_input_dirname_template(variable, rootpath, drs)
    valid_dirs = []

    for dirname_template in dirname_templates:
        # Find latest version if required
        if '[latestversion]' not in dirname_template:
            valid_dirs.append(dirname_template)
        else:
            part1, part2 = dirname_template.split('[latestversion]')
            part2 = part2.lstrip(os.sep)
            if os.path.exists(part1):
                list_versions = os.listdir(part1)
                list_versions.sort(reverse=True)
                if 'latest' in list_versions:
                    list_versions.insert(
                        0, list_versions.pop(list_versions.index('latest')))
                for version in list_versions:
                    dirname = os.path.join(part1, version, part2)
                    if os.path.isdir(dirname):
                        valid_dirs.append(dirname)
                        break
            else:
                raise IOError('Path {} does not exist'.format(part1))

    # Set the filename glob
    filename_glob = _get_filename(variable, drs)

    for dir_name in valid_dirs:
        # Find files
        files = find_files(dir_name, filename_glob)

        # Select files within the required time interval
        files = select_files(files, variable['start_year'],
                             variable['end_year'])
        all_files.extend(files)

    return all_files


def get_input_fx_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
    dirname_templates = get_input_fx_dirname_template(variable, rootpath, drs)
    fx_files = {}

    for j, dirname_template in zip(
            range(len(dirname_templates)), dirname_templates):
        # Find latest version if required
        if '[latestversion]' in dirname_template:
            part1, part2 = dirname_template.split('[latestversion]')
            part2 = part2.lstrip(os.sep)
            # root part1 could not exist at all
            if not os.path.exists(part1):
                fx_files[variable['fx_files'][j]] = None
                return fx_files
            list_versions = os.listdir(part1)
            list_versions.sort(reverse=True)
            if 'latest' in list_versions:
                list_versions.insert(
                    0, list_versions.pop(list_versions.index('latest')))
            for version in list_versions:
                if version == 'latest':
                    dirname = os.path.join(part1, version, part2)
                    if os.path.isdir(dirname):
                        break
                else:
                    dirname = os.path.join(part1, version, part2)
                    if os.path.isdir(dirname):
                        break
        else:
            dirname = dirname_template

        # Set the filename glob
        filename_glob = _get_fx_filename(variable, drs, j)

        # Find files
        fx_file_list = find_files(dirname, filename_glob)
        if fx_file_list:
            # Grab the first file only; fx vars should have a single file
            fx_files[variable['fx_files'][j]] = fx_file_list[0]
        else:
            # No files
            fx_files[variable['fx_files'][j]] = None

    return fx_files


def get_output_file(variable, preproc_dir):
    """Return the full path to the output (preprocessed) file"""
    cfg = get_project_config(variable['project'])

    outfile = os.path.join(
        preproc_dir,
        '{diagnostic}_{preprocessor}_{short_name}'.format(**variable),
        replace_tags(cfg['output_file'], variable) + '.nc')

    return outfile


def get_statistic_output_file(variable, statistic, preproc_dir):
    """Get multi model statistic filename depending on settings"""
    values = dict(variable)
    values['stat'] = statistic.title()

    template = os.path.join(
        preproc_dir,
        '{diagnostic}_{preprocessor}_{short_name}',
        'MultiModel{stat}_{field}_{short_name}_{start_year}-{end_year}.nc',
    )

    outfile = template.format(**values)

    return outfile
