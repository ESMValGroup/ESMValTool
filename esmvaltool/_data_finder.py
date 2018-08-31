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


def replace_tags(path, variable, fx_var=None):
    """Replace tags in the config-developer's file with actual values."""
    path = path.strip('/')

    tlist = re.findall(r'\[([^]]*)\]', path)

    paths = [path]
    for tag in tlist:
        original_tag = tag
        tag, _, _ = _get_caps_options(tag)

        if tag == 'var':
            replacewith = variable['short_name']
        elif tag == 'fx_var':
            replacewith = fx_var
        elif tag == 'latestversion':  # handled separately later
            continue
        elif tag == 'tier':
            replacewith = ''.join(('Tier', str(variable['tier'])))
        else:  # all other cases use the corresponding dataset dictionary key
            if tag in variable:
                replacewith = variable[tag]
            else:
                raise KeyError(
                    "Dataset key {} must be specified for project {}, check "
                    "your recipe entry".format(tag, variable['project']))

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


def _get_dirnames(input_type, variable, rootpath, drs):
    """Return a the full paths to input directories of input_type."""
    project = variable['project']

    root = get_rootpath(rootpath, project)
    input_dir = _select_drs(input_type, drs, project)

    dirnames = []
    for dirname_template in replace_tags(input_dir, variable):
        dirname = _resolve_latestversion(dirname_template)
        dirname = os.path.join(root, dirname)
        if os.path.exists(dirname):
            logger.debug("Found %s", dirname)
            dirnames.append(dirname)
        else:
            logger.debug("Skipping non-existent %s", dirname)

    return dirnames


def get_input_dirnames(variable, rootpath, drs):
    """Return a the full paths to input directories."""
    return _get_dirnames('input_dir', variable, rootpath, drs)


def get_input_fx_dirnames(variable, rootpath, drs):
    """Return a the full paths to fx file directories."""
    dirnames = []
    for fx_file in variable['fx_files']:
        # Need to reassign the mip so we can find sftlf/of
        # make a copy of variable -> new_variable for this
        new_variable = dict(variable)
        new_variable['mip'] = replace_mip_fx(fx_file)
        dirnames.extend(_get_dirnames('fx_dir', new_variable, rootpath, drs))

    return dirnames


def _get_filename_glob(variable, drs):
    input_file = _select_drs('input_file', drs, variable['project'])
    filename = replace_tags(input_file, variable)[0]
    return filename


def _get_fx_filename_glob(variable, drs, fx_var):
    input_file = _select_drs('fx_file', drs, variable['project'])
    filename = replace_tags(input_file, variable, fx_var=fx_var)[0]
    return filename


def get_input_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
    logger.debug("Looking for input files for variable %s for dataset %s",
                 variable['short_name'], variable['dataset'])
    dirnames = get_input_dirnames(variable, rootpath, drs)
    filename_glob = _get_filename_glob(variable, drs)

    all_files = []
    for dir_name in dirnames:
        files = find_files(dir_name, filename_glob)
        # Select files within the required time interval
        files = select_files(files, variable['start_year'],
                             variable['end_year'])
        all_files.extend(files)

    return all_files


def get_input_fx_filelist(variable, rootpath, drs):
    """Return the full path to input files."""
    fx_files = {}
    for fx_var in variable['fx_files']:
        logger.debug("Looking for fx_files files of type %s for dataset %s",
                     fx_var, variable['dataset'])
        dirnames = get_input_fx_dirnames(variable, rootpath, drs)
        if not dirnames:
            # No files
            fx_files[fx_var] = None
        else:
            # Set the filename glob
            filename_glob = _get_fx_filename_glob(variable, drs, fx_var)

            # Find all possible files
            all_files = [
                find_files(dir_name, filename_glob) for dir_name in dirnames
            ]
            # filter out empty entries
            all_files = [l for l in all_files if l]
            if not all_files:
                fx_files[fx_var] = None
            else:
                # Keep only the first entry
                fx_files[fx_var] = [fx_ls[0] for fx_ls in all_files][0]

    return fx_files


def get_output_file(variable, preproc_dir):
    """Return the full path to the output (preprocessed) file"""
    cfg = get_project_config(variable['project'])

    outfile = os.path.join(
        preproc_dir,
        '{diagnostic}_{preprocessor}_{short_name}'.format(**variable),
        replace_tags(cfg['output_file'], variable)[0] + '.nc')

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
