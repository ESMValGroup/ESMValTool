"""Functions for downloading climate data files"""
import logging
import os
import subprocess

from .._data_finder import get_start_end_year, select_files

logger = logging.getLogger(__name__)


def synda_search(variable):
    """Search files using synda."""
    query = {
        'model': variable.get('model'),
        'project': variable.get('project'),
        'cmor_table': variable.get('mip'),
        'ensemble': variable.get('ensemble'),
        'experiment': variable.get('exp'),
        'variable': variable.get('short_name'),
    }

    query = {facet: value for facet, value in query.items() if value}

    query = ("{}='{}'".format(facet, value) for facet, value in query.items())

    cmd = ['synda', 'search', '--file']
    cmd.extend(query)
    cmd = ' '.join(cmd)
    logger.debug("Running: %s", cmd)
    result = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    logger.debug('Result:\n%s', result.strip())

    files = (l.split()[-1] for l in result.split('\n') if l.startswith('new'))
    files = select_files(files, variable['start_year'], variable['end_year'])

    # filter partially overlapping files
    intervals = {get_start_end_year(name): name for name in files}
    files = []
    for (start, end), filename in intervals.items():
        for _start, _end in intervals:
            if start == _start and end == _end:
                continue
            if start >= _start and end <= _end:
                break
        else:
            files.append(filename)

    logger.debug("Selected files:\n%s", '\n'.join(files))

    return files


def synda_download(synda_name, dest_folder):
    """Download file using synda."""
    filename = '.'.join(synda_name.split('.')[-2:])
    local_file = os.path.join(dest_folder, filename)

    if not os.path.exists(local_file):
        cmd = [
            'synda', 'get', '--dest_folder={}'.format(dest_folder),
            '--verify_checksum', synda_name
        ]
        cmd = ' '.join(cmd)
        logger.debug("Running: %s", cmd)
        subprocess.check_call(cmd, shell=True)

    return local_file


def download(files, dest_folder):
    """Download files that are not available locally"""
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)

    local_files = []
    for name in files:
        local_file = synda_download(synda_name=name, dest_folder=dest_folder)
        local_files.append(local_file)

    return local_files
