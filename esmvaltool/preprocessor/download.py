import subprocess
import logging
import os

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
    logger.debug('Result: %s', result)

    files = []
    for line in result.split('\n'):
        if line.startswith('new'):
            synda_name = line.split()[-1]
            name = os.path.splitext(synda_name)[0]
            start, end = name.split('_')[-1].split('-')
            start_year, end_year = int(start[:4]), int(end[:4])
            if (start_year <= variable['end_year']
                    and end_year >= variable['start_year']):
                files.append(synda_name)

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
    os.makedirs(dest_folder, exist_ok=True)

    local_files = []
    for name in files:
        local_file = synda_download(synda_name=name, dest_folder=dest_folder)
        local_files.append(local_file)

    return local_files
