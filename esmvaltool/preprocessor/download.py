import subprocess
import logging
import os

from ..interface_scripts.data_finder import get_input_filename

logger = logging.getLogger(__name__)


def synda_search(model, variable):
    """Search files using synda."""
    mip = model['mip'] if 'mip' in model else variable.get('mip')
    query = {
        'model': model.get('name'),
        'project': model.get('project'),
        'cmor_table': mip,
        'ensemble': model.get('ensemble'),
        'experiment': model.get('exp'),
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
            filename = line.split()[-1]
            name = os.path.splitext(filename)[0]
            start, end = name.split('_')[-1].split('-')
            start_year, end_year = int(start[:4]), int(end[:4])
            if (start_year <= model['end_year']
                    and end_year >= model['start_year']):
                files.append(filename)

    return files


def synda_download(filename, dest_folder):
    """Download file using synda."""
    cmd = [
        'synda', 'get', '--dest_folder={}'.format(dest_folder),
        '--verify_checksum', filename
    ]
    cmd = ' '.join(cmd)
    logger.debug("Running: %s", cmd)
    subprocess.check_call(cmd, shell=True)


def download(files, model, variable, rootpath, drs):
    """Download files that are not available locally"""
    local_dir = os.path.dirname(
        get_input_filename(
            model=model, variable=variable, rootpath=rootpath, drs=drs))
    os.makedirs(local_dir, exist_ok=True)

    local_files = []
    for name in files:
        filename = name[name.index(variable['short_name'] + '_'):]
        local_file = os.path.join(local_dir, filename)
        local_files.append(local_file)
        if not os.path.exists(local_file):
            synda_download(filename=name, dest_folder=local_dir)

    return local_files
