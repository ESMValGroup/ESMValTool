"""
Call like:
python util/nml-utils/generateNML/generateNML.py --namelist <NAMELIST>

"""
import argparse
import os

import xmltodict
import sys

import logging

from collections import OrderedDict

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

NON_MODEL_LIST = [
    "OBS", "obs4mips", "OBS_gridfile", "reanalysis", "observation"
]

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
CONSOLE_HANDLER = logging.StreamHandler()
CONSOLE_HANDLER.setLevel(logging.DEBUG)
CONSOLE_HANDLER.setFormatter(
    logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(CONSOLE_HANDLER)


def get_modellines(infiles):
    modellines = list()
    for path_to_infile in infiles:
        if not path_to_infile:
            continue
        modellines.append(get_modelline(path_to_infile))
    return list(set(modellines))


def _get_times(filename):
    return tuple([
        item[0:4] for item in filename.split(".")[0].split('_')[-1].split('-')
    ])


def get_modelline(path_to_infile):
    d = dict()
    d['Amon'] = {'time_freq': 'mon', 'realm': 'atmos'}
    d['Omon'] = {'time_freq': 'mon', 'realm': 'ocean'}
    d['Lmon'] = {'time_freq': 'mon', 'realm': 'land'}
    d['AERmon'] = {'time_freq': 'mon', 'realm': 'aerosol'}

    parts = path_to_infile.split('/')
    start, end = _get_times(parts[-1])
    indices = [10, 8, 9, 10, 11]
    out = "ESGF_CMIP6 " + " ".join([
        parts[i] for i in indices
    ]) + " " + d[parts[13]]['time_freq'] + " " + d[parts[13]][
        'realm'] + " " + parts[13] + " " + parts[12] + " latest " + parts[
            15] + " " + start + " " + end + " CMIP6_template"
    logger.debug("Got this: %s", path_to_infile)
    logger.debug("Give this: %s", out)
    return out


def get_info_from_freva(**kwargs):
    facets = []
    freva_cmd = "freva"
    for key, value in kwargs.items():
        facets.append("{0}={1}".format(key, value))
    module = "module load cmip6-dicad/1.0"
    cmd = [
        "{0} &> /dev/null; {1} ".format(module, freva_cmd), "--databrowser",
        "project=cmip6"
    ] + facets
    cmd = " ".join(cmd)
    freva_out = subprocess.check_output(cmd, shell=True).decode().split('\n')
    logger.debug("Output of command '%s': '%s'", cmd, freva_out)
    return get_modellines(freva_out)


def _to_strings(content):
    if isinstance(content, list):
        out = list()
        for item in content:
            if isinstance(item, str) or isinstance(item, unicode):
                out.append(item)
            elif isinstance(item, OrderedDict) and '#text' in item.keys():
                out.append(item['#text'])
            else:
                logger.error("Wrong type %s", type(item))
                raise Exception
        return out
    else:
        logger.error("Wrong type %s", type(content))
        raise Exception


def get_available_dataset_info(requirements):
    out = list()
    for requirement in requirements:
        query = dict()
        for key, value in requirement.iteritems():
            query[key] = "'/{0}/'".format("|".join(_to_strings(value)))
        available_datasets = get_info_from_freva(**query)
        logger.debug("Available Datasets type '%s'", type(available_datasets))
        logger.debug("Available Datasets '%s'", available_datasets)
        out.append(available_datasets)
    return out


def _extract_non_models(liste):
    def _is_non_model(d):
        return any([item in NON_MODEL_LIST for item in d.split()])

    return [item for item in liste if _is_non_model(item)]


def get_namelist(namelist):

    _check_namelist(namelist)

    requirements_per_diagblock = get_namelist_diag_requirements(namelist)
    logger.debug("Requirements per diagblock '%s'", requirements_per_diagblock)
    available_datasets_per_diagblock = get_available_dataset_info(
        requirements_per_diagblock)
    if not available_datasets_per_diagblock:
        logger.debug("No data available")

    with open(namelist, 'r') as f:
        j = xmltodict.parse(f.read())

    if j['namelist']['MODELS'] is not None:
        j['namelist']['MODELS'] = ["{{ global_modelline }}"]

    logger.debug("DIAGNOSTICS %s", j['namelist']['DIAGNOSTICS']['diag'])
    logger.debug("Available Datasets %s", available_datasets_per_diagblock)

    if isinstance(j['namelist']['DIAGNOSTICS']['diag'], OrderedDict):
        non_models = _extract_non_models(
            j['namelist']['DIAGNOSTICS']['diag']['model'])
        logger.debug("Non models %s", non_models)
        j['namelist']['DIAGNOSTICS']['diag'][
            'model'] = available_datasets_per_diagblock[0] + non_models
    elif isinstance(j['namelist']['DIAGNOSTICS']['diag'], list):
        for i in range(len(j['namelist']['DIAGNOSTICS']['diag'])):
            non_models = _extract_non_models(
                j['namelist']['DIAGNOSTICS']['diag'][i]['model'])
            logger.debug("Non models %s", non_models)
            j['namelist']['DIAGNOSTICS']['diag'][i][
                'model'] = available_datasets_per_diagblock[i] + non_models
    else:
        raise Exception

    return xmltodict.unparse(j, pretty=True)


def _check_namelist(namelist):
    """Return True id namelist exist, else throw exception."""
    if not os.path.isfile(namelist):
        raise Exception
    return True


#def _get_variable_str(variable):
#    if isinstance(variable, str):
#        out = variable
#    elif isinstance(variable, list):
#        if variable:
#            try:
#                out = ";".join(variable)
#            except:
#                out = ";".join([item.__repr__() for item in variable])
#        else:
#            out = None
#    else:
#        try:
#            out = variable['#text']
#        except:
#            out = None
#    return out


def _get_experiments(modellines):
    """Get Experiments used in modellines."""
    return _get_unique_parts(modellines, flag=True)


def _get_cmor_table(modellines):
    """Get cmor_table used in modellines."""
    return _get_unique_parts(modellines, flag=False)


def _get_unique_parts(modellines, flag=True):
    valid_experiments = set([
        'historical', 'piControl', "1pctCO2", "esmFixClim1", "esmHistorical",
        "amip"
    ])
    valid_cmor_table = set(
        ['Amon', 'Omon', 'Lmon', 'aero', 'OImon', 'day', '3hr', 'AERmon'])

    match = list()
    out = list()

    if not isinstance(modellines, list):
        modellines = [modellines]
    for modelline in modellines:
        if not isinstance(modelline, str) and not isinstance(
                modelline, unicode):
            try:
                model_line_parts = modelline['#text'].split()
            except:
                logger.debug("Modelline is of type: '%s'", type(modelline))
                model_line_parts = []
        else:
            model_line_parts = modelline.split()

        if not model_line_parts:
            msg = "Empty modelline"
            logger.debug(msg)

        if any([item in NON_MODEL_LIST for item in model_line_parts]):
            continue

        if flag is True:
            match += list(
                set.intersection(set(model_line_parts), valid_experiments))
            target = "experiment"
        else:
            match += list(
                set.intersection(set(model_line_parts), valid_cmor_table))
            target = "cmor_table"

        if not match:
            msg = "Unknown {0} for modelline : {1}".format(target, modelline)
            logger.debug(msg)
        out.append(match)
    return list(set(match))


def get_namelist_diag_requirements(namelist):
    """Return requirements of each diagblock."""

    _check_namelist(namelist)

    out = list()

    with open(namelist, 'r') as f:
        j = xmltodict.parse(f.read())

    diagblocks = j['namelist']['DIAGNOSTICS']['diag']
    if not isinstance(diagblocks, list):
        diagblocks = [diagblocks]
    for diagblock in diagblocks:
        #variable = _get_variable_str(diagblock['variable'])
        if 'model' in diagblock.keys():
            experiment = _get_experiments(diagblock['model'])
            cmor_table = _get_cmor_table(diagblock['model'])
        else:
            experiment = []
            cmor_table = []
        out.append({
            'variable':
            diagblock['variable'] if isinstance(diagblock['variable'], list)
            else [diagblock['variable']],
            'experiment':
            experiment,
            'cmor_table':
            cmor_table
        })
    return out


def write_xml(filename, content, same_folder=False):
    head, tail = os.path.split(filename)
    if same_folder:
        newfilename = os.path.join(head, "new_" + tail)
    else:
        newfilename = "new_" + tail
    if os.path.exists(newfilename):
        logger.error("File %s already exists", newfilename)
        raise Exception
    with open(newfilename, 'w') as f:
        f.write(content)
    logger.debug("Created new file %s", newfilename)


def main():
    parser = argparse.ArgumentParser(
        description='Generate routine evaluation namelist.')
    parser.add_argument('--namelist', dest='namelist')

    args = parser.parse_args()

    write_xml(args.namelist, get_namelist(args.namelist))


if __name__ == "__main__":
    main()
