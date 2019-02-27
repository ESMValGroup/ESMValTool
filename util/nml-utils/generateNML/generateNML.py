"""
Call like:
python util/nml-utils/generateNML/generateNML.py --namelist <NAMELIST>

"""
import argparse
import os

from jinja2 import Template

import xmltodict
import sys
import yaml
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

T_MLINE = " ".join([
    '{{ project }}', '{{ name }}', '{{ product }}', '{{ institute }}',
    '{{ model }}', '{{ experiment }}', '{{ time_freq }}', '{{ realm }}',
    '{{ mip }}', '{{ ensemble }}', '{{ version }}', '{{ grid }}',
    '{{ start_year }}', '{{ end_year }}', '{{ ptid }}'
])


def get_modelline(**kwargs):
    d = dict()
    d['Amon'] = {'time_freq': 'mon', 'realm': 'atmos'}
    d['Omon'] = {'time_freq': 'mon', 'realm': 'ocean'}
    d['Lmon'] = {'time_freq': 'mon', 'realm': 'land'}

    kwargs.update({'version': 'latest', 'ptid': 'CMIP6_template'})
    if 'mip' in kwargs.keys():
        if kwargs['mip'] in d.keys():
            kwargs.update(d[kwargs['mip']])
    tt_mline = Template(T_MLINE)
    return tt_mline.render(**kwargs)

def get_info_from_freva(**kwargs):
    facets = []
    freva_cmd = "freva"
    #freva_cmd = "./mock_freva.sh"
    for key, value in kwargs.items():
        facets.append("{0}={1}".format(key,value))
    module = "module load cmip6-dicad/1.0"
    #cmd = "freva --databrowser project=cmip6 {0} --all-facets".format(" ".join(facets))
    cmd = ["{0} &> /dev/null; {1} ".format(module, freva_cmd), "--databrowser", "project=cmip6"] + facets + ["--all-facets"]
    cmd = " ".join(cmd)
    freva_out = subprocess.check_output(cmd, shell=True).decode()
    print(freva_out)
    print((freva_out.replace(":", ": [").replace("\n", "]\n")))
    return yaml.load((freva_out.replace(":", ": [").replace("\n", "]\n")))

def get_available_dataset_info(requirements):
    out = list()

    for diagnostic in requirements:
        for experiment in diagnostic['experiment']:
            available_datasets = get_info_from_freva(experiment=experiment)
            if len(available_datasets['model']) > 0:
                for model in available_datasets['model']:
                    for cmor_table in available_datasets['cmor_table']:
                        out.append({'model': model, 'experiment': experiment, 'cmor_table': cmor_table })
    return out

def get_namelist(namelist):


    _check_namelist(namelist)

    requirements_per_diagblock = get_namelist_diag_requirements(namelist)
    print(requirements_per_diagblock)
    available_datasets_per_diagblock = get_available_dataset_info(requirements_per_diagblock)
    if available_datasets_per_diagblock:
        print(type(available_datasets_per_diagblock[0]))
        print(available_datasets_per_diagblock[0])
    else:
        print("No data available")

    with open(namelist, 'r') as f:
        j = xmltodict.parse(f.read())

    if j['namelist']['MODELS'] is not None:
        j['namelist']['MODELS'] = ["{{ global_modelline }}"]

    for i in range(len(available_datasets_per_diagblock)):
        l = list()
        for item in available_datasets_per_diagblock[i]:
            if not isinstance(item, dict):
                print("Warning 1")
                continue
            l.append(get_modelline(**item))
        print(j['namelist']['DIAGNOSTICS']['diag'][0])
        j['namelist']['DIAGNOSTICS']['diag'][i]['model'] = l


    return xmltodict.unparse(j, pretty=True)


def _check_namelist(namelist):
    """Return True id namelist exist, else throw exception."""
    if not os.path.isfile(namelist):
        raise Exception
    return True


def _get_variable_str(variable):
    if isinstance(variable, str):
        out = variable
    elif isinstance(variable, list):
        if len(variable) != 0:
            try:
                out = ";".join(variable)
            except:
                out = ";".join([item.__repr__() for item in variable])
        else:
            out = None
    else:
        try:
            out = variable['#text']
        except:
            out = None
    return out


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
    valid_cmor_table = set([item.lower() for item in [
        'Amon','Omon','Lmon', 'aero', 'OImon', 'day', '3hr'
    ]])
    black_list = [
        "OBS", "obs4mips", "OBS_gridfile", "reanalysis", "observation"
    ]

    match = list()
    out = list()

    if not isinstance(modellines, list):
        modellines = [modellines]
    for modelline in modellines:
        if not isinstance(modelline, str):
            try:
                model_line_parts = modelline['#text'].split()
            except:
                print("Modelline is of type: {0}".format(type(modelline)))
                model_line_parts = []
        else:
            model_line_parts = modelline.split()

        if len(model_line_parts) == 0:
            msg = "Empty modelline"
            print(msg)

        if any([item in black_list for item in model_line_parts]):
            continue

        if flag is True:
            match += list(
                set.intersection(set(model_line_parts), valid_experiments))
            target = "experiment"
        else:
            match += list(
                set.intersection(set(model_line_parts), valid_cmor_table))
            target = "cmor_table"

        if len(match) == 0:
            msg = "Unknown {0} for modelline : {1}".format(target, modelline)
            print(msg)
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
    cnt = 0
    for diagblock in diagblocks:
        variable = _get_variable_str(diagblock['variable'])
        if 'model' in diagblock.keys():
            experiment = _get_experiments(diagblock['model'])
            cmor_table = _get_cmor_table(diagblock['model'])
        else:
            experiment = []
            cmor_table = []
        time_span = None
        out.append({'variables': diagblock['variable'] if isinstance(diagblock['variable'], list) else [diagblock['variable']], 'experiment': experiment, 'cmor_table': cmor_table})
        #print(
        #    "DiagBlock {0} of namelist {1} needs variable: {2}, experiment {3}, time_span {4}"
        #    .format(cnt, namelist, variable, experiment, time_span))
        cnt += 1
    return out


def main():
    parser = argparse.ArgumentParser(
        description='Generate routine evaluation namelist.')
    parser.add_argument('--namelist', dest='namelist')

    args = parser.parse_args()

    kwa = dict(args._get_kwargs())

    namelist = kwa['namelist']
    print(get_namelist(namelist))
    #print(get_template_string(namelist))
    ##requirements = get_namelist_diag_requirements(namelist)
    ##import json
    ##print(json.dumps(requirements, indent=4, separators=(',', ': ')))

if __name__ == "__main__":
    main()
