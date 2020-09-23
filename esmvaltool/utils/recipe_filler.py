import argparse
import os
from glob import glob

import yaml
from esmvalcore._config import read_config_developer_file
from esmvalcore.cmor.table import CMOR_TABLES

dataset_order = [
    'dataset', 'project', 'exp', 'mip', 'historical', 'ensemble', 'grid',
    'start_year', 'end_year'
]

# The base dictionairy (all wildcards):
base_dict = {
    #    '[institute]':     '*',
    #    '[dataset]':     '*',
    #    '[project]': '*',
    #    '[exp]':     '*',
    #    '[frequency]':    '*',
    #    '[ensemble]':     '*',
    #    '[mip]':     '*',
    #    '[modeling_realm]':     '*',
    #    '[short_name]': '*',
    #    '[grid]':     '*',
    #    '[latestversion]': 'latest',
    #    '[start_year]': '*',
    #    '[end_year]': '*',
    '{institute}': '*',
    '{dataset}': '*',
    '{project}': '*',
    '{exp}': '*',
    '{frequency}': '*',
    '{ensemble}': '*',
    '{mip}': '*',
    '{modeling_realm}': '*',
    '{short_name}': '*',
    '{grid}': '*',
    '{latestversion}': 'latest',
    '{start_year}': '*',
    '{end_year}': '*',
    '{activity}': '*',
    #
}


def rm_brackets(key):
    """
    Removes [ and  ] from a string
    """
    return key.replace('[', '').replace(']', '').replace('{',
                                                         '').replace('}', '')


def add_brackets(key):
    """
    Add { and } from to string
    """
    return ''.join(['{', key, '}'])


def get_site(CMIP='CMIP5'):
    """
    Get site (drs) from config-user.yml
    """
    config_yml = './config-user.yml'
    yamlconf = yaml.load(open(config_yml, 'r'))
    return yamlconf['drs'][CMIP]


def get_rootpath(CMIP='CMIP5'):
    """
    Get Basepath from config-user.yml
    """
    config_yml = './config-user.yml'
    yamlconf = yaml.load(open(config_yml, 'r'))
    return yamlconf['rootpath'][CMIP]


def get_input_dir(CMIP='CMIP5'):
    """
    Get input_dir from config-developer.yml
    """
    site = get_site(CMIP)
    yamlconf = read_config_developer_file()
    return yamlconf[CMIP]['input_dir'][site]


def get_input_file(CMIP='CMIP5'):
    """
    Get input_file from config-developer.yml
    """
    yamlconf = read_config_developer_file()
    return yamlconf[CMIP]['input_file']


def determine_basepath(CMIP='CMIP5'):
    """
    Determine a basepath
    """
    basepath = '/'.join([
        get_rootpath(CMIP=CMIP),
        get_input_dir(CMIP=CMIP),
        get_input_file(CMIP=CMIP)
    ])
    while basepath.find('//') > -1:
        basepath = basepath.replace('//', '/')
    return basepath


def list_all_files(file_dict, CMIP, debug=True):
    """
    List all files that match the dataset dictionairy.
    """
    basepath = determine_basepath(CMIP=CMIP)
    mip = file_dict['{mip}']
    short_name = file_dict['{short_name}']
    frequency = CMOR_TABLES['CMIP5'].get_variable(mip, short_name).frequency
    realm = CMOR_TABLES['CMIP5'].get_variable(mip, short_name).modeling_realm
    file_dict['{frequency}'] = frequency
    if len(realm) == 1:
        file_dict['{modeling_realm}'] = realm[0]

    # load all the files in the custom dict:
    new_path = basepath[:]
    for key, value in file_dict.items():
        new_path = new_path.replace(key, str(value))

    # Globs all the wildcards into a list of files.
    files = glob(new_path)
    return files


def get_time_range(fn_key):
    """
    Get a time range in years from the file name.
    """
    fn_key = fn_key.replace('.nc', '')
    start_time, end_time = fn_key.split('-')
    return start_time[:4], end_time[:4]


def file_to_recipe_dataset(fn, CMIP, file_dict):
    """
    Converts a filename to an recipe ready dataset.
    """
    # Add the obvious ones - ie the one you requested!
    output_dataset = {}
    output_dataset['project'] = CMIP
    for key, value in file_dict.items():
        if value == '*':
            continue
        key = rm_brackets(key)
        if key in dataset_order:
            output_dataset[key] = value

    # Split file name and base path into directory structure and filenames.
    basepath, basefile = os.path.split(determine_basepath(CMIP=CMIP))
    fnpath, fnfile = os.path.split(fn)

    # Assume basepath is separated by '/'
    basepath_split = basepath.split('/')
    fnpath_split = fnpath.split('/')

    # Check lengths of split lists.
    if len(basepath_split) != len(fnpath_split):
        assert "Paths do not match in length"

    # iterate through directory structure looking for useful bits.
    for base_key, fn_key in zip(basepath_split, fnpath_split):
        base_key = rm_brackets(base_key)
        if base_key not in dataset_order:
            continue
        output_dataset[base_key] = fn_key

    # Some of the key words include the splitting character '_' !
    basefile = basefile.replace('short_name', 'shortname')
    basefile = basefile.replace('start_year', 'startyear')
    basefile = basefile.replace('end_year', 'endyear')

    # Assume filename is separated by '_'
    basefile_split = basefile.split('_')
    fnfile_split = fnfile.split('_')

    # Check lengths of split lists.
    if len(basefile_split) != len(fnfile_split):
        assert "File names do not match in length!"

    # iterate through directory structure looking for useful bits.
    for base_key, fn_key in zip(basefile_split, fnfile_split):
        base_key = rm_brackets(base_key)
        # Check to make sure you're not adding a CMIP5 dataset into a CMIP6.
        if base_key == 'ensemble':
            if CMIP == 'CMIP5' and len(fn_key) != 6:
                return None
            if CMIP == 'CMIP6' and len(fn_key) != 8:
                return None

        # Special case for start end time

        if base_key == '*.nc':
            start_year, end_year = get_time_range(fn_key)
            output_dataset['start_year'] = start_year
            output_dataset['end_year'] = end_year
            assert 0

        output_dataset[base_key] = fn_key

    return output_dataset


def condense_times(add_datasets):
    """
    Condense the files to start_year and end_year
    Takes a list of fn:outs and joins together the ones with matching times
    """
    datasets = []
    seen = set()
    new_l = []
    for d in add_datasets:
        d.pop("shortname", None)
        d.pop("ensemble*.nc", None)
        t = tuple(d.items())
        if t not in seen:
            seen.add(t)
            datasets.append(d)

    return datasets


def get_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Tool to add additional_datasets to a recipe.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('recipe', help='Path or name of the yaml recipe file')
    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.path.dirname(__file__),
                                             'config-user.yml'),
                        help='Config file')

    parser.add_argument('-o', '--output', help='Output recipe')

    args = parser.parse_args()
    return args


def parse_recipe_to_dicts(recipe):
    """
    Parses a recipe's variables into a dictionary of dictionairies.
    """
    yamlrecipe = yaml.load(open(recipe, 'r'))

    output_dicts = {}

    # For each diagnostic in the recipe:
    for diag in yamlrecipe['diagnostics']:
        # For each variable in the diagnostic:
        for variable, var_dict in yamlrecipe['diagnostics'][diag][
                'variables'].items():
            new_dict = base_dict.copy()
            # for each key in the variable:
            for var_key, var_value in var_dict.items():
                var_key = add_brackets(var_key)
                if var_key in new_dict:
                    new_dict[var_key] = var_value
            output_dicts[(diag, variable)] = new_dict
    return output_dicts


def add_datasets_into_recipe(additional_datasets, recipe):
    """
    Add the datasets into a new recipe.

    """
    with open(recipe, 'r') as yamlfile:
        cur_yaml = yaml.safe_load(yamlfile)
        for diag_var, add_dat in additional_datasets.items():
            if 'additional_datasets' in cur_yaml['diagnostics']:
                cur_yaml['diagnostics'][
                    diag_var[0]]['additional_datasets'].extend(add_dat)
            else:
                cur_yaml['diagnostics'][
                    diag_var[0]]['additional_datasets'] = add_dat
    if cur_yaml:
        with open(recipe, 'w') as yamlfile:
            yaml.safe_dump(cur_yaml, yamlfile)


def run():
    """
    Run the `recipe_filler` program,.

    This tool takes a recipe without datasets provided, and adds an
    additional_datasets section for each diagnostics for each variable.

    It looks at whether any specific requirements are needed for the variable,
    and adds them.

    This won't work for derived objects.

    This won't check whether any datasets are provided already.
    """
    # Get arguments
    args = get_args()
    recipe = args.recipe

    # Convert the recipe into a yaml dict.
    recipe_dicts = parse_recipe_to_dicts(recipe)

    # Create a list of additional_datasets for each diagnostic/variable.
    additional_datasets = {}
    for (diag, variable), recipe_dict in recipe_dicts.items():
        files_recipe_dicts = {}
        recipe_dict['{short_name}'] = variable
        for CMIP in ['CMIP5']:
            files = list_all_files(recipe_dict, CMIP)
            add_datasets = []
            for fn in sorted(files):
                out = file_to_recipe_dataset(fn, CMIP, recipe_dict)
                if out is None:
                    continue
                add_datasets.append(out)
        additional_datasets[(diag, variable)] = condense_times(add_datasets)

    # add datasets to recipe as additional_datasets
    add_datasets_into_recipe(additional_datasets, recipe)


if __name__ == "__main__":
    run()
