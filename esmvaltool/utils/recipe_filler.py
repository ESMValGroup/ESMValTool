"""Fill in recipe."""
import argparse
import itertools
import os
import shutil
from glob import glob

import yaml
from esmvalcore._config import read_config_developer_file
from esmvalcore.cmor.table import CMOR_TABLES

dataset_order = [
    'dataset', 'project', 'exp', 'mip', 'historical', 'ensemble', 'grid',
    'start_year', 'end_year'
]

# cmip eras
cmip_eras = ["CMIP5", "CMIP6"]

# The base dictionairy (all wildcards):
base_dict = {
    'institute': '*',
    'dataset': '*',
    'project': '*',
    'exp': '*',
    'frequency': '*',
    'ensemble': '*',
    'mip': '*',
    'modeling_realm': '*',
    'short_name': '*',
    'grid': '*',
    'latestversion': 'latest',
    'start_year': '*',
    'end_year': '*',
    'activity': '*',
}


def get_site_rootpath(cmip_era):
    """Get site (drs) from config-user.yml."""
    config_yml = get_args().config_file
    yamlconf = yaml.load(open(config_yml, 'r'))
    return yamlconf['drs'][cmip_era], yamlconf['rootpath'][cmip_era]


def get_input_dir(cmip_era):
    """Get input_dir from config-developer.yml."""
    site = get_site_rootpath(cmip_era)[0]
    yamlconf = read_config_developer_file()
    return yamlconf[cmip_era]['input_dir'][site]


def get_input_file(cmip_era):
    """Get input_file from config-developer.yml."""
    yamlconf = read_config_developer_file()
    return yamlconf[cmip_era]['input_file']


def determine_basepath(cmip_era):
    """Determine a basepath."""
    basepath = os.path.join(
        get_site_rootpath(cmip_era)[1], get_input_dir(cmip_era),
        get_input_file(cmip_era))
    while basepath.find('//') > -1:
        basepath = basepath.replace('//', '/')

    return basepath


def filter_years(files, start_year, end_year):
    """Filter out files that are outside time range."""
    valid_files = []
    available_years = {}
    all_files_roots = [("").join(fil.split("_")[0:-1]) for fil in files]
    for fil in files:
        available_years[("").join(fil.split("_")[0:-1])] = []
    for fil in files:
        available_years[("").join(fil.split("_")[0:-1])].append(
            fil.split("_")[-1].strip(".nc").split("-"))

    for root, yr_list in available_years.items():
        yr_list = list(itertools.chain.from_iterable(yr_list))
        actual_years = []
        for year in yr_list:
            if len(year) == 4:
                actual_years.append(int(year))
            else:
                actual_years.append(int(year[0:4]))
        actual_years = sorted(list(set(actual_years)))
        if actual_years[0] <= start_year and actual_years[-1] >= end_year:
            idx = all_files_roots.index(root)
            valid_files.append(files[idx])

    return valid_files


def list_all_files(file_dict, cmip_era):
    """List all files that match the dataset dictionary."""
    basepath = determine_basepath(cmip_era)
    mip = file_dict['mip']
    short_name = file_dict['short_name']
    frequency = CMOR_TABLES[cmip_era].get_variable(mip, short_name).frequency
    realms = CMOR_TABLES[cmip_era].get_variable(mip, short_name).modeling_realm
    file_dict['frequency'] = frequency
    new_path = basepath[:]

    # could have multiple realms
    all_files = []
    for realm in realms:
        file_dict['modeling_realm'] = realm

        # load all the files in the custom dict
        for key, value in file_dict.items():
            new_path = new_path.replace('{' + key + '}', str(value))

        # Globs all the wildcards into a list of files.
        files = glob(new_path)
        all_files.extend(files)

    # filter time
    if all_files:
        all_files = filter_years(all_files, file_dict["start_year"],
                                 file_dict["end_year"])

    return all_files


def file_to_recipe_dataset(fn, cmip_era, file_dict):
    """Convert a filename to an recipe ready dataset."""
    # Add the obvious ones - ie the one you requested!
    output_dataset = {}
    output_dataset['project'] = cmip_era
    for key, value in file_dict.items():
        if value == '*':
            continue
        if key in dataset_order:
            output_dataset[key] = value

    # Split file name and base path into directory structure and filenames.
    _, basefile = os.path.split(determine_basepath(cmip_era))
    _, fnfile = os.path.split(fn)

    # Some of the key words include the splitting character '_' !
    basefile = basefile.replace('short_name', 'shortname')
    basefile = basefile.replace('start_year', 'startyear')
    basefile = basefile.replace('end_year', 'endyear')

    # Assume filename is separated by '_'
    basefile_split = [key.replace("{", "") for key in basefile.split('_')]
    basefile_split = [key.replace("}", "") for key in basefile_split]
    fnfile_split = fnfile.split('_')

    # iterate through directory structure looking for useful bits.
    for base_key, fn_key in zip(basefile_split, fnfile_split):
        if base_key == '*.nc':
            fn_key = fn_key.replace('.nc', '')
            start_year, end_year = fn_key.split('-')
            output_dataset['start_year'] = start_year
            output_dataset['end_year'] = end_year
        elif base_key == "ensemble*.nc":
            output_dataset['ensemble'] = fn_key
        elif base_key == "grid*.nc":
            output_dataset['grid'] = fn_key
        elif base_key not in ["shortname", "ensemble*.nc", "*.nc"]:
            output_dataset[base_key] = fn_key

    return output_dataset


def remove_duplicates(add_datasets):
    """Remove accidental duplicates."""
    datasets = []
    seen = set()
    for dataset in add_datasets:
        tup_dat = tuple(dataset.items())
        if tup_dat not in seen:
            seen.add(tup_dat)
            datasets.append(dataset)

    return datasets


def get_args():
    """Parse command line arguments."""
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
    """Parse a recipe's variables into a dictionary of dictionairies."""
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
                if var_key in new_dict:
                    new_dict[var_key] = var_value
            output_dicts[(diag, variable)] = new_dict
    return output_dicts


def add_datasets_into_recipe(additional_datasets, output_recipe):
    """Add the datasets into a new recipe."""
    with open(output_recipe, 'r') as yamlfile:
        cur_yaml = yaml.safe_load(yamlfile)
        for diag_var, add_dat in additional_datasets.items():
            if add_dat:
                if 'additional_datasets' in cur_yaml['diagnostics']:
                    cur_yaml['diagnostics'][diag_var[0]]['variables'][
                        diag_var[1]]['additional_datasets'].extend(add_dat)
                else:
                    cur_yaml['diagnostics'][diag_var[0]]['variables'][
                        diag_var[1]]['additional_datasets'] = add_dat
    if cur_yaml:
        with open(output_recipe, 'w') as yamlfile:
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
    input_recipe = args.recipe
    output_recipe = args.output

    # Convert the recipe into a yaml dict.
    recipe_dicts = parse_recipe_to_dicts(input_recipe)

    # Create a list of additional_datasets for each diagnostic/variable.
    additional_datasets = {}
    for (diag, variable), recipe_dict in recipe_dicts.items():
        recipe_dict['short_name'] = variable
        for cmip_era in cmip_eras:
            files = list_all_files(recipe_dict, cmip_era)
            add_datasets = []
            for fn in sorted(files):
                out = file_to_recipe_dataset(fn, cmip_era, recipe_dict)
                if out is None:
                    continue
                add_datasets.append(out)
            additional_datasets[(diag, variable, cmip_era)] = \
                remove_duplicates(add_datasets)

    # add datasets to recipe as additional_datasets
    shutil.copyfile(input_recipe, output_recipe, follow_symlinks=True)
    add_datasets_into_recipe(additional_datasets, output_recipe)


if __name__ == "__main__":
    run()
