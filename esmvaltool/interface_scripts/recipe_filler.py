#!/usr/bin/ipython
from glob import glob

import argparse
import os
import sys
import yaml
import numpy as np

import itertools

from  esmvalcore._config import read_config_developer_file
#from  esmvalcore._config import read_config_user_file

dataset_order = ['dataset', 'project', 'exp', 'mip', 'historical', 'ensemble', 'grid', 'start_year', 'end_year']

# The base dictionairy (all wildcards):
base_dict  = {
    '[institute]':     '*',
    '[dataset]':     '*',
    '[project]': '*',
    '[exp]':     '*',
    '[frequency]':    '*',
    '[ensemble]':     '*',
    '[mip]':     '*',
    '[modeling_realm]':     '*',
    '[short_name]': '*',
    '[grid]':     '*',
    '[latestversion]': '*',
    '[start_year]': '*',
    '[end_year]': '*'
    }

class hashabledict(dict):
    """
    A hashable dict that can be used as a index of another dict.
    """
    def __hash__(self):
        return hash(tuple(sorted(self.items())))

def rm_brackets(key):
    """
    Removes [ and  ] from a string
    """
    return key.replace('[', '').replace(']','')


def add_brackets(key):
    """
    Add [ and  ] from to string
    """
    return ''.join(['[', key, ']'])


def get_site(CMIP = 'CMIP5'):
    """
    Get site (drs) from config-user.yml
    """
    config_yml = './config-user.yml'
    yamlconf = yaml.load(open(config_yml, 'r'))
    return yamlconf['drs'][CMIP]


def get_rootpath(CMIP = 'CMIP5'):
    """
    Get Basepath from config-user.yml
    """
    config_yml = './config-user.yml'
    yamlconf = yaml.load(open(config_yml, 'r'))
    return yamlconf['rootpath'][CMIP]


def get_input_dir(CMIP = 'CMIP5'):
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
    basepath = '/'.join([get_rootpath(CMIP = CMIP),
                         get_input_dir(CMIP = CMIP),
                         get_input_file(CMIP = CMIP)])
    while basepath.find('//') > -1:
        basepath = basepath.replace('//', '/')
    return basepath


def list_all_files(file_dict, CMIP,debug = False):
    """
    List all files that match the dataset dictionairy.
    """
    basepath = determine_basepath(CMIP = CMIP)
    if debug:
        print(CMIP, "base path: ", basepath)

    # load all the files in the custom dict:
    new_path =  basepath[:]
    for key, value in file_dict.items():
        new_path = new_path.replace(key, value)


    # Globs all the wildcards into a list of files.
    if debug:
        print(CMIP, "New path: ", new_path)
    files = glob(new_path)
    return files


def convert_dict_to_recipe(file_dict, leadingspaces=0):
    """"
    Converts a dict into an additional_datasets line.
    """
    order = ['dataset', 'project', 'exp', 'mip', 'historical', 'ensemble',
             'grid', 'start_year', 'end_year']
    out = ''.join([' ' for itn in range(leadingspaces)])
    out = ''.join([out, '- {'])
    for key in order:
        if key in file_dict:
            out = ''.join([out, key, ': ', file_dict[key], ', '])
    out = ''.join([out, '}'])
    #print ("convert_dict_to_recipe:", out)
    return out


def get_time_range(fn_key):
    """
    Get a time range in years from the file name.
    """
    fn_key = fn_key.replace('.nc', '')
    start_time, end_time = fn_key.split('-')
    return start_time[:4], end_time[:4]


def file_to_recipe_dataset(fn, CMIP, file_dict, debug = False):
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
    basepath, basefile = os.path.split(determine_basepath(CMIP = CMIP))
    fnpath, fnfile = os.path.split(fn)

    # ###############
    # Look at directory structure
    if debug:
        print("basepath", basepath)
        print("fnpath", fnpath)

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

    # ###############
    # Look at file naming structure
    if debug:
        print("basefile", basefile)
        print("fnfile", fnfile)

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
                if debug:
                    print('Not a CMIP5 dataset:', fn, (fn_key))
                return None
            if CMIP == 'CMIP6' and len(fn_key) != 8:
                if debug:
                    print('Not a CMIP6 dataset:', fn, (fn_key))
                return None

        # Special case for start end time
        if base_key == '*.nc':
            start_year, end_year = get_time_range(fn_key)
            output_dataset['start_year'] = start_year
            output_dataset['end_year'] = end_year

        output_dataset[base_key] = fn_key

    return output_dataset


def condense_times(files_recipe_dicts):
    """
    Condense the files to start_year and end_year
    Takes a list of fn:outs and joins together the ones with matching times
    """
    match_by = dataset_order[:-2]

    groups = {}
    for fn, output_dataset in files_recipe_dicts.items():
        # tmp dict is a hashabledict, so it can be used as a index of
        # another dict.
        tmp_dict = hashabledict(output_dataset.copy())

        # Pop the start and end year
        start_year = tmp_dict.pop('start_year')
        end_year = tmp_dict.pop('end_year')
        nc = tmp_dict.pop('*.nc')

        # Group the dictionairies by everything except start and end year.
        if tmp_dict in groups.keys():
            groups[tmp_dict].append(int(start_year))
            groups[tmp_dict].append(int(end_year))
        else:
            groups[tmp_dict] = [int(start_year), int(end_year)]

    # Add the start and end year back in
    condense_times_dict = {}
    for output_dataset, years in groups.items():
        output_dataset['start_year'] = str(np.min(years))
        output_dataset['end_year'] = str(np.max(years))
        recipe_string = convert_dict_to_recipe(output_dataset)
        condense_times_dict[recipe_string] = output_dataset
    return condense_times_dict


def get_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Tool to add additional_datasets to a recipe.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('recipe', help='Path or name of the yaml recipe file')
    parser.add_argument(
        '-c',
        '--config-file',
        default=os.path.join(os.path.dirname(__file__), 'config-user.yml'),
        help='Config file')
    
    parser.add_argument(
        '-o',
        '--output',
        help='Output recipe')

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
        for variable, var_dict in yamlrecipe['diagnostics'][diag]['variables'].items():
            new_dict = base_dict.copy()
            new_dict['[short_name]'] = variable

            # for each key in the variable:
            for var_key, var_value in var_dict.items():
                var_key = add_brackets(var_key)
                if var_key in new_dict:
                    new_dict[var_key] = var_value
            output_dicts[(diag, variable)] = new_dict
    return output_dicts


def calc_leading_spaces(line):
    """
    Calculates the number of leading spaces in the line.
    """
    if line.find('\t')>-1:
        assert("This tool doesn't work with tabs! please use spaces.")
    return sum(1 for _ in itertools.takewhile(lambda c: c == ' ', line))

#
# def find_diag_in_recipe(recipe, diag, debug=True):
#     """
#     Find the line and line number of the diagnostic pair.
#
#     recipe is the loaded list of lines.
#     """
#     diag_indent = 0
#     diag_line = 0
#     for line_number, line in enumerate(recipe):
#         if debug: print(line_number, '\t', diag, variable)
#         if diag in line:
#             diag_indent = line.find(diag)
#             diag_line = line_number
#             if debug: print("Found diag:", line_number, diag, line)
#             continue
#         if not diag_indent:
#             # haven't found the diagnostic yet.
#             continue
#

def add_datasets_into_recipe(additional_datasets, debug=False):
    """
    Add the datasets into a new recipe.

    additional_datasets: dict
        This dict is
                key: (diagnostic name, variable name/short_name)
                value: condensed_times_dict
                condensed_times_dict is a dict:
                    key: string to add into recipe
                    value: dictionary of things to put into the recipe.

    This won't work if the variable isn't named as header in the recipe.
    ie:

      mfo:
        exp: historical

    should work, but

      mfo_diagnostic:
        short_name: mfo
        exp: historical

    probably won't work.

    short_name
    """
    recipe_fn = get_args().recipe
    recipe = open(recipe_fn, 'r')
    recipe_text = recipe.readlines()
    print (recipe_text)
    for (diag, variable), condensed_times_dict in additional_datasets.items():
        diag_indent = 0
        tmp_recipe_text = recipe_text[:]
        for line_number, line in enumerate(tmp_recipe_text):
            if debug: print(line_number, '\t', diag, variable)
            if diag in line:
                diag_indent = line.find(diag)
                if debug: print("Found diag:", line_number, diag, line)
                continue
            if not diag_indent:
                # haven't found the diagnostic yet.
                continue

            leading_spaces = calc_leading_spaces(line)
            if debug: print("Found leading_spaces:", line_number, line, leading_spaces)

            if leading_spaces > diag_indent:
                if debug:
                    print("Intended:",  line_number, leading_spaces, ">",
                          diag_indent)

                var_indent = 0
                if variable in line:
                    if debug:
                        print("Found variable:",  line_number, variable, line)

                    var_indent = line.find(variable)
                    white_space = ''.join([' ' for itr in range(var_indent)])
                    add_datasets_txt = white_space
                    add_datasets_txt += '  additional_datasets:\n'

                    for recipe_strings in sorted(condensed_times_dict):
                        add_datasets_txt += white_space + '    ' + recipe_strings + '\n'
                    #print('Adding:', add_datasets_txt, )
                    recipe_text.insert(line_number+1, add_datasets_txt)

                # Removing datsets specifiers which appear twice:
                if var_indent < leading_spaces:
                    #Iteratre though dataset fiels to check if they exist.
                    for dataset_field in dataset_order:
                        start_char = leading_spaces
                        end_char = leading_spaces+len(dataset_field)
                        line_dataset_field = line[start_char:end_char]
                        if line_dataset_field == dataset_field:
                            # Add 1 as we added an entry above.
                            if line == recipe_text[line_number+1]:
                                del recipe_text[line_number+1]

            if leading_spaces <= diag_indent:
                break

    try:
        output_file = get_args().output
    except: pass

    if output_file is None:
        output_file = recipe_fn.replace('.yml', '_filled.yml')


    print("Saving new recipe to:", output_file)
    f = open(output_file, "w")
    recipe_text = "".join(recipe_text)
    f.write(recipe_text)
    f.close()


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
        for CMIP in ['CMIP5', 'CMIP6']:
            files = list_all_files(recipe_dict, CMIP)
            for fn in sorted(files):
                out = file_to_recipe_dataset(fn, CMIP, recipe_dict)
                if out is None: continue
                files_recipe_dicts[fn] = out
        condensed_times_dict = condense_times(files_recipe_dicts)
        additional_datasets[(diag, variable)] = condensed_times_dict

        #print(diag, variable, convert_dict_to_recipe(condensed_times_dict))

    # We've now got all the files that fit the requirements of the recipe.
    # Sending it back into a new recipe.
    add_datasets_into_recipe(additional_datasets)


if __name__ == "__main__":
    run()
