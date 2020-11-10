"""
Fill in a blank recipe with additional datasets.

Tool to obtain a set of additional datasets when given a blank recipe.
The blank recipe should contain, to the very least, a list of diagnostics
each with their variable(s). Example of minimum settings:

diagnostics:
  diagnostic:
    variables:
      ta:
        mip: Amon
        start_year: 1850
        end_year: 1900

Note that the tool will exit if any of these minimum settings are missing!

Key features:

- you can add as many variable parameters as are needed; if not added, the tool
  will use the "*" wildcard and find all available combinations;
- you can restrict the number of datasets to be looked for with the `dataset:`
  key for each variable, pass a list of datasets as value, e.g.
  `dataset: [MPI-ESM1-2-LR, MPI-ESM-LR]`;
- `start_year` and `end_year` are mandatory and are used to filter out the
  datasets that don't have data in the interval;
- `config-user: rootpath: CMIPX` may be a list, rootpath lists are supported;

Caveats:

- the tool doesn't yet work for derived variables;
- operation restricted to CMIP data.

Have fun!
"""
import argparse
import itertools
import logging
import os
import shutil
import sys
from glob import glob

import yaml
from esmvalcore._config import (configure_logging, read_config_developer_file,
                                read_config_user_file)
from esmvalcore.cmor.table import CMOR_TABLES

logger = logging.getLogger(__name__)

HEADER = r"""
______________________________________________________________________
          _____ ____  __  ____     __    _ _____           _
         | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
         |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
         | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
         |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
______________________________________________________________________

""" + __doc__

dataset_order = [
    'dataset', 'project', 'exp', 'mip', 'ensemble', 'grid',
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
    'start_year': '*',
    'end_year': '*',
    'activity': '*',
}


def get_site_rootpath(cmip_era):
    """Get site (drs) from config-user.yml."""
    config_yml = get_args().config_file
    with open(config_yml, 'r') as yamf:
        yamlconf = yaml.safe_load(yamf)
    drs = yamlconf['drs'][cmip_era]
    rootdir = yamlconf['rootpath'][cmip_era]
    logger.debug("%s root directory %s" % (cmip_era, rootdir))
    if drs == 'default' and 'default' in yamlconf['rootpath']:
        rootdir = yamlconf['rootpath']['default']
        logger.debug(f"Using drs default and "
                     f"default: %s data directory" % rootdir)

    return drs, rootdir


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
    if isinstance(get_site_rootpath(cmip_era)[1], list):
        rootpaths = get_site_rootpath(cmip_era)[1]
    else:
        rootpaths = [get_site_rootpath(cmip_era)[1]]
    basepaths = []
    for rootpath in rootpaths:
        if get_input_dir(cmip_era) != os.path.sep:
            basepath = os.path.join(rootpath, get_input_dir(cmip_era),
                                    get_input_file(cmip_era))
        else:
            basepath = os.path.join(rootpath,
                                    get_input_file(cmip_era))
        while basepath.find('//') > -1:
            basepath = basepath.replace('//', '/')
        basepaths.append(basepath)
    logger.debug(f"We will look for files of patterns %s" % basepaths)

    return basepaths


def filter_years(files, start_year, end_year):
    """
    Filter out files that are outside time range.

    Nifty function that takes a list of files and two years
    as arguments; it will build a series of filter dictionaries
    and check if data is available for the entire interval;
    it will return a single file per dataset, the first file
    in the list of files that cover the specified interval.
    """
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


def _resolve_latestversion(dirname_template):
    """Resolve the 'latestversion' tag."""
    if '{latestversion}' not in dirname_template:
        return dirname_template

    # Find latest version
    part1, part2 = dirname_template.split('{latestversion}')
    part2 = part2.lstrip(os.sep)
    part1_contents = glob(part1)
    if part1_contents:
        versions = os.listdir(part1_contents[0])
        versions.sort(reverse=True)
        for version in ['latest'] + versions:
            dirname = os.path.join(part1, version, part2)
            if glob(dirname):
                return dirname

    return dirname_template


def list_all_files(file_dict, cmip_era):
    """List all files that match the dataset dictionary."""
    mip = file_dict['mip']
    short_name = file_dict['short_name']
    try:
        frequency = CMOR_TABLES[cmip_era].get_variable(mip,
                                                       short_name).frequency
        realms = CMOR_TABLES[cmip_era].get_variable(mip,
                                                    short_name).modeling_realm
    except AttributeError:
        logger.warning(f"Could not find %s CMOR table "
                       f"for variable %s with mip %s" %
                       (cmip_era, short_name, mip))
        return []
    file_dict['frequency'] = frequency

    basepaths = determine_basepath(cmip_era)
    all_files = []

    for basepath in basepaths:
        new_path = basepath[:]

        # could have multiple realms
        for realm in realms:
            file_dict['modeling_realm'] = realm

            # load all the files in the custom dict
            for key, value in file_dict.items():
                new_path = new_path.replace('{' + key + '}', str(value))
            new_path = _resolve_latestversion(new_path)
            if new_path.startswith("~"):
                new_path = os.path.expanduser(new_path)
                if not new_path.startswith(os.sep):
                    logger.error(f"Could not expand ~ to user home dir "
                                 f"please expand it in the config user file!")
                    sys.exit(1)
                logger.warning("Expanding path to %s" % new_path)

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
    basefiles = determine_basepath(cmip_era)
    _, fnfile = os.path.split(fn)

    for basefile in basefiles:
        _, basefile = os.path.split(basefile)
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
    """
    Remove accidental duplicates.

    Close to 0% chances this will ever be used.
    May be used when there are actual duplicates in data
    storage, we've seen these before, but seldom.
    """
    datasets = []
    seen = set()

    for dataset in add_datasets:
        tup_dat = tuple(dataset.items())
        if tup_dat not in seen:
            seen.add(tup_dat)
            datasets.append(dataset)

    return datasets


def check_recipe(recipe_dict):
    """Perform a quick recipe check for mandatory fields."""
    do_exit = False
    if "diagnostics" not in recipe_dict:
        logger.error("Recipe missing diagnostics section.")
        do_exit = True
    for diag_name, diag in recipe_dict["diagnostics"].items():
        if "variables" not in diag:
            logger.error("Diagnostic %s missing variables.", diag_name)
            do_exit = True
        for var_name, var_pars in diag["variables"].items():
            if "mip" not in var_pars:
                logger.error("Variable %s missing mip.", var_name)
                do_exit = True
            if "start_year" not in var_pars:
                logger.error("Variable %s missing start_year.", var_name)
                do_exit = True
            if "end_year" not in var_pars:
                logger.error("Variable %s missing end_year.", var_name)
                do_exit = True
    if do_exit:
        logger.error("Please fix the issues in recipe and rerun. Exiting.")
        sys.exit(1)


def check_config_file(user_config_file):
    """Perform a quick recipe check for mandatory fields."""
    do_exit = False
    if "rootpath" not in user_config_file:
        logger.error("Config file missing rootpath section.")
        do_exit = True
    if "drs" not in user_config_file:
        logger.error("Config file missing drs section.")
        do_exit = True
    for proj in cmip_eras:
        if proj not in user_config_file["rootpath"].keys():
            logger.error("Config file missing rootpath for %s" % proj)
            do_exit = True
        if proj not in user_config_file["drs"].keys():
            logger.error("Config file missing drs for %s" % proj)
            do_exit = True
    if do_exit:
        logger.error("Please fix issues in config file and rerun. Exiting.")
        sys.exit(1)


def parse_recipe_to_dicts(yamlrecipe):
    """Parse a recipe's variables into a dictionary of dictionairies."""
    output_dicts = {}
    for diag in yamlrecipe['diagnostics']:
        for variable, var_dict in yamlrecipe['diagnostics'][diag][
                'variables'].items():
            new_dict = base_dict.copy()
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


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('recipe', help='Path/name of yaml pilot recipe file')
    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.environ["HOME"], '.esmvaltool',
                                             'config-user.yml'),
                        help='User configuration file')

    parser.add_argument('-o',
                        '--output',
                        default=os.path.join(os.getcwd(),
                                             'recipe_autofilled.yml'),
                        help='Output recipe, default recipe_autofilled.yml')

    args = parser.parse_args()
    return args


def run():
    """Run the `recipe_filler` tool. Help in __doc__ and via --help."""
    # Get arguments
    args = get_args()
    input_recipe = args.recipe
    output_recipe = args.output
    cmip_eras = ["CMIP5", "CMIP6"]

    # read the config file
    config_user = read_config_user_file(args.config_file,
                                        'recipe_filler',
                                        options={})

    # configure logger
    run_dir = os.path.join(config_user['output_dir'], 'recipe_filler')
    if not os.path.isdir(run_dir):
        os.makedirs(run_dir)
    log_files = configure_logging(output_dir=run_dir,
                                  console_log_level=config_user['log_level'])
    logger.info(HEADER)
    logger.info("Using user configuration file: %s" % args.config_file)
    logger.info("Using pilot recipe file: %s" % input_recipe)
    logger.info("Writing filled out recipe to: %s" % output_recipe)
    logger.info("Writing program log files to:\n%s", "\n".join(log_files))

    # check config user file
    check_config_file(config_user)

    # parse recipe
    with open(input_recipe, 'r') as yamlfile:
        yamlrecipe = yaml.safe_load(yamlfile)
        check_recipe(yamlrecipe)
        recipe_dicts = parse_recipe_to_dicts(yamlrecipe)

    # Create a list of additional_datasets for each diagnostic/variable.
    additional_datasets = {}
    for (diag, variable), recipe_dict in recipe_dicts.items():
        logger.info("Looking for data for variable %s in diagnostic %s",
                    variable, diag)
        new_datasets = []
        if "short_name" not in recipe_dict:
            recipe_dict['short_name'] = variable

        # filter on user request
        if isinstance(recipe_dict['dataset'], list):
            datasets = recipe_dict['dataset']
        else:
            datasets = [recipe_dict['dataset']]
        if recipe_dict['project'] != "*":
            cmip_eras = [recipe_dict['project']]

        logger.info("Seeking data for datasets: %s", str(datasets))
        for dataset in datasets:
            recipe_dict['dataset'] = dataset
            logger.info("Analyzing dataset: %s", dataset)
            for cmip_era in cmip_eras:
                files = list_all_files(recipe_dict, cmip_era)
                add_datasets = []
                for fn in sorted(files):
                    logger.info("Data directory: %s", os.path.dirname(fn))
                    out = file_to_recipe_dataset(fn, cmip_era, recipe_dict)
                    logger.info("New recipe entry: %s", out)
                    if out is None:
                        continue
                    add_datasets.append(out)
                new_datasets.extend(add_datasets)
        additional_datasets[(diag, variable, cmip_era)] = \
            remove_duplicates(new_datasets)

    # add datasets to recipe as additional_datasets
    shutil.copyfile(input_recipe, output_recipe, follow_symlinks=True)
    add_datasets_into_recipe(additional_datasets, output_recipe)
    logger.info("Finished recipe filler. Go get some science done now!")


if __name__ == "__main__":
    run()
