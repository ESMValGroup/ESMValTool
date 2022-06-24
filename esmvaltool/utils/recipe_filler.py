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

- you can add as many variable parameters as are needed; if not added, the
  tool will use the "*" wildcard and find all available combinations;
- you can restrict the number of datasets to be looked for with the `dataset:`
  key for each variable, pass a list of datasets as value, e.g.
  `dataset: [MPI-ESM1-2-LR, MPI-ESM-LR]`;
- you can specify a pair of experiments eg `exp: [rcp26, rcp85]`
  for each variable; this will look for each available dataset per experiment
  and assemble an aggregated data stretch from each experiment; equivalent to
  esmvaltool's syntax of multiple experiments; this option needs an ensemble
  to be declared explicitly; it will return no entry if there are gaps in data
- `start_year` and `end_year` are mandatory and are used to filter out the
  datasets that don't have data in the interval; if you want all possible years
  hence no filtering on years just use "*" for start and end years;
- `config-user: rootpath: CMIPX` may be a list, rootpath lists are supported;

Caveats:

- the tool doesn't yet work for derived variables;
- operation restricted to CMIP data.

Have fun!
"""
import argparse
import datetime
import itertools
import logging
import logging.config
import os
import shutil
import time
import warnings
from glob import glob
from pathlib import Path

import yaml
from ruamel.yaml import YAML

import esmvalcore
from esmvalcore.cmor.table import CMOR_TABLES, read_cmor_tables

logger = logging.getLogger(__name__)

CFG = {}


# standard libs from esmvalcore ported here to avoid private func import
def load_config_developer(cfg_file=None):
    """Load the config developer file and initialize CMOR tables."""
    cfg_developer = read_config_developer_file(cfg_file)
    for key, value in cfg_developer.items():
        CFG[key] = value
    read_cmor_tables(CFG)


def _purge_file_handlers(cfg: dict) -> None:
    """Remove handlers with filename set.

    This is used to remove file handlers which require an output
    directory to be set.
    """
    cfg['handlers'] = {
        name: handler
        for name, handler in cfg['handlers'].items()
        if 'filename' not in handler
    }
    prev_root = cfg['root']['handlers']
    cfg['root']['handlers'] = [
        name for name in prev_root if name in cfg['handlers']
    ]


def _update_stream_level(cfg: dict, level=None):
    """Update the log level for the stream handlers."""
    handlers = cfg['handlers']

    for handler in handlers.values():
        if level is not None and 'stream' in handler:
            if handler['stream'] in ('ext://sys.stdout', 'ext://sys.stderr'):
                handler['level'] = level.upper()


def _get_log_files(cfg: dict, output_dir: str = None) -> list:
    """Initialize log files for the file handlers."""
    log_files = []

    handlers = cfg['handlers']

    for handler in handlers.values():
        filename = handler.get('filename', None)

        if filename:
            if not os.path.isabs(filename):
                handler['filename'] = os.path.join(output_dir, filename)
            log_files.append(handler['filename'])

    return log_files


def configure_logging(cfg_file: str = None,
                      output_dir: str = None,
                      console_log_level: str = None) -> list:
    """Configure logging.

    Parameters
    ----------
    cfg_file : str, optional
        Logging config file. If `None`, defaults to `configure-logging.yml`
    output_dir : str, optional
        Output directory for the log files. If `None`, log only to the console.
    console_log_level : str, optional
        If `None`, use the default (INFO).

    Returns
    -------
    log_files : list
        Filenames that will be logged to.
    """
    if cfg_file is None:
        cfg_loc = Path(esmvalcore.__file__ + "esmvalcore")
        # TODO change to new location of config module in 2.3.0
        cfg_file = cfg_loc.parents[0] / '_config' / 'config-logging.yml'

    cfg_file = Path(cfg_file).absolute()

    with open(cfg_file) as file_handler:
        cfg = yaml.safe_load(file_handler)

    if output_dir is None:
        _purge_file_handlers(cfg)

    log_files = _get_log_files(cfg, output_dir=output_dir)
    _update_stream_level(cfg, level=console_log_level)

    logging.config.dictConfig(cfg)
    logging.Formatter.converter = time.gmtime
    logging.captureWarnings(True)

    return log_files


def read_config_developer_file(cfg_file=None):
    """Read the developer's configuration file."""
    if cfg_file is None:
        cfg_loc = Path(esmvalcore.__file__ + "esmvalcore")
        cfg_file = cfg_loc.parents[0] / 'config-developer.yml'

    with open(cfg_file, 'r') as file:
        cfg = yaml.safe_load(file)

    return cfg


def _normalize_path(path):
    """Normalize paths.

    Expand ~ character and environment variables and convert path to absolute.

    Parameters
    ----------
    path: str
        Original path

    Returns
    -------
    str:
        Normalized path
    """
    if path is None:
        return None
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def read_config_user_file(config_file, folder_name, options=None):
    """Read config user file and store settings in a dictionary."""
    if not config_file:
        config_file = '~/.esmvaltool/config-user.yml'
    config_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(config_file)))
    # Read user config file
    if not os.path.exists(config_file):
        print(f"ERROR: Config file {config_file} does not exist")

    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    # DEPRECATED: remove in v2.4
    for setting in ('write_plots', 'write_netcdf'):
        if setting in cfg:
            msg = (
                f"Using '{setting}' in {config_file} is deprecated and will "
                "be removed in ESMValCore version 2.4. For diagnostics "
                "that support this setting, it should be set in the "
                "diagnostic script section of the recipe instead. "
                f"Remove the setting from {config_file} to get rid of this "
                "warning message.")
            print(f"Warning: {msg}")
            warnings.warn(DeprecationWarning(msg))

    if options is None:
        options = dict()
    for key, value in options.items():
        cfg[key] = value
        # DEPRECATED: remove in v2.4
        if key in ('write_plots', 'write_netcdf'):
            msg = (
                f"Setting '{key}' from the command line is deprecated and "
                "will be removed in ESMValCore version 2.4. For diagnostics "
                "that support this setting, it should be set in the "
                "diagnostic script section of the recipe instead.")
            print(f"Warning: {msg}")
            warnings.warn(DeprecationWarning(msg))

    # set defaults
    defaults = {
        'compress_netcdf': False,
        'exit_on_warning': False,
        'output_file_type': 'png',
        'output_dir': 'esmvaltool_output',
        'auxiliary_data_dir': 'auxiliary_data',
        'save_intermediary_cubes': False,
        'remove_preproc_dir': True,
        'max_parallel_tasks': None,
        'run_diagnostic': True,
        'profile_diagnostic': False,
        'config_developer_file': None,
        'drs': {},
        # DEPRECATED: remove default settings below in v2.4
        'write_plots': True,
        'write_netcdf': True,
    }

    for key in defaults:
        if key not in cfg:
            logger.info(
                "No %s specification in config file, "
                "defaulting to %s", key, defaults[key])
            cfg[key] = defaults[key]

    cfg['output_dir'] = _normalize_path(cfg['output_dir'])
    cfg['auxiliary_data_dir'] = _normalize_path(cfg['auxiliary_data_dir'])

    cfg['config_developer_file'] = _normalize_path(
        cfg['config_developer_file'])

    for key in cfg['rootpath']:
        root = cfg['rootpath'][key]
        if isinstance(root, str):
            cfg['rootpath'][key] = [_normalize_path(root)]
        else:
            cfg['rootpath'][key] = [_normalize_path(path) for path in root]

    # insert a directory date_time_recipe_usertag in the output paths
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join((folder_name, now))
    cfg['output_dir'] = os.path.join(cfg['output_dir'], new_subdir)

    # create subdirectories
    cfg['preproc_dir'] = os.path.join(cfg['output_dir'], 'preproc')
    cfg['work_dir'] = os.path.join(cfg['output_dir'], 'work')
    cfg['plot_dir'] = os.path.join(cfg['output_dir'], 'plots')
    cfg['run_dir'] = os.path.join(cfg['output_dir'], 'run')

    # Read developer configuration file
    load_config_developer(cfg['config_developer_file'])

    return cfg


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
    'dataset', 'project', 'exp', 'mip', 'ensemble', 'grid', 'start_year',
    'end_year'
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


def _get_site_rootpath(cmip_era):
    """Get site (drs) from config-user.yml."""
    config_yml = get_args().config_file
    with open(config_yml, 'r') as yamf:
        yamlconf = yaml.safe_load(yamf)
    drs = yamlconf['drs'][cmip_era]
    rootdir = yamlconf['rootpath'][cmip_era]
    logger.debug("%s root directory %s", cmip_era, rootdir)
    if drs == 'default' and 'default' in yamlconf['rootpath']:
        rootdir = yamlconf['rootpath']['default']
        logger.debug("Using drs default and "
                     "default: %s data directory", rootdir)

    return drs, rootdir


def _get_input_dir(cmip_era):
    """Get input_dir from config-developer.yml."""
    site = _get_site_rootpath(cmip_era)[0]
    yamlconf = read_config_developer_file()

    return yamlconf[cmip_era]['input_dir'][site]


def _get_input_file(cmip_era):
    """Get input_file from config-developer.yml."""
    yamlconf = read_config_developer_file()
    return yamlconf[cmip_era]['input_file']


def _determine_basepath(cmip_era):
    """Determine a basepath."""
    if isinstance(_get_site_rootpath(cmip_era)[1], list):
        rootpaths = _get_site_rootpath(cmip_era)[1]
    else:
        rootpaths = [_get_site_rootpath(cmip_era)[1]]
    basepaths = []
    for rootpath in rootpaths:
        if _get_input_dir(cmip_era) != os.path.sep:
            basepath = os.path.join(rootpath, _get_input_dir(cmip_era),
                                    _get_input_file(cmip_era))
        else:
            basepath = os.path.join(rootpath, _get_input_file(cmip_era))
        basepath = basepath.replace('//', '/')
        basepaths.append(basepath)
    logger.debug("We will look for files of patterns %s", basepaths)

    return basepaths


def _overlapping_datasets(files, all_years, start_year, end_year):
    """Process overlapping datasets and check for avail data in time range."""
    valid_files = []
    ay_sorted = sorted(all_years)
    if ay_sorted[0] <= start_year and ay_sorted[-1] >= end_year:
        yr_pairs = sorted(
            [all_years[i:i + 2] for i in range(0, len(all_years), 2)])
        yr_pairs = list(k for k, _ in itertools.groupby(yr_pairs))
        d_y = [
            yr_pairs[j][1] - yr_pairs[j + 1][0]
            for j in range(len(yr_pairs) - 1)
        ]
        gaps = [c for c in d_y if c < -1]
        if not gaps:
            valid_files = files
            logger.info("Contiguous data from multiple experiments.")
        else:
            logger.warning("Data from multiple exps has >1 year gaps! ")
            logger.debug("Start %s/end %s requested - "
                         "files covering %s found.",
                         start_year, end_year, yr_pairs)

    return valid_files


def filter_years(files, start_year, end_year, overlap=False):
    """
    Filter out files that are outside requested time range.

    Nifty function that takes a list of files and two years
    as arguments; it will build a series of filter dictionaries
    and check if data is available for the entire interval;
    it will return a single file per dataset, the first file
    in the list of files that cover the specified interval;
    optional argument `overlap` used if multiple experiments are
    used and overlap between datasets is present.

    Parameters
    ----------
    files: list
        A list of files that need filtering by requested time range.

    start_year: int
        Integer start year of requested range.

    end_year: int
        Integer end year of requested range.

    overlap: bool
        Flag if datasets overlap; defaults to False.

    Returns
    -------
    list
        List of files which have been identified as falling in
        the requested time range; if multiple files within time range
        per dataset, the first file will be returned.

    """
    valid_files = []
    available_years = {}

    if start_year == "*" and end_year == "*":
        return files

    if not files:
        return valid_files

    all_files_roots = [("").join(fil.split("_")[0:-1]) for fil in files]
    for fil in files:
        available_years[("").join(fil.split("_")[0:-1])] = []
    for fil in files:
        available_years[("").join(fil.split("_")[0:-1])].append(
            fil.split("_")[-1].strip(".nc").split("-"))

    all_years = []
    for root, yr_list in available_years.items():
        actual_years = []
        yr_list = list(itertools.chain.from_iterable(yr_list))
        for year in yr_list:
            if len(year) == 4:
                actual_years.append(int(year))
            else:
                actual_years.append(int(year[0:4]))
        actual_years = sorted(actual_years)
        all_years.extend(actual_years)
        if not overlap:
            actual_years = sorted(list(set(actual_years)))
            if actual_years[0] <= start_year and actual_years[-1] >= end_year:
                idx = all_files_roots.index(root)
                valid_files.append(files[idx])

    # multiple experiments to complete each other
    if overlap:
        valid_files = _overlapping_datasets(files, all_years, start_year,
                                            end_year)

    if not valid_files:
        logger.warning("No data found to fully cover start "
                       "%s / end %s as requested!", start_year, end_year)

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
    """
    List all files that match the dataset dictionary.

    Function that returnes all files that are determined by a
    file_dict dictionary; file_dict is keyed on usual parameters
    like `dataset`, `project`, `mip` etc; glob.glob is used
    to find files; speedup is achieved by replacing wildcards
    with values from CMOR tables.

    Parameters
    ----------
    file_dict: dict
        Dictionary to hold dataset specifications.

    cmip_era: str
        Either CMIP5 or CMIP6.

    Returns
    -------
    list:
        List of found files.

    """
    mip = file_dict['mip']
    short_name = file_dict['short_name']
    try:
        frequency = CMOR_TABLES[cmip_era].get_variable(mip,
                                                       short_name).frequency
        realms = CMOR_TABLES[cmip_era].get_variable(mip,
                                                    short_name).modeling_realm
    except AttributeError:
        logger.warning("Could not find %s CMOR table "
                       "for variable %s with mip %s",
                       cmip_era, short_name, mip)
        return []
    file_dict['frequency'] = frequency

    basepaths = _determine_basepath(cmip_era)
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
                    raise ValueError(
                        "Could not expand ~ to user home dir "
                        "please expand it in the config user file!")
                logger.info("Expanding path to %s", new_path)

            # Globs all the wildcards into a list of files.
            files = glob(new_path)
            all_files.extend(files)
    if not all_files:
        logger.warning("Could not find any file for data specifications.")

    return all_files


def _file_to_recipe_dataset(fn_path, cmip_era, file_dict):
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
    basefiles = _determine_basepath(cmip_era)
    _, fnfile = os.path.split(fn_path)

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
            elif base_key == "shortname":
                pass
            else:
                output_dataset[base_key] = fn_key
    if "exp" in file_dict:
        if isinstance(file_dict["exp"], list):
            output_dataset["exp"] = file_dict["exp"]

    return output_dataset


def _remove_duplicates(add_datasets):
    """
    Remove accidental duplicates.

    Close to 0% chances this will ever be used.
    May be used when there are actual duplicates in data
    storage, we've seen these before, but seldom.
    """
    datasets = []
    seen = set()

    for dataset in add_datasets:
        orig_exp = dataset["exp"]
        dataset["exp"] = str(dataset["exp"])
        tup_dat = tuple(dataset.items())
        if tup_dat not in seen:
            seen.add(tup_dat)
            dataset["exp"] = orig_exp
            datasets.append(dataset)

    return datasets


def _check_recipe(recipe_dict):
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
            if "exp" in var_pars:
                if isinstance(var_pars["exp"],
                              list) and "ensemble" not in var_pars:
                    logger.error("Asking for experiments list for ")
                    logger.error("variable %s - you need to ", var_name)
                    logger.error("define an ensemble for this case.")
                    do_exit = True
    if do_exit:
        raise ValueError("Please fix the issues in recipe and rerun")


def _check_config_file(user_config_file):
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
            logger.error("Config file missing rootpath for %s", proj)
            do_exit = True
        if proj not in user_config_file["drs"].keys():
            logger.error("Config file missing drs for %s", proj)
            do_exit = True
    if do_exit:
        raise ValueError("Please fix issues in config file and rerun")


def _parse_recipe_to_dicts(yamlrecipe):
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


def _add_datasets_into_recipe(additional_datasets, output_recipe):
    """Add the datasets into a new recipe."""
    yaml = YAML()
    yaml.default_flow_style = False
    with open(output_recipe, 'r') as yamlfile:
        cur_yaml = yaml.load(yamlfile)
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
            yaml.dump(cur_yaml, yamlfile)


def _find_all_datasets(recipe_dict, cmip_eras):
    """Find all datasets explicitly."""
    datasets = []
    for cmip_era in cmip_eras:
        if cmip_era == "CMIP6":
            activity = "CMIP"
        else:
            activity = ""
        drs, site_path = _get_site_rootpath(cmip_era)
        if drs in ["default", "SMHI"]:
            logger.info("DRS is %s; filter on dataset disabled.", drs)
            datasets = ["*"]
        else:
            if not isinstance(site_path, list):
                site_path = [site_path]
            for site_pth in site_path:
                if drs in ["BADC", "DKRZ", "CP4CDS"]:
                    institutes_path = os.path.join(site_pth, activity)
                elif drs in ["ETHZ", "RCAST"]:
                    exp = recipe_dict["exp"][0]
                    if exp == "*":
                        exp = "piControl"  # all institutes have piControl
                    mip = recipe_dict["mip"]
                    var = recipe_dict["short_name"]
                    institutes_path = os.path.join(site_pth, exp, mip, var)

                if not os.path.isdir(institutes_path):
                    logger.warning("Path to data %s "
                                   "does not exist; will look everywhere.",
                                   institutes_path)
                    datasets = ["*"]
                    return datasets

                institutes = os.listdir(institutes_path)
                if drs in ["BADC", "DKRZ", "CP4CDS"]:
                    for institute in institutes:
                        datasets.extend(
                            os.listdir(os.path.join(institutes_path,
                                                    institute)))
                else:
                    datasets.extend(institutes)

    return datasets


def _get_exp(recipe_dict):
    """Get the correct exp as list of single or multiple exps."""
    if isinstance(recipe_dict["exp"], list):
        exps_list = recipe_dict["exp"]
        logger.info("Multiple %s experiments requested", exps_list)
    else:
        exps_list = [recipe_dict["exp"]]
        logger.info("Single %s experiment requested", exps_list)

    return exps_list


def _get_datasets(recipe_dict, cmip_eras):
    """Get the correct datasets as list if needed."""
    if recipe_dict["dataset"] == "*":
        datasets = _find_all_datasets(recipe_dict, cmip_eras)
        return datasets
    if isinstance(recipe_dict['dataset'], list):
        datasets = recipe_dict['dataset']
        logger.info("Multiple %s datasets requested", datasets)
    else:
        datasets = [recipe_dict['dataset']]
        logger.info("Single %s dataset requested", datasets)

    return datasets


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


def _get_timefiltered_files(recipe_dict, exps_list, cmip_era):
    """Obtain all files that correspond to requested time range."""
    # multiple experiments allowed, complement data from each exp
    if len(exps_list) > 1:
        files = []
        for exp in exps_list:
            recipe_dict["exp"] = exp
            files.extend(list_all_files(recipe_dict, cmip_era))
        files = filter_years(files,
                             recipe_dict["start_year"],
                             recipe_dict["end_year"],
                             overlap=True)
        recipe_dict["exp"] = exps_list

    else:
        files = list_all_files(recipe_dict, cmip_era)
        files = filter_years(files, recipe_dict["start_year"],
                             recipe_dict["end_year"])

    return files


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
    logger.info("Using user configuration file: %s", args.config_file)
    logger.info("Using pilot recipe file: %s", input_recipe)
    logger.info("Writing filled out recipe to: %s", output_recipe)
    log_files = "\n".join(log_files)
    logger.info("Writing program log files to:\n%s", log_files)

    # check config user file
    _check_config_file(config_user)

    # parse recipe
    with open(input_recipe, 'r') as yamlfile:
        yamlrecipe = yaml.safe_load(yamlfile)
        _check_recipe(yamlrecipe)
        recipe_dicts = _parse_recipe_to_dicts(yamlrecipe)

    # Create a list of additional_datasets for each diagnostic/variable.
    additional_datasets = {}
    for (diag, variable), recipe_dict in recipe_dicts.items():
        logger.info("Looking for data for "
                    "variable %s in diagnostic %s", variable, diag)
        new_datasets = []
        if "short_name" not in recipe_dict:
            recipe_dict['short_name'] = variable
        elif recipe_dict['short_name'] == "*":
            recipe_dict['short_name'] = variable

        # adjust cmip era if needed
        if recipe_dict['project'] != "*":
            cmip_eras = [recipe_dict['project']]

        # get datasets depending on user request; always a list
        datasets = _get_datasets(recipe_dict, cmip_eras)

        # get experiments depending on user request; always a list
        exps_list = _get_exp(recipe_dict)

        # loop through datasets
        for dataset in datasets:
            recipe_dict['dataset'] = dataset
            logger.info("Seeking data for dataset: %s", dataset)
            for cmip_era in cmip_eras:
                files = _get_timefiltered_files(recipe_dict, exps_list,
                                                cmip_era)

                # assemble in new recipe
                add_datasets = []
                for fn in sorted(files):
                    fn_dir = os.path.dirname(fn)
                    logger.info("Data directory: %s", fn_dir)
                    out = _file_to_recipe_dataset(fn, cmip_era, recipe_dict)
                    logger.info("New recipe entry: %s", out)
                    if out is None:
                        continue
                    add_datasets.append(out)
                new_datasets.extend(add_datasets)
        additional_datasets[(diag, variable, cmip_era)] = \
            _remove_duplicates(new_datasets)

    # add datasets to recipe as additional_datasets
    shutil.copyfile(input_recipe, output_recipe, follow_symlinks=True)
    _add_datasets_into_recipe(additional_datasets, output_recipe)
    logger.info("Finished recipe filler. Go get some science done now!")


if __name__ == "__main__":
    run()
