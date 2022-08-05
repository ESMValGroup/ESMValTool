import os
import argparse
import itertools
import logging
import pyesgf.search as pys
import re
import yaml

from esmvalcore.esgf._search import get_esgf_facets
from esmvalcore.esgf.facets import DATASET_MAP, FACETS


# Intial testing done with an input type like:
# Input "recipe" file
# datasets:
#  - {dataset: BCC-ESM1, project: CMIP6, exp: historical, ensemble: r1i1p1f1, mip: Amon, grid: gn, short_name: tas}
#  - {dataset: CanESM2, project: CMIP5, exp: historical, ensemble: r1i1p1, mip: Amon, short_name: tas}

# set up logging
logger = logging.getLogger(__name__)


def load_esgf_pyclient_config():
    """Load a basic esgf-pyclient configuration."""
    cfg = {
        'logon': {
            'interactive': False,
            'bootstrap': True,
        },
        'search_connection': {
            'url': 'http://esgf-node.llnl.gov/esg-search',
            'distrib': True,
            'timeout': 120,
            'cache': '~/.esmvaltool/cache/pyesgf-search-results',
            'expire_after': 86400,
        },
    }

    return cfg


def get_connection():
    """Connect to ESGF."""
    cfg = load_esgf_pyclient_config()
    connection = pys.SearchConnection(**cfg["search_connection"])
    return connection


def list_cmip6_grids():
    connection = get_connection()
    context = connection.new_context(project='CMIP6', latest=True)
    return list(context.facet_counts['grid_label'])


def get_facet_counts(facets):
    connection = get_connection()
    esgf_facets = get_esgf_facets(facets)
    context = connection.new_context(**esgf_facets, latest=True)
    return context.facet_counts


def get_facet_options(variable, facet):
    project = variable['project']
    connection = get_connection()
    esgf_facets = get_esgf_facets(variable)
    if isinstance(facet, list) and len(facet) == 1:
        facet = facet[0]
    esgf_facet = FACETS[project][facet]
    context = connection.new_context(**esgf_facets,
                                     facets=[esgf_facet],
                                     latest=True)
    result = list(context.facet_counts[esgf_facet])
    if facet == 'dataset':
        reverse_dataset_map = {
            p: {v: k
                for k, v in DATASET_MAP[p].items()}
            for p in DATASET_MAP
        }
        for i, dataset in enumerate(result):
            result[i] = reverse_dataset_map[project].get(dataset, dataset)

    result.sort()

    return result


def format_dataset(project, dataset, ensemble, dataset_dict, grid=None):
    ds = dataset_dict["dataset"]
    exp = dataset_dict["exp"]
    mip = dataset_dict["mip"]
    short_name = dataset_dict["short_name"]
    if project == 'CMIP5':
        if ds == dataset:
            logger.info(
                " - {dataset: %s, project: CMIP5, exp: %s, ensemble: %s, mip: %s, short_name: %s}",
                ds, exp, ensemble, mip, short_name)
    elif project == 'CMIP6':
        grid = dataset_dict["grid"]
        if ds == dataset:
            logger.info(
                " - {dataset: %s, project: CMIP6, exp: %s, ensemble: %s, mip: %s, short_name: %s, grid: %s}",
                ds, exp, ensemble, mip, short_name, grid)


def search_variables(variables, facet):
    """Search datasets that have all variables."""
    result = None
    for variable in variables:
        options = get_facet_options(variable, facet)
        if result is None:
            result = set(options)
        else:
            result = result.intersection(options)

    return sorted(result)


def group_ensembles(ensembles):
    ensembles = [
        tuple(int(i) for i in re.findall(r'\d+', ens)) for ens in ensembles
    ]

    for i in range(len(ensembles[0])):
        ensembles = _group_ensembles(ensembles, i)
    ensembles.sort()

    groups = []
    for ens in ensembles:
        txt = ''
        for name, value in zip('ripf', ens):
            txt += name
            if value[0] == value[1]:
                txt += f"{value[0]}"
            else:
                txt += f"({value[0]}:{value[1]})"
        groups.append(txt)

    return groups


def _group_ensembles(ensembles, i=0):

    def order(ens):
        prefix, suffix = ens[:i], ens[i + 1:]
        return (prefix, suffix, ens[i])

    def grouper(ens):
        prefix, suffix = ens[:i], ens[i + 1:]
        return (prefix, suffix)

    grouped_ensembles = []
    ensembles = sorted(ensembles, key=order)
    for (prefix, suffix), bunch in itertools.groupby(ensembles, key=grouper):
        bunch = list(bunch)
        prev = bunch[0][i]
        groups = [[prev]]
        for ensemble in bunch[1:]:
            if ensemble[i] == prev + 1:
                prev += 1
            else:
                groups[-1].append(prev)
                prev = ensemble[i]
                groups.append([prev])
        groups[-1].append(prev)
        result = []
        for group in groups:
            item = prefix + (group, ) + suffix
            result.append(item)
        grouped_ensembles.extend(result)

    return grouped_ensembles


def _expand_dataset(dataset):
    """Load a recipe or a dict containing datasets info."""
    grids = []
    # expand from dataset dict
    project = dataset["project"]
    exps = dataset["exp"]
    if not isinstance(exps, list):
        exps = [exps]
    mip = dataset["mip"]
    if not isinstance(mip, list):
        mips = [mip]
    if "grid" in dataset:
        grids = dataset["grid"]
        if not isinstance(grids, list):
            grids = [grids]
    short_names = dataset["short_name"]
    if not isinstance(short_names, list):
        short_names = [short_names]

    return project, exps, mips, grids, short_names


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--recipe', required=False,
                        help='Path/name of yaml pilot recipe file')
    parser.add_argument('-o',
                        '--output',
                        required=False,
                        default=os.path.join(os.getcwd(),
                                             'recipe_autofilled.yml'),
                        help='Output recipe, default recipe_autofilled.yml')

    args = parser.parse_args()

    return args


def _parse_recipe_to_dicts(yamlrecipe):
    """Parse a recipe's variables into a dictionary of dictionairies."""
    datasets = yamlrecipe["datasets"]

    return datasets


def _assemble_variable(project, exps, mips, grids, short_names):
    """Get a dict for variable."""
    variable = {
        'project': project,
        'exp': exps,
        'mip': mips,
        'short_name': short_names,
        'grid': grids,
    }

    facets = []
    to_remove = []

    # remove all values of ["*"] or None
    # and add to facets if ["*"]
    for key, val in variable.items():
        if not val:
            to_remove.append(key)
        elif isinstance(val, list) and val[0] == "*":
            to_remove.append(key)
            facets.append(key)
    if to_remove:
        for key in to_remove:
            del variable[key]

    return variable, facets


def _get_results(project, exps,
                 mips, grids,
                 short_names):
    """Get ESGF search results."""
    variables = []
    all_facets = ["dataset"]
    results = []

    for exp in exps:
        for mip in mips:
            for short_name in short_names:
                variable, facets = _assemble_variable(project, exps,
                                                      mips, grids,
                                                      short_names)
                variables.append(variable)
                all_facets.extend(facets)

    logger.info("We will use the following facets %s for search", str(facets))
    datasets = sorted(search_variables(variables, ['dataset']),
                                       key=str.lower)
    for dataset in datasets:
        if dataset == 'fio-esm':
            dataset = 'FIO-ESM'
        for variable in variables:
            variable['dataset'] = dataset
        if not all_facets:
            results.append([dataset])
        else:
            for facet in all_facets:
                result = search_variables(variables, facet)
                results.append(result)

    return results, all_facets


def search_datasets_esgf(dataset_list=None):
    """Search the ESGF database for requested datasets."""
    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    logging.getLogger().setLevel('info'.upper())

    # load inputs
    args = get_args()
    input_recipe = args.recipe
    output_recipe = args.output
    if input_recipe is not None:
        logger.info("Running from the command line.")
        with open(input_recipe, 'r') as yamlfile:
            yamlrecipe = yaml.safe_load(yamlfile)
            dataset_list = _parse_recipe_to_dicts(yamlrecipe)
    else:
        logger.info("Running via module import.")

    logger.info("Found the following datasets on ESGF from your recipe:")

    if dataset_list:
        for dataset_dict in dataset_list:
            results = []

            # expand from dataset dict
            project, exps, mips, grids, short_names = \
                _expand_dataset(dataset_dict)

            # run the thing
            if project == 'CMIP5':
                results, facets = _get_results(project, exps,
                                               mips, grids,
                                               short_names)
            elif project == 'CMIP6':
                for grid in grids:
                    results, facets = _get_results(project, exps,
                                                   mips, grids,
                                                   short_names)

            print(facets)
            if results:
                for result in results:
                    dataset, ensemble, grid = result
                    format_dataset(project, dataset, ensemble,
                                   dataset_dict, grid)


if __name__ == "__main__":
    search_datasets_esgf()
