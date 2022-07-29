import itertools
import logging
import pyesgf.search as pys
import re

from esmvalcore.esgf._search import get_esgf_facets
from esmvalcore.esgf.facets import DATASET_MAP, FACETS


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


def format_dataset(project, dataset, ensemble, grid=None):
    if project == 'CMIP5':
        print(
            "  - {mip: Amon, project: CMIP5, exp: [historical, rcp85], dataset: "
            + dataset + ", ensemble: " + ensemble + "}")
    elif project == 'CMIP6':
        print(
            "  - {mip: Amon, project: CMIP6, exp: [historical, ssp585], dataset: "
            + dataset + ", ensemble: " + ensemble + ", grid: " + grid + "}")


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


def test_group_ensembles_cmip5():

    ensembles = [
        "r1i1p1",
        "r2i1p1",
        "r3i1p1",
        "r4i1p1",
        "r1i2p1",
    ]
    groups = group_ensembles(ensembles)
    expected = ['r1i2p1', 'r(1:4)i1p1']
    print(groups)
    print(expected)
    assert groups == expected


def test_group_ensembles_cmip6():

    ensembles = [
        "r1i1p1f1",
        "r4i1p1f1",
        "r3i1p2f1",
        "r4i1p2f1",
        "r3i1p1f1",
    ]
    groups = group_ensembles(ensembles)
    expected = ['r1i1p1f1', 'r(3:4)i1p(1:2)f1']
    print(groups)
    print(expected)
    assert groups == expected


if __name__ == '__main__':
    # test_group_ensembles_cmip5()
    # test_group_ensembles_cmip6()

    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    logging.getLogger().setLevel('info'.upper())

    project = 'CMIP5'
    variables = []
    for exp in 'historical', 'rcp85':
        for short_name in 'tas', 'pr':
            variable = {
                'project': project,
                'exp': exp,
                'mip': 'Amon',
                'short_name': short_name,
            }
            variables.append(variable)

    datasets = sorted(search_variables(variables, 'dataset'), key=str.lower)
    for dataset in datasets:
        for variable in variables:
            variable['dataset'] = dataset
        ensembles = search_variables(variables, 'ensemble')
        for ensemble in group_ensembles(ensembles):
            if dataset == 'fio-esm':
                dataset = 'FIO-ESM'
            format_dataset(project, dataset, ensemble)

    project = 'CMIP6'
    results = []
    for grid in ['gn', 'gr', 'gr1']:
        variables = []
        for exp in 'historical', 'ssp585':
            for short_name in 'tas', 'pr':
                variable = {
                    'project': project,
                    'exp': exp,
                    'mip': 'Amon',
                    'short_name': short_name,
                    'grid': grid,
                }
                variables.append(variable)

        datasets = sorted(search_variables(variables, 'dataset'),
                          key=str.lower)
        for dataset in datasets:
            for variable in variables:
                variable['dataset'] = dataset
            ensembles = search_variables(variables, 'ensemble')
            for ensemble in group_ensembles(ensembles):
                results.append([dataset, ensemble, grid])

    results.sort(key=lambda i: (i[0].lower(
    ), tuple(int(i) for i in re.findall(r'\d+', i[1]))))
    for dataset, ensemble, grid in results:
        format_dataset(project, dataset, ensemble, grid)
