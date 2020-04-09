"""Obtain dataset attributes from the dataset or its filename

Where to get the information for the matching from
Valid values: 'attributes', 'filename'.

Use 'attributes' to match the datasets by information from
attributes, specifically 'realization', 'physics_version',
'initaliazation_method', 'experiment', 'model_id', and
'parent_experiment_rip'. For filename, relevant parts of
the filenames are used (using regular expressions for
e.g. the ensemble part of the filename). If information is
missing, an exception is raised.
Both options can be given in a list if the first option
fails, the next one is tried. If both fail, an exception
is raised.

The default is [attributes, filename]: first by
attributes, and use the filename as a fallback.
info_from: [attributes, filename]
"""
import logging
import re
import pandas as pd


ATTRIBUTES_KEYS = {
    'experiment': ["experiment_id"],
    'model': ["model_id", "source_id"],
    'realization': ["realization", "realization_index"],
    'initialization': ["initialization_method", "initialization_index"],
    'physics': ["physics_version", "physics_index"],
    'prip': ["parent_experiment_rip", "parent_variant_label"],
    'var':  ["variable_id"],
}
ATTRIBUTES_DEFAULT = {
    'experiment': "",
    'model': "",
    'realization': 0,
    'initialization': 0,
    'physics': 0,
    'prip': None,
    'var': '',
}
FILENAME_PATTERNS = {
    'esmval': """^CMIP\\d_\
(?P<model>[-A-Za-z0-9]+)_\
(?P<mip>[A-Za-z]+)_\
(?P<experiment>[A-Za-z0-9]+)_\
r(?P<realization>\\d+)\
i(?P<initialization>\\d+)\
p(?P<physics>\\d+)_\
(?P<var>[a-z]+)_\
.*\\.nc$\
""",
    'cmip5': """^\
(?P<var>[a-z]+)_\
(?P<mip>[A-Za-z]+)_\
(?P<model>[-A-Za-z0-9]+)_\
(?P<experiment>[A-Za-z0-9]+)_\
r(?P<realization>\\d+)\
i(?P<initialization>\\d+)\
p(?P<physics>\\d+)_\
.*\\.nc$\
""",
    'cmip6': """^\
(?P<var>[a-z]+)_\
(?P<mip>[A-Za-z]+)_\
(?P<model>[-A-Za-z0-9]+)_\
(?P<experiment>[A-Za-z0-9]+)_\
r(?P<realization>\\d+)\
i(?P<initialization>\\d+)\
p(?P<physics>\\d+)\
f\\d+_\
gn_\
.*\\.nc$\
""",
}

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


def _get_attributes_from_cube(attrs, attributes):
    """DUMMY"""
    data, found = {}, True

    data['experiment'] = None
    keys = attributes['experiment']
    if isinstance(keys, str):
        keys = [keys]
    for key in keys:
        if key in attrs:
            data['experiment'] = attrs[key]
            break

    data['model'] = None
    keys = attributes['model']
    if isinstance(keys, str):
        keys = [keys]
    for key in keys:
        if key in attrs:
            data['model'] = attrs[key]
            break

    keys = attributes['var']
    if isinstance(keys, str):
        keys = [keys]
    for key in keys:
        if key in attrs:
            data['var'] = attrs[key]
            break

    data['prip'] = None
    keys = attributes['prip']
    if isinstance(keys, str):
        keys = [keys]
    for key in keys:
        if key in attrs:
            prip = attrs[key]
            match = re.search(r'r(\d+)i(\d+)p(\d+)(f\d+)?', prip)
            if match:
                data['prip'] = int(match.group(1)), int(match.group(2)), int(match.group(3))
                break

    for key in ['realization', 'physics', 'initialization']:
        data[key] = 0
        keys = attributes[key]
        if isinstance(keys, str):
            keys = [keys]
        for k in keys:
            if k in attrs:
                data[key] = int(attrs[k])
                break

    return data, found


def _get_attributes_from_filename(filename, patterns):
    """DUMMY"""
    data, found = {}, True

    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            break
    else:  # no match
        return {}, False

    data['experiment'] = match.group('experiment')
    data['model'] = match.group('model')
    data['var'] = match.group('var')
    for key in ['realization', 'physics', 'initialization']:
        try:
            data[key] = int(match.group(key))
        except IndexError:
            pass

    return data, found


def get_single(cube, path, info_from, attributes, patterns):
    """Obtain attributes for a single cube"""
    attrs = cube.attributes
    filename = path.name
    data = {}
    # Run in reversed, so the information from the most valuable asset
    # is added last (in `data.update(result)`)
    for info in reversed(info_from):
        if info == 'attributes':
            result, _ = _get_attributes_from_cube(attrs, attributes)
        elif info == 'filename':
            result, _ = _get_attributes_from_filename(filename, patterns)
        else:
            raise ValueError("incorrect argument for 'info_from': "
                             "should be 'attributes' or 'filename'")
        data.update(result)

    if 'experiment' not in data:
        raise KeyError("missing experiment info")
    if 'model' not in data:
        raise KeyError("missing model info")
    if not ('realization' in data and 'physics' in data and 'initialization' in data):
        raise KeyError("missing ensemble RIP info")

    return data


def get(cubes, paths, info_from=('attributes', 'filename'),
        attributes=None, filename_pattern=None):
    """Match dataset experiments, historical with future experiments

    Anything where the experiment does not match 'historical' is
    assumed to be a future experiment (rcp, ssp).

    Parameters
    ----------

    - match_by: 'ensemble' or 'model'

    - info_from: 'attributes' or 'filename', or both in a list or tuple

    - attributes: dict, or None

        The attribute names for the following info. This information
        is given as a dict, with the keys below. The default values
        are given after each key, and is used when attributes=None.
        - experiment: "experiment"
        - model: "model_id"
        - realization: "realization"
        - initialization: "initialization"
        - physics: "physics_version"

    """

    if attributes is None:
        attributes = ATTRIBUTES_KEYS
    if filename_pattern is None:
        filename_pattern = FILENAME_PATTERNS.values()
    elif isinstance(filename_pattern, str):
        filename_pattern = [filename_pattern]
    if isinstance(info_from, str):
        info_from = [info_from]

    dataset = []
    for cube, path in zip(cubes, paths):
        default = ATTRIBUTES_DEFAULT.copy()
        attrs = get_single(cube, path, info_from=info_from,
                           attributes=attributes, patterns=filename_pattern)
        default.update(attrs)
        default.update({'cube': cube, 'path': path})
        dataset.append(default)

    dataset = pd.DataFrame(dataset)
    dataset['experiment'] = dataset['experiment'].str.lower()

    return dataset
