"""Run the preprocessor."""

import copy
import logging
from collections import OrderedDict

import iris

from .derive import derive
from .mask import mask
from .multimodel import multi_model_mean
from .reformat import fix_data, fix_file, fix_metadata
from .regrid import vinterp as extract_levels
from .regrid import regrid
from .time_area import area_average as average_region
from .time_area import area_slice as extract_region
from .time_area import time_slice as extract_time
from .time_area import seasonal_mean

logger = logging.getLogger(__name__)

# TODO: review preprocessor functions
SINGLE_MODEL_FUNCTIONS = {
    # File reformatting/CMORization
    'fix_file': fix_file,
    # Load cube from file
    'load': iris.load,
    # Derive variable
    'derive': derive,
    # Metadata reformatting/CMORization
    'fix_metadata': fix_metadata,
    # Time extraction
    'extract_time': extract_time,
    # Data reformatting/CMORization
    'fix_data': fix_data,
    # Level extraction
    'extract_levels': extract_levels,
    # Regridding
    'regrid': regrid,
    # Masking
    'mask': mask,
    # Region selection
    'extract_region': extract_region,
    # Grid-point operations
    'average_region': average_region,
    # 'average_zone': average_zone,
    # 'cross_section': cross_section,
    # Time operations
    # 'annual_cycle': annual_cycle,
    # 'diurnal_cycle': diurnal_cycle,
    'seasonal_mean': seasonal_mean,
    # Save to file
    'save': iris.save,
}

MULTI_MODEL_FUNCTIONS = {
    'load': iris.load,
    'multi_model_mean': multi_model_mean,
    'save': iris.save,
}


def _as_ordered_dict(cfg, order):
    """Copy cfg to an OrderedDict using the specified order."""
    ordered_cfg = OrderedDict()
    for step in order:
        if step in cfg:
            ordered_cfg[step] = cfg[step]

    return ordered_cfg


def preprocess_multi_model(files, cfg):
    """Run multi-model preprocessor functions.

    Parameters
    ----------
    files: iterable over str objects
        Input file paths.
    cfg: dict
        Dictionary with as keys preprocessor functions and as values the
        arguments taken by those functions.

    """

    steps = [
        # TODO: add more steps as they become available
        'load',
        'multi_model_mean',
        'save',
    ]

    cfg = _as_ordered_dict(cfg=cfg, order=steps)
    _preprocess_files(files=files, cfg=cfg, functions=MULTI_MODEL_FUNCTIONS)


def preprocess_single_model(files, cfg):
    """Run single-model preprocessor functions.

    Parameters
    ----------
    files: iterable over str objects
        List of input file paths.
    cfg: dict
        Dictionary with as keys preprocessor functions and as values the
        arguments taken by those functions.

    """
    steps = [
        # TODO: add more steps as they become available
        'fix_file',
        'load',
        'fix_metadata',
        'select_time',
        'fix_data',
        'extract_level',
        'regrid',
        'mask',
        'save',
    ]

    cfg = _as_ordered_dict(cfg=cfg, order=steps)
    _preprocess_files(files=files, cfg=cfg, functions=SINGLE_MODEL_FUNCTIONS)


def preprocess_single_model_derived(files_list, cfg_list):
    """Preprocess a derived variable.

    Parameters
    ----------
    files_list: iterable over iterable over str
        List containing a list of files for each input variable.
    cfg_list: iterable over dicts
        List of dictionaries with as keys preprocessor functions and as values
        the arguments taken by those functions.

    """
    output_files = []
    for i, input_files in enumerate(files_list):
        preprocess_single_model(files=input_files, cfg=cfg_list[i])
        output_files.append(cfg_list[i]['save']['target'])

    preprocess_single_model(files=output_files, cfg=cfg_list[-1])


def _preprocess_files(files, cfg, functions):
    """Run preprocessor"""
    cfg = copy.deepcopy(cfg)

    # Apply file fixes
    if 'fix_file' in cfg:
        args = cfg.pop('fix_file')
        function = functions['fix_file']
        for file in files:
            function(file, **args)

    # Load one or more cubes
    args = cfg.pop('load')
    function = functions['load']
    cubes = function(files, **args)

    # Run preprocessor
    _preprocess_cubes(cubes=cubes, cfg=cfg, functions=functions)


def _preprocess_cubes(cubes, cfg, functions):
    """Apply preprocessor steps defined in `cfg` to `cubes`."""
    cfg = copy.deepcopy(cfg)

    while True:
        try:
            step, args = cfg.popitem(0)
        except KeyError:
            return

        function = functions[step]
        cubes = function(cubes, **args)
