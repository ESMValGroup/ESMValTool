"""Preprocessor module."""
import logging
import os

from iris.cube import Cube

from .._task import AbstractTask
from ._derive import derive
from ._download import download
from ._io import cleanup, extract_metadata, load_cubes, save_cubes, concatenate
from ._mask import mask_fillvalues, mask_landocean
from ._multimodel import multi_model_statistics
from ._reformat import fix_data, fix_file, fix_metadata, cmor_check_data, \
    cmor_check_metadata
from ._regrid import vinterp as extract_levels
from ._regrid import regrid
from ._time_area import area_average as average_region
from ._time_area import area_slice as extract_region
from ._time_area import time_slice as extract_time
from ._time_area import seasonal_mean

logger = logging.getLogger(__name__)

__all__ = [
    'download',
    # File reformatting/CMORization
    'fix_file',
    # Load cube from file
    'load_cubes',
    # Derive variable
    'derive',
    # Metadata reformatting/CMORization
    'fix_metadata',
    # Concatenate all cubes in one
    'concatenate',
    'cmor_check_metadata',
    # Time extraction
    'extract_time',
    # Data reformatting/CMORization
    'fix_data',
    # Level extraction
    'extract_levels',
    # Regridding
    'regrid',
    # Masking
    'mask_landocean',
    'mask_fillvalues',
    # Region selection
    'extract_region',
    # Grid-point operations
    'average_region',
    # 'average_zone': average_zone,
    # 'cross_section': cross_section,
    # Time operations
    # 'annual_cycle': annual_cycle,
    # 'diurnal_cycle': diurnal_cycle,
    'seasonal_mean',
    'multi_model_statistics',
    'cmor_check_data',
    # Save to file
    'save_cubes',
    'cleanup',
    'extract_metadata',
]

DEFAULT_ORDER = tuple(__all__)
assert set(DEFAULT_ORDER).issubset(set(globals()))

MULTI_MODEL_FUNCTIONS = {
    'multi_model_statistics',
    'mask_fillvalues',
    'extract_metadata',
}
assert MULTI_MODEL_FUNCTIONS.issubset(set(DEFAULT_ORDER))

# Preprocessor functions that take a list instead of a file/Cube as input.
_LIST_INPUT_FUNCTIONS = MULTI_MODEL_FUNCTIONS | {
    'download',
    'load_cubes',
    'concatenate',
    'derive',
    'save_cubes',
    'cleanup',
}
assert _LIST_INPUT_FUNCTIONS.issubset(set(DEFAULT_ORDER))

# Preprocessor functions that return a list instead of a file/Cube.
_LIST_OUTPUT_FUNCTIONS = MULTI_MODEL_FUNCTIONS | {
    'download',
    'load_cubes',
    'save_cubes',
    'cleanup',
}
assert _LIST_OUTPUT_FUNCTIONS.issubset(set(DEFAULT_ORDER))


def split_settings(settings, step):
    """Split settings, using step as a separator."""
    before = {}
    for _step in DEFAULT_ORDER:
        if _step == step:
            break
        if _step in settings:
            before[_step] = settings[_step]
    after = {
        k: v
        for k, v in settings.items() if not (k == step or k in before)
    }
    return before, after


def _get_multi_model_settings(all_settings, step):
    """Select settings for multi model step"""
    for settings in all_settings.values():
        if step in settings:
            return {step: settings[step]}


def _group_input(in_files, out_files):
    """Group a list of input files by output file."""
    grouped_files = {}

    def get_matching(in_file):
        """Find the output file which matches input file best."""
        in_chunks = os.path.basename(in_file).split('_')
        score = 0
        fname = []
        for out_file in out_files:
            out_chunks = os.path.basename(out_file).split('_')
            tmp = sum(c in out_chunks for c in in_chunks)
            if tmp > score:
                score = tmp
                fname = [out_file]
            elif tmp == score:
                fname.append(out_file)
        if not fname:
            logger.warning(
                "Unable to find matching output file for input file %s",
                in_file)
        return fname

    # Group input files by output file
    for in_file in in_files:
        for out_file in get_matching(in_file):
            if out_file not in grouped_files:
                grouped_files[out_file] = []
            grouped_files[out_file].append(in_file)

    return grouped_files


def preprocess_multi_model(input_files, all_settings, order, debug=False):
    """Run preprocessor on multiple models for a single variable."""
    # Group input files by output file
    all_items = _group_input(input_files, all_settings)
    logger.debug("Processing %s", all_items)

    # List of all preprocessor steps used
    steps = [
        step for step in order
        if any(step in settings for settings in all_settings.values())
    ]
    # Find multi model steps
    # This assumes that the multi model settings are the same for all models
    multi_model_steps = [
        step for step in steps if step in MULTI_MODEL_FUNCTIONS
    ]
    # Append a dummy multi model step if the final step is not multi model
    dummy_step = object()
    if steps[-1] not in MULTI_MODEL_FUNCTIONS:
        multi_model_steps.append(dummy_step)

    # Process
    for step in multi_model_steps:
        multi_model_settings = _get_multi_model_settings(all_settings, step)
        # Run single model steps
        for name in all_settings:
            settings, all_settings[name] = split_settings(
                all_settings[name], step)
            all_items[name] = preprocess(all_items[name], settings, order,
                                         debug)
        if step is not dummy_step:
            # Run multi model step
            multi_model_items = [
                item for name in all_items for item in all_items[name]
            ]
            all_items = {}
            result = preprocess(multi_model_items, multi_model_settings, order,
                                debug)
            for item in result:
                if isinstance(item, Cube):
                    name = item.attributes['_filename']
                    if name not in all_items:
                        all_items[name] = []
                    all_items[name].append(item)
                else:
                    all_items[item] = [item]

    return [filename for name in all_items for filename in all_items[name]]


def preprocess(items, settings, order, debug=False):
    """Run preprocessor"""
    steps = (step for step in order if step in settings)
    for step in steps:
        logger.debug("Running preprocessor step %s", step)
        function = globals()[step]
        args = settings[step]

        if step in _LIST_INPUT_FUNCTIONS:
            logger.debug("Running %s(%s, %s)", function.__name__, items, args)
            result = [function(items, **args)]
        else:
            result = []
            for item in items:
                logger.debug("Running %s(%s, %s)", function.__name__, item,
                             args)
                result.append(function(item, **args))

        if step in _LIST_OUTPUT_FUNCTIONS:
            items = tuple(item for subitem in result for item in subitem)
        else:
            items = tuple(result)

        if debug:
            logger.debug("Result %s", items)
            cubes = [item for item in items if isinstance(item, Cube)]
            save_cubes(cubes, debug=debug, step=step)

    return items


class PreprocessingTask(AbstractTask):
    """Task for running the preprocessor"""

    def __init__(self,
                 settings,
                 output_dir,
                 ancestors=None,
                 input_files=None,
                 order=DEFAULT_ORDER,
                 debug=None):
        """Initialize"""
        super(PreprocessingTask, self).__init__(
            settings=settings, output_dir=output_dir, ancestors=ancestors)
        self.order = list(order)
        self.debug = debug
        self._input_files = input_files

    def _run(self, input_files):
        # If input_data is not available from ancestors and also not
        # specified in self.run(input_data), use default
        if not self.ancestors and not input_files:
            input_files = self._input_files
        output_files = preprocess_multi_model(
            input_files, self.settings, self.order, debug=self.debug)
        return output_files

    def __str__(self):
        """Get human readable description."""
        settings = dict(self.settings)
        self.settings = {
            os.path.basename(k): v
            for k, v in self.settings.items()
        }

        txt = "{}:\norder: {}\n{}".format(
            self.__class__.__name__,
            self.order,
            super(PreprocessingTask, self).str(),
        )

        self.settings = settings

        if self._input_files is not None:
            txt += '\ninput_files: {}'.format(self._input_files)
        return txt
