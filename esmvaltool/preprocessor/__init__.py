"""Preprocessor module."""
import logging
import os
from collections import OrderedDict

from iris.cube import Cube

from ..task import AbstractTask
from ._derive import derive
from ._download import download
from ._io import load_cubes, save_cubes
from ._mask import mask_fillvalues, mask_landocean
from ._multimodel import multi_model_statistics
from ._reformat import fix_data, fix_file, fix_metadata, cmor_check_data
from ._regrid import vinterp as extract_levels
from ._regrid import regrid
from ._time_area import area_average as average_region
from ._time_area import area_slice as extract_region
from ._time_area import time_slice as extract_time
from ._time_area import seasonal_mean

logger = logging.getLogger(__name__)

# TODO: review preprocessor functions
PREPROCESSOR_FUNCTIONS = {
    'download': download,
    # File reformatting/CMORization
    'fix_file': fix_file,
    # Load cube from file
    'load': load_cubes,
    # Metadata reformatting/CMORization
    'fix_metadata': fix_metadata,
    # Time extraction
    'extract_time': extract_time,
    # Data reformatting/CMORization
    'fix_data': fix_data,
    # Derive variable
    'derive': derive,
    # Level extraction
    'extract_levels': extract_levels,
    # Regridding
    'regrid': regrid,
    # Masking
    'mask_landocean': mask_landocean,
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
    'multi_model_statistics': multi_model_statistics,
    'mask_fillvalues': mask_fillvalues,
    'cmor_check_data': cmor_check_data,
    # Save to file
    'save': save_cubes,
}

DEFAULT_ORDER = (
    # TODO: add more steps as they become available
    'download',
    'fix_file',
    'load',
    'derive',
    'fix_metadata',
    'extract_time',
    'fix_data',
    'extract_levels',
    'regrid',
    'mask_landocean',
    'mask_fillvalues',
    'extract_region',
    'average_region',
    'seasonal_mean',
    'multi_model_statistics',
    'cmor_check_data',
    'save',
)

FIXED_DEFAULT_ORDER = list(DEFAULT_ORDER)
assert set(DEFAULT_ORDER) == set(PREPROCESSOR_FUNCTIONS)

MULTI_MODEL_FUNCTIONS = {
    'multi_model_statistics',
    'mask_fillvalues',
}
assert MULTI_MODEL_FUNCTIONS.issubset(set(PREPROCESSOR_FUNCTIONS))

# Preprocessor functions that take a list instead of a file/Cube as input.
_LIST_INPUT_FUNCTIONS = {
    'download',
    'load',
    'derive',
    'mask_fillvalues',
    'multi_model_statistics',
    'save',
}
assert _LIST_INPUT_FUNCTIONS.issubset(set(PREPROCESSOR_FUNCTIONS))

# Preprocessor functions that return a list instead of a file/Cube.
_LIST_OUTPUT_FUNCTIONS = {
    'download',
    'load',
    'mask_fillvalues',
    'multi_model_statistics',
    'save',
}
assert _LIST_OUTPUT_FUNCTIONS.issubset(set(PREPROCESSOR_FUNCTIONS))


def _split_settings(settings, step):
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


def _as_ordered_dict(settings, order):
    """Copy settings to an OrderedDict using the specified order."""
    ordered_settings = OrderedDict()
    for step in order:
        if step in settings:
            ordered_settings[step] = settings[step]

    return ordered_settings


def _group_input(in_files, out_files):
    """Group a list of input files by output file."""
    grouped_files = {}

    def get_matching(in_file):
        """Find the output file which matches input file best."""
        in_chunks = os.path.basename(in_file).split('_')
        score = 0
        fname = None
        for out_file in out_files:
            out_chunks = os.path.basename(out_file).split('_')
            tmp = sum(c in out_chunks for c in in_chunks)
            if tmp > score:
                score = tmp
                fname = out_file
        if not fname:
            logger.warning(
                "Unable to find matching output file for input file %s",
                in_file)
        return fname

    # Group input files by output file
    for in_file in in_files:
        out_file = get_matching(in_file)
        if out_file:
            if out_file not in grouped_files:
                grouped_files[out_file] = []
            grouped_files[out_file].append(in_file)

    return grouped_files


def preprocess_multi_model(input_files, all_settings, order, debug=False):
    """Run preprocessor on multiple models for a single variable."""
    # Group input files by output file
    all_items = _group_input(input_files, all_settings)
    logger.debug("Processing %s", all_items)

    # Define multi model steps
    # This assumes that the multi model settings are the same for all models
    multi_model_steps = [
        step
        for step in DEFAULT_ORDER if (step in MULTI_MODEL_FUNCTIONS and any(
            step in settings for settings in all_settings.values()))
    ]
    final_step = object()
    multi_model_steps.append(final_step)

    # Process
    for step in multi_model_steps:
        multi_model_settings = _get_multi_model_settings(all_settings, step)
        # Run single model steps
        for name in all_settings:
            settings, all_settings[name] = _split_settings(
                all_settings[name], step)
            settings = _as_ordered_dict(settings, order)
            all_items[name] = preprocess(
                items=all_items[name], settings=settings, debug=debug)
        if step is not final_step:
            # Run multi model step (all_items should be cubes by now)
            multi_model_items = [
                cube for name in all_items for cube in all_items[name]
            ]
            all_items = {}
            result = preprocess(
                items=multi_model_items,
                settings=multi_model_settings,
                debug=debug)
            for cube in result:
                filename = cube.attributes['_filename']
                if filename not in all_items:
                    all_items[filename] = []
                all_items[filename].append(cube)

    return [filename for name in all_items for filename in all_items[name]]


def preprocess(items, settings, debug=False):
    """Run preprocessor"""
    if debug:
        step_names = list({setting for setting, arg in settings.items()})
        step_names = sorted(step_names, key=FIXED_DEFAULT_ORDER.index)
    for step, args in settings.items():
        logger.debug("Running preprocessor step %s", step)
        function = PREPROCESSOR_FUNCTIONS[step]

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
            idx = step_names.index(step)
            enum_step = str(idx).zfill(2) + '_' + step
            save_cubes(cubes, debug=debug, step=enum_step)

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
