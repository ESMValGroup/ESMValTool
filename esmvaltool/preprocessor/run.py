"""Run the preprocessor."""
import logging
import os
from collections import OrderedDict

import iris

from esmvaltool.preprocessor.derive import get_required
from esmvaltool.task import AbstractTask

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


def load_cubes(files, mode='merge', **kwargs):
    """Load iris cubes from files"""
    if mode == 'merge':
        cubes = iris.load(files, **kwargs)
        return cubes

    if mode == 'concatenate':
        cubes = iris.load_raw(files, **kwargs)
        cubes.concatenate()
        return cubes

    raise NotImplementedError("mode={} not supported".format(mode))


def save_cubes(cubes, **args):
    """Save iris cubes to file."""
    filename = args['target']
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    iris.save(cubes, **args)
    return filename


# TODO: review preprocessor functions
FUNCTIONS = {
    # File reformatting/CMORization
    'fix_file': fix_file,
    # Load cube from file
    'load': iris.load_cube,
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
    'multi_model_mean': multi_model_mean,
    # Save to file
    'save': save_cubes,
}

DEFAULT_ORDER = (
    # TODO: add more steps as they become available
    'fix_file',
    'load',
    'fix_metadata',
    'extract_time',
    'fix_data',
    'derive',
    'extract_levels',
    'regrid',
    'mask',
    'extract_region',
    'average_region',
    'seasonal_mean',
    'multi_model_mean',
    'save',
)

assert set(DEFAULT_ORDER) == set(FUNCTIONS)

# Preprocessor functions that take a CubeList instead of a Cube as input.
_LIST_INPUT_FUNCTIONS = {
    'load',
    'derive',
    'multi_model_mean',
    'save',
}

assert _LIST_INPUT_FUNCTIONS.issubset(set(FUNCTIONS))


def _as_ordered_dict(settings, order):
    """Copy settings to an OrderedDict using the specified order."""
    ordered_settings = OrderedDict()
    for step in order:
        if step in settings:
            ordered_settings[step] = settings[step]

    return ordered_settings


def _split_settings(step, settings):
    """Split settings in before and after `step`"""
    index = DEFAULT_ORDER.index(step)
    before_order = list(DEFAULT_ORDER[:index])
    after_order = list(DEFAULT_ORDER[index:])
    before = {(k, v) for k, v in settings if k in before_order}
    after = {(k, v) for k, v in settings if k in after_order}
    return before, after


def get_single_model_task(settings, short_name=None):
    """Get a task for preprocessing a single model"""

    if 'derive' in settings:
        if not short_name:
            raise ValueError("Cannot determine input variables.")
        before_settings, after_settings = _split_settings(
            step='derive', settings=settings)

        # create tasks that should be done before derive step
        before_tasks = []
        for input_short_name in get_required(short_name):
            before_settings['load']['constraints'] = input_short_name
            task = PreprocessingTask(settings=before_settings)
            before_tasks.append(task)

        # create complete task
        task = PreprocessingTask(
            settings=after_settings, ancestors=before_tasks)
    else:
        task = PreprocessingTask(settings=settings)

    return task


def get_multi_model_task(settings):
    """Get a task for preprocessing multiple models"""


class PreprocessingTask(AbstractTask):
    """Task for running the preprocessor"""

    def __init__(self, settings, order=DEFAULT_ORDER, ancestors=None):
        """Initialize"""
        super(PreprocessingTask, self).__init__(settings, ancestors)
        self.order = list(order)

    def _run(self, input_data):
        settings = _as_ordered_dict(settings=self.settings, order=self.order)
        # If input_data is not available from ancestors and also not
        # specified in self.run(input_data), try to get it from settings
        if not self.ancestors and not input_data:
            if 'load' in settings and 'uris' in settings['load']:
                input_data = settings['load'].pop('uris')
        self.output_data = preprocess(input_data, settings)

    def __str__(self):
        txt = "{}:\norder: {}\n{}".format(
            self.__class__.__name__,
            self.order,
            super(PreprocessingTask, self).str(),
        )
        return txt


def preprocess(items, settings):
    """Run preprocessor"""
    for step, args in settings.items():
        function = FUNCTIONS[step]
        if step in _LIST_INPUT_FUNCTIONS:
            items = (items, )
        items = tuple(function(item, **args) for item in items)

    return items
