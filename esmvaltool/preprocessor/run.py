"""Run the preprocessor."""

import copy
import logging
from collections import OrderedDict

import iris
from iris.cube import CubeList

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


# TODO: review preprocessor functions
FUNCTIONS = {
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
    'save': iris.save,
}

DEFAULT_ORDER = [
    # TODO: add more steps as they become available
    'fix_file',
    'load',
    'fix_metadata',
    'select_time',
    'fix_data',
    'derive',
    'extract_level',
    'regrid',
    'mask',
    'multi_model_mean',
    'save',
]

_MULTI_CUBE_FUNCTIONS = {
    'derive',
    'multi_model_mean',
}


def _as_ordered_dict(settings, order):
    """Copy settings to an OrderedDict using the specified order."""
    ordered_settings = OrderedDict()
    for step in order:
        if step in settings:
            ordered_settings[step] = settings[step]

    return ordered_settings


def _get_fix_settings(short_name, model, project, mip=None):
    """Get fixes preprocessor configuration"""
    fixcfg = {
        'project': project,
        'model': model,
        'short_name': short_name,
    }

    settings = {}
    settings['fix_file'] = dict(fixcfg)
    if mip:
        fixcfg['mip'] = mip
    settings['fix_metadata'] = dict(fixcfg)
    settings['fix_data'] = dict(fixcfg)

    return settings


def get_default_settings(files, settings, **kwargs):

    cfg = {}

    # Configure loading
    cfg['load'] = {
        'uris': files,
        'method': 'concatenate',
    }
    if 'short_name' in kwargs:
        cfg['load']['constraint'] = kwargs['short_name']

    # Configure fixes
    if any(function.startswith('fix_') for function in settings):
        fixcfg = _get_fix_settings(
            short_name=kwargs['short_name'],
            model=kwargs['model'],
            project=kwargs['project'],
            mip=kwargs.get('mip', None))
        cfg.update(fixcfg)

    # Configure time extraction
    if 'start_year' in kwargs and 'end_year' in kwargs:
        cfg['extract_time'] = {
            'yr1': kwargs['start_year'],
            'yr2': kwargs['end_year'] + 1,
            'mo1': 1,
            'mo2': 1,
            'd1': 1,
            'd2': 1,
        }

    # Override with settings from preprocessor profile
    cfg.update(copy.deepcopy(settings))

    return cfg


def get_multi_model_task(files_list, settings_list, **kwargs):
    """Get a task for preprocessing multiple models"""


def _split_settings(step, settings):
    """Split settings in before and after `step`"""
    index = DEFAULT_ORDER.index(step)
    before_order = list(DEFAULT_ORDER[:index])
    after_order = list(DEFAULT_ORDER[index:])
    before = {(k, v) for k, v in settings if k in before_order}
    after = {(k, v) for k, v in settings if k in after_order}
    return before, after


def get_single_model_task(files, settings, **kwargs):
    """Get a task for preprocessing a single model"""
    settings = get_default_settings(files, settings, **kwargs)

    if 'derive' in settings:
        before_settings, after_settings = _split_settings(
            step='derive', settings=settings)

        # create tasks that should be done before derive step
        before_tasks = []
        for short_name in get_required(kwargs['short_name']):
            before_settings['load']['constraint'] = short_name
            task = PreprocessingTask(settings=before_settings)
            before_tasks.append(task)

        # create complete task
        task = PreprocessingTask(
            settings=after_settings, ancestors=before_tasks)
    else:
        task = PreprocessingTask(settings=settings)

    return task


class PreprocessingTask(AbstractTask):
    """Task for running the preprocessor"""

    def __init__(self, settings, ancestors=None):
        """Initialize"""
        super(PreprocessingTask, self).__init__(settings, ancestors)

        self.order = list(DEFAULT_ORDER)

    def _run(self, input_data):
        settings = _as_ordered_dict(settings=self.settings, order=self.order)
        self.output_data = _preprocess(input_data, settings)


def _preprocess(input_data, settings):
    """Run preprocessor"""
    settings = copy.deepcopy(settings)

    # Apply file fixes
    if 'fix_file' in settings:
        fix_args = settings.pop('fix_file')
        fix = FUNCTIONS['fix_file']
        for file in input_data:
            fix(file, **fix_args)

    if 'load' in settings:
        # Load cubes from file if configured
        load_args = settings.pop('load')
        load = FUNCTIONS['load']
        cubes = load(input_data, **load_args)
    else:
        # Else assume input_data is a CubeList
        cubes = input_data

    # Remove save arguments from settings and store for later use
    save_args = settings.pop('save', None)

    # Run preprocessor
    cubes = _preprocess_cubes(cubes=cubes, settings=settings)

    # Return a CubeList if no save function is configured
    if save_args is None:
        return cubes
    # Else save cubes to file and return the configured filename
    save = FUNCTIONS['save']
    save(cubes, **save_args)
    filename = save_args['target']
    return filename


def _preprocess_cubes(cubes, settings):
    """Apply preprocessor steps defined in `settings` to `cubes`."""
    settings = copy.deepcopy(settings)

    while True:
        try:
            step, args = settings.popitem(0)
        except KeyError:
            return cubes

        function = FUNCTIONS[step]
        if function in _MULTI_CUBE_FUNCTIONS:
            cubes = function(cubes, **args)
        else:
            cubes = CubeList(function(cube, **args) for cube in cubes)
