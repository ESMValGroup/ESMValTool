"""Run the preprocessor."""
import logging
import os
from collections import OrderedDict

import iris

from ..task import AbstractTask
from ._derive import derive, get_required
from ._download import download
from ._mask import mask_fillvalues, mask_landocean
from ._multimodel import multi_model_mean
from ._reformat import fix_data, fix_file, fix_metadata
from ._regrid import vinterp as extract_levels
from ._regrid import regrid
from ._time_area import area_average as average_region
from ._time_area import area_slice as extract_region
from ._time_area import time_slice as extract_time
from ._time_area import seasonal_mean

iris.FUTURE.netcdf_promote = True
iris.FUTURE.netcdf_no_unlimited = True

logger = logging.getLogger(__name__)


def load_cubes(files, mode='merge', **kwargs):
    """Load iris cubes from files"""
    if mode == 'ordered':
        # load cubes in order of files
        cubes = []
        for path in files:
            logger.debug("Loading %s", path)
            cube = iris.load_cube(path, **kwargs)
            cubes.append(cube)
        return cubes

    if mode == 'merge':
        logger.debug("Loading and merging:\n%s", "\n".join(files))
        cubes = iris.load(files, **kwargs)
        return cubes

    if mode == 'concatenate':
        logger.debug("Loading and concatenating:\n%s", "\n".join(files))
        cubes = iris.load_raw(files, **kwargs)
        iris.util.unify_time_units(cubes)
        cubes = cubes.concatenate()
        return cubes

    raise NotImplementedError("mode={} not supported".format(mode))


def _save_cubes(cubes, **args):
    """Save iris cube to file."""
    filename = args['target']

    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    if (os.path.exists(filename)
            and all(cube.has_lazy_data() for cube in cubes)):
        logger.info("Not saving cubes %s to %s to avoid data loss. "
                    "The cube is probably unchanged.", cubes, filename)
    else:
        logger.debug("Saving cubes %s to %s", cubes, filename)
        iris.save(cubes, **args)

    return filename


def save_cubes(cubes, **args):
    """Save iris cubes to file."""
    filenames = []

    if 'target' in args:
        filename = _save_cubes(cubes, **args)
        filenames.append(filename)
    elif 'targets' in args:
        n_targets, n_cubes = len(args['targets']), len(cubes)
        if n_targets < n_cubes:
            logger.error("Only %s targets specified  for %s cubes:\n"
                         "targets: %s\n"
                         "cubes: %s", n_targets, n_cubes, args['targets'],
                         cubes)
            raise ValueError("Unable to save all cubes.")

        for cube, target in zip(cubes, args.pop('targets')):
            filename = _save_cubes([cube], target=target, **args)
            filenames.append(filename)
    else:
        raise ValueError("No target(s) specified.")

    return filenames


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
    'multi_model_mean': multi_model_mean,
    'mask_fillvalues': mask_fillvalues,
    # Save to file
    'save': save_cubes,
}

DEFAULT_ORDER = (
    # TODO: add more steps as they become available
    'download',
    'fix_file',
    'load',
    'fix_metadata',
    'extract_time',
    'fix_data',
    'derive',
    'extract_levels',
    'regrid',
    'mask_landocean',
    'mask_fillvalues',
    'extract_region',
    'average_region',
    'seasonal_mean',
    'multi_model_mean',
    'save',
)
assert set(DEFAULT_ORDER) == set(PREPROCESSOR_FUNCTIONS)

MULTI_MODEL_FUNCTIONS = {
    'multi_model_mean',
    'mask_fillvalues',
}
assert MULTI_MODEL_FUNCTIONS.issubset(set(PREPROCESSOR_FUNCTIONS))

# Preprocessor functions that take a list instead of a file/Cube as input.
_LIST_INPUT_FUNCTIONS = {
    'download',
    'load',
    'derive',
    'mask_fillvalues',
    'multi_model_mean',
    'save',
}
assert _LIST_INPUT_FUNCTIONS.issubset(set(PREPROCESSOR_FUNCTIONS))

# Preprocessor functions that return a list instead of a file/Cube.
_LIST_OUTPUT_FUNCTIONS = {
    'download',
    'load',
    'mask_fillvalues',
    'save',
}
assert _LIST_OUTPUT_FUNCTIONS.issubset(set(PREPROCESSOR_FUNCTIONS))


def select_single_model_settings(settings):
    """Select settings that belong to a single model function."""
    settings = {
        step: settings[step]
        for step in settings if step not in MULTI_MODEL_FUNCTIONS
    }

    return settings


def select_multi_model_settings(settings):
    """Select settings that belong to a multi model function."""
    settings = {
        step: settings[step]
        for step in settings if step in MULTI_MODEL_FUNCTIONS
    }

    return settings


def _as_ordered_dict(settings, order):
    """Copy settings to an OrderedDict using the specified order."""
    ordered_settings = OrderedDict()
    for step in order:
        if step in settings:
            ordered_settings[step] = settings[step]

    return ordered_settings


# TODO: remove the two functions below once derived variables implemented
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


class PreprocessingTask(AbstractTask):
    """Task for running the preprocessor"""

    def __init__(self,
                 settings,
                 order=DEFAULT_ORDER,
                 ancestors=None,
                 input_data=None,
                 debug=None):
        """Initialize"""
        super(PreprocessingTask, self).__init__(settings, ancestors)
        self.order = list(order)
        self.debug = debug
        self._input_data = input_data

    def _run(self, input_data):
        settings = _as_ordered_dict(settings=self.settings, order=self.order)
        # If input_data is not available from ancestors and also not
        # specified in self.run(input_data), use default
        if not self.ancestors and not input_data:
            input_data = self._input_data
        output_data = preprocess(input_data, settings, debug=self.debug) 
        return output_data

    def __str__(self):
        """Get human readable description."""
        txt = "{}:\norder: {}\n{}".format(
            self.__class__.__name__,
            self.order,
            super(PreprocessingTask, self).str(),
        )
        if self._input_data is not None:
            txt += '\ninput_data: {}'.format(self._input_data)
        return txt


def preprocess(items, settings, debug=None):
    """Run preprocessor"""
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
            cubes = [
                item for item in items if isinstance(item, iris.cube.Cube)
            ]
            if cubes:
                targets = []
                for path in debug['paths']:
                    target = os.path.join(path, step + '.nc')
                    if os.path.exists(target):
                        target = os.path.splitext(target)[0] + '_1.nc'
                    targets.append(target)
                save_cubes(cubes, targets=targets)

    return items
