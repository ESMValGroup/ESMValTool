"""Preprocessor module."""
import copy
import inspect
import logging

import six
from iris.cube import Cube, CubeList

from .._provenance import TrackedFile
from .._task import BaseTask
from ._area import (average_region, extract_named_regions, extract_region,
                    zonal_means)
from ._derive import derive
from ._download import download
from ._io import (_get_debug_filename, cleanup, concatenate, load, save,
                  write_metadata)
from ._mask import (mask_above_threshold, mask_below_threshold,
                    mask_fillvalues, mask_inside_range, mask_landsea,
                    mask_landseaice, mask_outside_range)
from ._multimodel import multi_model_statistics
from ._reformat import (cmor_check_data, cmor_check_metadata, fix_data,
                        fix_file, fix_metadata)
from ._regrid import extract_levels, regrid
from ._time import (annual_mean, extract_month, extract_season, extract_time,
                    regrid_time, seasonal_mean, time_average)
from ._volume import (average_volume, depth_integration, extract_trajectory,
                      extract_transect, extract_volume)

logger = logging.getLogger(__name__)

__all__ = [
    'download',
    # File reformatting/CMORization
    'fix_file',
    # Load cubes from file
    'load',
    # Derive variable
    'derive',
    # Metadata reformatting/CMORization
    'fix_metadata',
    # Concatenate all cubes in one
    'concatenate',
    'cmor_check_metadata',
    # Time extraction
    'extract_time',
    'extract_season',
    'extract_month',
    # Data reformatting/CMORization
    'fix_data',
    # Level extraction
    'extract_levels',
    # Mask landsea (fx or Natural Earth)
    'mask_landsea',
    # Mask landseaice, sftgif only
    'mask_landseaice',
    # Regridding
    'regrid',
    # Masking missing values
    'mask_fillvalues',
    'mask_above_threshold',
    'mask_below_threshold',
    'mask_inside_range',
    'mask_outside_range',
    # Region selection
    'extract_region',
    'extract_volume',
    'extract_trajectory',
    'extract_transect',
    # 'average_zone': average_zone,
    # 'cross_section': cross_section,
    'multi_model_statistics',
    # Grid-point operations
    'extract_named_regions',
    'depth_integration',
    'average_region',
    'average_volume',
    # Time operations
    # 'annual_cycle': annual_cycle,
    # 'diurnal_cycle': diurnal_cycle,
    'zonal_means',
    'annual_mean',
    'seasonal_mean',
    'time_average',
    'regrid_time',
    'cmor_check_data',
    # Save to file
    'save',
    'cleanup',
]

DEFAULT_ORDER = tuple(__all__)

# The order of intial and final steps cannot be configured
INITIAL_STEPS = DEFAULT_ORDER[:DEFAULT_ORDER.index('fix_data') + 1]
FINAL_STEPS = DEFAULT_ORDER[DEFAULT_ORDER.index('cmor_check_data'):]

MULTI_MODEL_FUNCTIONS = {
    'multi_model_statistics',
    'mask_fillvalues',
}


def _get_itype(step):
    """Get the input type of a preprocessor function."""
    function = globals()[step]
    itype = inspect.getargspec(function).args[0]
    return itype


def check_preprocessor_settings(settings):
    """Check preprocessor settings."""
    # The inspect functions getargspec and getcallargs are deprecated
    # in Python 3, but their replacements are not available in Python 2.
    # TODO: Use the new Python 3 inspect API
    for step in settings:
        if step not in DEFAULT_ORDER:
            raise ValueError(
                "Unknown preprocessor function '{}', choose from: {}".format(
                    step, ', '.join(DEFAULT_ORDER)))

        function = function = globals()[step]
        argspec = inspect.getargspec(function)
        args = argspec.args[1:]
        # Check for invalid arguments
        invalid_args = set(settings[step]) - set(args)
        if invalid_args:
            raise ValueError(
                "Invalid argument(s): {} encountered for preprocessor "
                "function {}. \nValid arguments are: [{}]".format(
                    ', '.join(invalid_args), step, ', '.join(args)))

        # Check for missing arguments
        defaults = argspec.defaults
        end = None if defaults is None else -len(defaults)
        missing_args = set(args[:end]) - set(settings[step])
        if missing_args:
            raise ValueError(
                "Missing required argument(s) {} for preprocessor "
                "function {}".format(missing_args, step))
        # Final sanity check in case the above fails to catch a mistake
        try:
            inspect.getcallargs(function, None, **settings[step])
        except TypeError:
            logger.error(
                "Wrong preprocessor function arguments in "
                "function '%s'", step)
            raise


def _check_multi_model_settings(products):
    """Check that multi dataset settings are identical for all products."""
    multi_model_steps = (step for step in MULTI_MODEL_FUNCTIONS
                         if any(step in p.settings for p in products))
    for step in multi_model_steps:
        reference = None
        for product in products:
            settings = product.settings.get(step)
            if settings is None:
                continue
            elif reference is None:
                reference = product
            elif reference.settings[step] != settings:
                raise ValueError(
                    "Unable to combine differing multi-dataset settings for "
                    "{} and {}, {} and {}".format(
                        reference.filename, product.filename,
                        reference.settings[step], settings))


def _get_multi_model_settings(products, step):
    """Select settings for multi model step"""
    _check_multi_model_settings(products)
    settings = {}
    exclude = set()
    for product in products:
        if step in product.settings:
            settings = product.settings[step]
        else:
            exclude.add(product)
    return settings, exclude


def _run_preproc_function(function, items, kwargs):
    """Run preprocessor function."""
    msg = "{}({}, {})".format(function.__name__, items, kwargs)
    logger.debug("Running %s", msg)
    try:
        return function(items, **kwargs)
    except Exception:
        logger.error("Failed to run %s", msg)
        raise


def preprocess(items, step, **settings):
    """Run preprocessor"""
    logger.debug("Running preprocessor step %s", step)
    function = globals()[step]
    itype = _get_itype(step)

    result = []
    if itype.endswith('s'):
        result.append(_run_preproc_function(function, items, settings))
    else:
        for item in items:
            result.append(_run_preproc_function(function, item, settings))

    items = []
    for item in result:
        if isinstance(item,
                      (PreprocessorFile, Cube, six.string_types)):
            items.append(item)
        else:
            items.extend(item)

    return items


def get_step_blocks(steps, order):
    """Group steps into execution blocks."""
    blocks = []
    prev_step_type = None
    for step in order[order.index('load') + 1:order.index('save')]:
        if step in steps:
            step_type = step in MULTI_MODEL_FUNCTIONS
            if step_type is not prev_step_type:
                block = []
                blocks.append(block)
            prev_step_type = step_type
            block.append(step)
    return blocks


class PreprocessorFile(TrackedFile):
    """Preprocessor output file."""

    def __init__(self, attributes, settings, ancestors=None):
        super(PreprocessorFile, self).__init__(attributes['filename'],
                                               attributes, ancestors)

        self.settings = copy.deepcopy(settings)
        if 'save' not in self.settings:
            self.settings['save'] = {}
        self.settings['save']['filename'] = self.filename

        self.files = [a.filename for a in ancestors or ()]

        self._cubes = None
        self._prepared = False

    def check(self):
        """Check preprocessor settings."""
        check_preprocessor_settings(self.settings)

    def apply(self, step, debug=False):
        """Apply preprocessor step to product."""
        if step not in self.settings:
            raise ValueError(
                "PreprocessorFile {} has no settings for step {}".format(
                    self, step))
        self.cubes = preprocess(self.cubes, step, **self.settings[step])
        if debug:
            logger.debug("Result %s", self.cubes)
            filename = _get_debug_filename(self.filename, step)
            save(self.cubes, filename)

    def prepare(self):
        """Apply preliminary file operations on product."""
        if not self._prepared:
            for step in DEFAULT_ORDER[:DEFAULT_ORDER.index('load')]:
                if step in self.settings:
                    self.files = preprocess(self.files, step,
                                            **self.settings[step])
            self._prepared = True

    @property
    def cubes(self):
        """Cubes."""
        if self.is_closed:
            self.prepare()
            self._cubes = preprocess(self.files, 'load',
                                     **self.settings.get('load', {}))
        return self._cubes

    @cubes.setter
    def cubes(self, value):
        self._cubes = value

    def save(self):
        """Save cubes to disk."""
        if self._cubes is not None:
            self.files = preprocess(self._cubes, 'save',
                                    **self.settings['save'])
            self.files = preprocess(self.files, 'cleanup',
                                    **self.settings.get('cleanup', {}))

    def close(self):
        """Close the file."""
        self.save()
        self._cubes = None

    @property
    def is_closed(self):
        """Check if the file is closed."""
        return self._cubes is None

    def _initialize_entity(self):
        """Initialize the entity representing the file."""
        super(PreprocessorFile, self)._initialize_entity()
        settings = {
            'preprocessor:' + k: str(v)
            for k, v in self.settings.items()
        }
        self.entity.add_attributes(settings)


# TODO: use a custom ProductSet that raises an exception if you try to
# add the same Product twice


def _apply_multimodel(products, step, debug):
    """Apply multi model step to products."""
    settings, exclude = _get_multi_model_settings(products, step)

    logger.debug("Applying %s to\n%s", step, '\n'.join(
        str(p) for p in products - exclude))
    result = preprocess(products - exclude, step, **settings)
    products = set(result) | exclude

    if debug:
        for product in products:
            logger.debug("Result %s", product.filename)
            if not product.is_closed:
                for cube in product.cubes:
                    logger.debug("with cube %s", cube)

    return products


class PreprocessingTask(BaseTask):
    """Task for running the preprocessor"""

    def __init__(
            self,
            products,
            ancestors=None,
            name='',
            order=DEFAULT_ORDER,
            debug=None,
            write_ncl_interface=False,
    ):
        """Initialize"""
        super(PreprocessingTask, self).__init__(ancestors=ancestors, name=name)
        _check_multi_model_settings(products)
        self.products = set(products)
        self.order = list(order)
        self.debug = debug
        self.write_ncl_interface = write_ncl_interface

    def _intialize_product_provenance(self):
        """Initialize product provenance."""
        for product in self.products:
            product.initialize_provenance(self.activity)

        # Hacky way to initialize the multi model products as well.
        step = 'multi_model_statistics'
        input_products = [p for p in self.products if step in p.settings]
        if input_products:
            statistic_products = input_products[0].settings[step].get(
                'output_products', {}).values()
            for product in statistic_products:
                product.initialize_provenance(self.activity)

    def _run(self, _):
        """Run the preprocessor."""
        self._intialize_product_provenance()

        steps = {
            step
            for product in self.products for step in product.settings
        }
        blocks = get_step_blocks(steps, self.order)
        for block in blocks:
            logger.debug("Running block %s", block)
            if block[0] in MULTI_MODEL_FUNCTIONS:
                for step in block:
                    self.products = _apply_multimodel(self.products, step,
                                                      self.debug)
            else:
                for product in self.products:
                    logger.debug("Applying single-model steps to %s", product)
                    for step in block:
                        if step in product.settings:
                            product.apply(step, self.debug)
                    if block == blocks[-1]:
                        product.close()

        for product in self.products:
            product.close()
        metadata_files = write_metadata(self.products,
                                        self.write_ncl_interface)
        return metadata_files

    def __str__(self):
        """Get human readable description."""
        order = [
            step for step in self.order
            if any(step in product.settings for product in self.products)
        ]
        products = '\n\n'.join(str(p) for p in self.products)
        txt = "{}:\norder: {}\n{}\n{}".format(
            self.__class__.__name__,
            order,
            products,
            super(PreprocessingTask, self).str(),
        )
        return txt
