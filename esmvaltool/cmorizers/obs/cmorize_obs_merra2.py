"""ESMValTool CMORizer for MERRA2 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/

Last access
    20191129

Download and processing instructions
    - For download instructions see the download script `download_merra2.sh`.

"""
import logging
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime, timedelta
from os import cpu_count
from pathlib import Path
from warnings import catch_warnings, filterwarnings

import glob
import iris
import os
import numpy as np

from esmvalcore.cmor.table import CMOR_TABLES
from esmvalcore.preprocessor import daily_statistics

from . import utilities as utils

logger = logging.getLogger(__name__)




def _extract_variable(in_files, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input files '%s'",
                var['short_name'], ', '.join(in_files))
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

#    cube = _load_cube(in_files, var)
    var_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var['raw']))
    cube_list = iris.load_raw(in_files)#,constraint=var_constraint)
    selected = [c for c in cube_list if c.var_name == var['raw']]
    selected = iris.cube.CubeList(selected)

    drop_attrs = ['History', 'Filename', 'Comment', 'RangeBeginningDate', 'RangeEndingDate', 'GranuleID', 'ProductionDateTime']
    for c in selected:
        for attr in drop_attrs:
            c.attributes.pop(attr)
    # Needed before merging
    from iris.util import unify_time_units
    unify_time_units(selected)
    # TODO remove strange attributes to time coordinates 
    # selected[1].coord('time')

    import IPython;IPython.embed()
    selected.concatenate_cube()



    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    _fix_units(cube, definition)

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    cube = _fix_coordinates(cube, definition)

    logger.debug("Saving cube\n%s", cube)
    utils.save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        local_keys=['positive'],
    )
    logger.info("Finished CMORizing %s", ', '.join(in_files))



def cmorization(in_dir, out_dir, cfg, config_user):
    """Run CMORizer for MERRA-2."""
    cfg['attributes']['comment'] = cfg['attributes']['comment'].strip().format(
        year=datetime.now().year)
    cfg.pop('cmor_table')

    year = '2018'

    for short_name, var in cfg['variables'].items():
        if 'short_name' not in var:
            var['short_name'] = short_name
        # Now get list of files
        filepattern = os.path.join(in_dir, var['file'].format(year=year))
        in_files = glob.glob(filepattern)
        _extract_variable(in_files, var, cfg, out_dir)

