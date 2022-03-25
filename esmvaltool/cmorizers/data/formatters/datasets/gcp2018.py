"""ESMValTool CMORizer for GCP2018 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.icos-cp.eu/GCP/2018

Last access
    20210908

Download and processing instructions
    Download the following file: '2018 Global Budget v1.0'
"""

import logging
import os
import warnings
from datetime import datetime

import iris
import numpy as np
import pandas as pd
from cf_units import Unit

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_coords(data_table):
    """Extract coordinates."""
    time_units = Unit('days since 1850-01-01 00:00:00')
    times = [datetime(year=year, month=6, day=15) for year in data_table.index]
    time_dim = iris.coords.DimCoord(time_units.date2num(times),
                                    var_name='time',
                                    standard_name='time',
                                    long_name='time',
                                    units=time_units)
    time_dim.guess_bounds()
    lat_dim = iris.coords.DimCoord([0.0],
                                   bounds=[[-90.0, 90.0]],
                                   var_name='lat',
                                   standard_name='latitude',
                                   long_name='latitude',
                                   units=Unit('degrees_north'))
    lon_dim = iris.coords.DimCoord([180.0],
                                   bounds=[[0.0, 360.0]],
                                   var_name='lon',
                                   standard_name='longitude',
                                   long_name='longitude',
                                   units=Unit('degrees_east'))
    return [(time_dim, 0), (lat_dim, 1), (lon_dim, 2)]


def _extract_variable(variable_name, var, cfg, data_table, out_dir):
    """Extract variable."""
    # Set correct header for data_frame and remove lines that do not include
    # final data (indicated by NaNs)
    header = data_table.iloc[cfg['header_line']]
    data_table = data_table[cfg['header_line'] + 1:]
    data_table.columns = header
    data_table = data_table.dropna()

    # Coordinates
    coords = _get_coords(data_table)

    # Data
    if variable_name == 'fgco2':
        new_data = data_table['ocean sink'].values
    elif variable_name == 'nbp':
        new_data = (data_table['land sink'].values -
                    data_table['land-use change emissions'].values)
    elif variable_name == 'nbp_residual':
        new_data = (
            data_table['fossil emissions excluding carbonation'].values -
            data_table['atmospheric growth'].values -
            data_table['ocean sink'].values -
            data_table['land-use change emissions'].values)
    else:
        raise NotImplementedError(
            f"Derivation of '{variable_name}' not possible yet")
    for _ in range(2):
        new_data = np.expand_dims(new_data, -1)
    new_units = Unit('Gt yr-1')
    if var.get('area'):
        new_data /= var['area']
        new_units = Unit('Gt yr-1 m-2')
    cube = iris.cube.Cube(new_data.astype(np.float32),
                          dim_coords_and_dims=coords,
                          units=new_units)

    # Fix units
    short_name = var.pop('short_name', variable_name)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix metadata
    attrs = dict(cfg['attributes'])
    version_suffix = var.pop('version_suffix', None)
    if version_suffix is not None:
        attrs['version'] += f'-{version_suffix}'
    attrs.update(var)
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    filepath = os.path.join(in_dir, cfg['filename'])
    logger.info("Reading '%s'", filepath)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Unknown extension is not supported',
            category=UserWarning,
            module='openpyxl',
        )
        data_table = pd.read_excel(filepath,
                                   sheet_name='Global Carbon Budget',
                                   index_col=0)

    # Run the cmorization
    for (variable_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", variable_name)
        _extract_variable(variable_name, var, cfg, data_table, out_dir)
