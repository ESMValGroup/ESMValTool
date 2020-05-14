"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
   http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid

Last access
   20190510

Download and processing instructions
   Download and unpack the sithick files from:
   https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/v2.1/hiday/

   And the grid info files from:
   https://pscfiles.apl.washington.edu/zhang/PIOMAS/utilities/grid.dat
   https://pscfiles.apl.washington.edu/zhang/PIOMAS/utilities/grid.dat.pop

   Other variables provided by PIOMAS are not supported, but extending support
   should be achievable for most of them just modifying the config file
"""

import logging
import os
import glob
import collections

import numpy as np
from cf_units import Unit
import iris
from iris.coords import AuxCoord, DimCoord

from esmvaltool.cmorizers.obs.utilities import save_variable, set_global_atts

logger = logging.getLogger(__name__)

# read in CMOR configuration

NX = 360
NY = 120


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']

    cmorizer = PIOMAS(cfg, in_dir, out_dir)
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    cmorizer.prepare_grid_info()
    cmorizer.cmorize()


Coords = collections.namedtuple('Coords', 'lat lon')


class PIOMAS():
    """Cmorizer class for PIOMAS."""

    def __init__(self, cfg, in_dir, out_dir):
        self.cfg = cfg
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.scalar_coords = None
        self.vector_coords = None
        self.areacello = None

    def prepare_grid_info(self):
        """Read grid information."""
        grids = np.loadtxt(os.path.join(
            self.in_dir, self.cfg['custom']['scalar_file']))
        grids = grids.reshape(2, NY, NX)
        self.scalar_coords = self._create_lat_lon_coords(
            grids[1, ...], grids[0, ...])

        grids = np.loadtxt(os.path.join(
            self.in_dir, self.cfg['custom']['vector_file']))
        grids = grids.reshape(7, NY, NX)
        self.vector_coords = self._create_lat_lon_coords(
            grids[1, ...], grids[0, ...])

        # Area in m2
        self.areacello = grids[2, ...] * grids[3, ...] * 1e6

    def cmorize(self):
        """Cmorize available data."""
        # run the cmorization
        for var, vals in self.cfg['variables'].items():
            var_info = self.cfg['cmor_table'].get_variable(vals['mip'], var)
            self.cfg['attributes']['mip'] = vals['mip']

            if vals['type'] == 'scalar':
                coords = self.scalar_coords
            else:
                coords = self.vector_coords
            if var == "areacello":
                cube = self._create_areacello(coords, var_info)
                set_global_atts(cube, self.cfg['attributes'])
                save_variable(
                    cube, var_info.short_name, self.out_dir,
                    self.cfg['attributes'],
                )
            else:
                self._cmorize_var(var_info, vals, coords)

    def _cmorize_var(self, var_info, vals, coords):
        file_expression = os.path.join(
            self.in_dir, '{0}.H????'.format(vals['raw']))
        for file_path in glob.glob(file_expression):
            cube = PIOMAS._create_cube(
                PIOMAS._read_binary_file(file_path),
                coords,
                int(file_path[-4:]),
                var_info,
                vals['units']
            )
            set_global_atts(cube, self.cfg['attributes'])
            save_variable(
                cube, var_info.short_name, self.out_dir, self.cfg['attributes']
            )

    @staticmethod
    def _create_lat_lon_coords(lat, lon):
        lon_coord = AuxCoord(
            lon,
            standard_name='longitude',
            var_name='lon',
            units='degrees_east'
        )

        lat_coord = AuxCoord(
            lat,
            standard_name='latitude',
            var_name='lat',
            units='degrees_north'
        )
        return Coords(lat_coord, lon_coord)

    @staticmethod
    def _create_cube(data, coords, year, var_info, raw_units):
        time_coord = DimCoord(
            np.arange(0, data.shape[0]),
            standard_name='time',
            var_name='time',
            units=Unit('days since {}-01-01'.format(year), calendar='noleap'),
        )

        cube = iris.cube.Cube(
            data,
            standard_name=var_info.standard_name,
            var_name=var_info.short_name,
            units=raw_units,
        )
        cube.add_dim_coord(time_coord, 0)
        cube.add_aux_coord(coords.lon, (1, 2))
        cube.add_aux_coord(coords.lat, (1, 2))
        return cube

    def _create_areacello(self, coords, var_info):
        cube = iris.cube.Cube(
            self.areacello,
            standard_name=var_info.standard_name,
            var_name=var_info.short_name,
            units='m2',
        )
        cube.add_aux_coord(coords.lon, (0, 1))
        cube.add_aux_coord(coords.lat, (0, 1))
        return cube

    @staticmethod
    def _read_binary_file(data_path, vector=False):
        fd_data = open(data_path, 'rb')
        data = np.fromfile(fd_data, dtype=np.dtype('f'), count=-1)
        days = data.shape[0] // NX // NY
        data = data.reshape(days, NY, NX)
        if vector:
            return data[0:days:2, ...], data[1:days:2, ...]
        return data
