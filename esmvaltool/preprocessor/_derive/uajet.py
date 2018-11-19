"""Derivation of variable `uajet`."""

import cf_units
import iris
from iris import Constraint
import numpy as np

from ._derived_variable_base import DerivedVariableBase

# Constants (Southern hemisphere at 850 hPa)
LAT = [-80.0, -30.0]
PLEV = 85000


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `uajet`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'ua',
            'field': 'T3{frequency}'
        }]
    }

    def calculate(self, cubes):
        """Compute latitude of maximum meridional wind speed."""
        # Load cube, extract correct region and perform zonal mean
        ua_cube = cubes.extract_strict(Constraint(name='eastward_wind'))
        ua_cube = ua_cube.interpolate(
            [('air_pressure', PLEV)], scheme=iris.analysis.Linear())
        ua_cube = ua_cube.extract(
            iris.Constraint(latitude=lambda cell: LAT[0] <= cell <= LAT[1]))
        ua_cube = ua_cube.collapsed('longitude', iris.analysis.MEAN)

        # Calculate maximum jet position
        uajet_vals = []
        for time_slice in ua_cube.slices(['latitude']):
            ua_data = time_slice.data

            # Get maximum ua and corresponding index
            idx_max_ua = np.argmax(ua_data)
            slc = slice(idx_max_ua - 1, idx_max_ua + 2)

            # Perform 2nd degree polynomial fit to get maximum jet position
            x_vals = ua_data[slc]
            y_vals = time_slice.coord('latitude').points[slc]
            polyfit = np.polyfit(x_vals, y_vals, 2)
            polynom = np.poly1d(polyfit)
            uajet_vals.append(polynom(np.max(ua_data)))

        uajet_cube = iris.cube.Cube(
            uajet_vals,
            units=cf_units.Unit('degrees_north'),
            dim_coords_and_dims=[(ua_cube.coord('time'), 0)],
            attributes={
                'plev': PLEV,
                'lat_range_0': LAT[0],
                'lat_range_1': LAT[1]
            })

        return uajet_cube
