import numpy as np
import iris
import os
import json

iris.FUTURE.cell_datetime_objects = True
iris.FUTURE.netcdf_promote = True

FIELD_TYPES = {
    'T3M': ('time', 'air_pressure', 'latitude', 'longitude'),
    'tas': ('longitude', 'latitude', 'time', 'height2m'),
}

DIM_COORD_UNITS = {
    # 'time': 'days since 1950-01-01 00:00:00',
    'time': 'days since 2000-1-1',
    # 'air_pressure': 'hPa',
    'air_pressure': 'Pa',
    # 'longitude': 'degrees_east',
    'longitude': 'degrees',
    # 'latitude': 'degrees_north',
    'latitude': 'degrees',
}


class CMORCheck(object):

    def __init__(self, cube):
        self.cube = cube
        self.field_type = cube.var_name
        self._errors = list()
        self._cmor_tables_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cmip6-cmor-tables',
                                                'Tables')
        self._load_coord_information()

    def _load_coord_information(self):
        json_data = open(os.path.join(self._cmor_tables_folder, 'CMIP6_coordinate.json')).read()
        self.json_data = json.loads(json_data)

    def check(self):
        self._check_rank()

        if len(self._errors) > 0:
            for error in self._errors:
                print(error)
            raise CMORCheckError('There were errors in variable {0}'.format(self.cube.standard_name))

    def _check_rank(self):
        # Field_type is like T3m or T3Om
        dim_names = FIELD_TYPES[self.field_type]

        # Check number of dim_coords matches rank required
        rank = len(dim_names)
        dim_coords = self.cube.coords(dim_coords=True)
        if len(dim_coords) != rank:
            self.report_error('Coordinate rank does not match')

        # Check names of dimensions
        for coord_name in dim_names:
            if not self.cube.coords(coord_name):
                self.report_error('Coordinate {} does not exist', coord_name)

        # Check metadata for dim_coords
        for coord in dim_coords:
            # Check units
            if str(coord.units) != DIM_COORD_UNITS[coord.name()]:
                self.report_error('Units for {0} are {1}, not {2}', coord.name(), coord.units,
                                  DIM_COORD_UNITS[coord.name()])
            # Check calendar for time coordinates
            if coord.name() == 'time':
                if not coord.units.calendar:
                    self.report_error('Coord time does not contain a calendar')
            # Check interval
                # Examining cell_methods
            # Check monotonicity
            if not coord.is_monotonic():
                self.report_error('Coord {} is not monotonic', coord.name())
            if coord.name() == 'air_pressure':
                # Given we know coord is monotonic, can just compare first two points
                if coord.points[0] < coord.points[1]:
                    self.report_error('Coord air_pressure is not decreasing')
            # Check ranges
            if coord.name() == 'longitude':
                if np.any(np.logical_or(coord.points < 0.0, coord.points >= 360.0)):
                    self.report_error('Coord longitude is out of bounds')
            if coord.name() == 'latitude':
                if np.any(np.logical_or(coord.points < -90.0, coord.points > 90.0)):
                    self.report_error('Coord latitude is out of bounds')

    def report_error(self, message, *args):
        self._errors.append(message.format(*args))


class CMORCheckError(Exception):
    pass


def main():
    data_folder = '/Users/nube/esmval_data'
    cube = iris.load_cube(os.path.join(data_folder, 'ETHZ_CMIP5/historical/Amon/ta/CMCC-CESM/r1i1p1',
                                       'ta_Amon_CMCC-CESM_historical_r1i1p1_200001-200212.nc'))
    checker = CMORCheck(cube)
    checker.field_type = 'T3M'
    checker.check()

    checker.check()

if __name__ == '__main__':
    main()
