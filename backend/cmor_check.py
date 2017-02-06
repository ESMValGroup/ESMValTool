import numpy as np
import iris
iris.FUTURE.cell_datetime_objects = True
iris.FUTURE.netcdf_promote = True

FIELD_TYPES = {
'T3M': ('time', 'air_pressure', 'latitude', 'longitude'),
}

DIM_COORD_UNITS = {
#'time': 'days since 1950-01-01 00:00:00',
'time': 'days since 2000-1-1',
#'air_pressure': 'hPa',
'air_pressure': 'Pa',
#'longitude': 'degrees_east',
'longitude': 'degrees',
#'latitude': 'degrees_north',
'latitude': 'degrees',
}


class CMORCheck(object):

    def __init__(self, cube):
        self.cube = cube
        self.field_type = None

    def check(self):
        try:
            self._check_rank()
        except CMORCheckError,e:
            print e
            return False
        return True

    def _check_rank(self):
        # Field_type is like T3m or T3Om
        dim_names = FIELD_TYPES[self.field_type]

        # Check number of dim_coords matches rank required
        rank = len(dim_names)
        dim_coords = self.cube.coords(dim_coords=True)
        if len(dim_coords) != rank:
            raise CMORCheckError('Coordinate rank does not match')

        # Check names of dimensions
        for coord_name in dim_names:
            if not self.cube.coords(coord_name):
                raise CMORCheckError('Coordinate {} does not exist'.format(coord_name))

        # Check metadata for dim_coords
        for coord in dim_coords:
            # Check units
            if str(coord.units) != DIM_COORD_UNITS[coord.name()]:
                raise CMORCheckError('Units for {0} are {1}, not {2}'.format(coord.name(), coord.units, DIM_COORD_UNITS[coord.name()]))
            # Check calendar for time coordinates
            if coord.name() == 'time':
                if not coord.units.calendar:
                    raise CMORCheckError('Coord time does not contain a calendar')
            # Check interval
                # Examining cell_methods
            # Check monotonicity
            if not coord.is_monotonic():
                raise CMORCheckError('Coord {} is not monotonic'.format(coord.name()))
            if coord.name() == 'air_pressure':
                # Given we know coord is monotonic, can just compare first two points
                if coord.points[0] < coord.points[1]:
                    raise CMORCheckError('Coord air_pressure is not decreasing')
            # Check ranges
            if coord.name() == 'longitude':
                if np.any(np.logical_or(coord.points < 0.0, coord.points >= 360.0)):
                    raise CMORCheckError('Coord longitude is out of bounds')
            if coord.name() == 'latitude':
                if np.any(np.logical_or(coord.points < -90.0, coord.points > 90.0)):
                    raise CMORCheckError('Coord latitude is out of bounds')




class CMORCheckError(Exception):
    pass

if __name__ == '__main__':
    cube = iris.load_cube('/home/paul/ESMValTool/data/BACKEND-DATA/ETHZ_CMIP5/'
                          'historical/Amon/ta/CMCC-CESM/r1i1p1/'
                          'ta_Amon_CMCC-CESM_historical_r1i1p1_200001-200212.nc')
    #print(cube)
    #print(cube.coords())
    checker = CMORCheck(cube)
    checker.field_type = 'T3M'
    print(checker.check())
