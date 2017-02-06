import iris

FIELD_TYPES = {
'T3M': ('time', 'air_pressure', 'latitude', 'longitude')
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
        rank = len(dim_names)
        dim_coords = self.cube.coords(dim_coords=True)
        # Check number of dim_coords matches rank required
        if len(dim_coords) != rank:
            raise CMORCheckError('Coordinate rank does not match')
        # Check names of dimensions
        for coord_name in dim_names:
            if not self.cube.coords(coord_name):
                raise CMORCheckError('Coordinate {} does not exist'.format(coord_name))

class CMORCheckError(Exception):
    pass

if __name__ == '__main__':
    cube = iris.load_cube('/home/paul/ESMValTool/data/BACKEND-DATA/ETHZ_CMIP5/'
                          'historical/Amon/ta/CMCC-CESM/r1i1p1/'
                          'ta_Amon_CMCC-CESM_historical_r1i1p1_200001-200212.nc')
    checker = CMORCheck(cube)
    checker.field_type = 'T3M'
    print(checker.check())
