import iris


class CMORCheck(object):

    def __init__(self, cube):
        self.cube = cube
        self.field_type = None

    def check(self):
        try:
            self._check_rank()
        except CMORCheckError:
            return False
        return True

    def _check_rank(self):
        # Field_type is like T3m or T3Om
        rank = int(self.field_type[1])
        dim_coords = self.cube.coords(dim_coord=True)
        if len(dim_coords) != rank:
            raise CMORCheckError('Coordinates does not match')


class CMORCheckError(Exception):
    pass

if __name__ == '__main__':
    cube = iris.load('/Users/nube/esmval_data/ETHZ_CMIP5/historical/Amon/ta/CMCC-CESM/r1i1p1/'
                     'ta_Amon_CMCC-CESM_historical_r1i1p1_200001-200212.nc')
    checker = CMORCheck(cube)
    checker.field_type = 'T3M'