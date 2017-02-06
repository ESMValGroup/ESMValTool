import iris


class CMORCheck(object):

    def __init__(self, cube):
        self.cube = cube
        self.field_type = None
        self._errors = list()

    def check(self):
        self._check_rank()
        self._is_correct()


    def _check_rank(self):
        # Field_type is like T3m or T3Om
        rank = int(self.field_type[1])
        dim_coords = self.cube.coords(dim_coords=True)
        if len(dim_coords) != rank:
            self.report_error('Coordinates does not match')

    def report_error(self, message, *args):
        self._errors.append(message.format(*args))

    def _is_correct(self):
        if len(self._errors) > 0:
            for error in self._errors:
                print(error)
            raise CMORCheckError('There were errors in varible {0}'.format(self.cube.standard_name))


class CMORCheckError(Exception):
    pass

if __name__ == '__main__':
    cube = iris.load_cube('/Users/nube/esmval_data/ETHZ_CMIP5/historical/Amon/ta/CMCC-CESM/r1i1p1/'
                          'ta_Amon_CMCC-CESM_historical_r1i1p1_200001-200212.nc')
    checker = CMORCheck(cube)
    checker.field_type = 'T3M'
    checker.check()
