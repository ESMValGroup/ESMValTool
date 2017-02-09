import numpy as np
import iris
import iris.exceptions
import os
import json
import warnings

iris.FUTURE.cell_datetime_objects = True
iris.FUTURE.netcdf_promote = True


class CMORTable(object):

    # Dictionary to map CMIP5 variable names to CMIP6
    CMIP5_to_CMIP6 = {
        'sic': 'siconc',
        'tro3': 'o3',
    }

    def __init__(self, table, var_name):
        table = self._translate_table_name(table)
        cwd = os.path.dirname(os.path.realpath(__file__))
        self._cmor_tables_folder = os.path.join(cwd, 'cmip6-cmor-tables', 'Tables')
        self._cmor_tables_file = 'CMIP6_{}.json'.format(table)
        self._load_variable_information(var_name)

    def _translate_table_name(self, table):
        if table == 'OImon':
            table = 'SImon'
        return table

    def _load_coord_information(self):
        table_file = os.path.join(self._cmor_tables_folder,
                                  'CMIP6_coordinate.json')
        with open(table_file) as inf:
            json_data = inf.read()
        self._coord = json.loads(json_data)

        # Fill up coordinate axes with CMOR metadata
        self.coords = {}
        for var_name in self.var['dimensions'].split():
            if var_name in self._generic_levels:
                coord = 'generic_level'
                axis = 'Z'
            else:
                coord = self._coord['axis_entry'][var_name]
                axis = coord['axis']
                if not axis:
                    axis = 'none'

            if axis not in self.coords:
                self.coords[axis] = coord
                # Don't look here!
                if axis == 'T':
                    self.coords[axis]['units'] = u'days since 1850-1-1 00:00:00'
            else:
                print('axis {} already exists'.format(axis))

    def _load_variable_information(self, var_name):
        table_file = os.path.join(self._cmor_tables_folder,
                                  self._cmor_tables_file)
        with open(table_file) as inf:
            json_data = inf.read()
        self._var = json.loads(json_data)
        self._generic_levels = self._var['Header']['generic_levels'].split()
        if var_name in self.CMIP5_to_CMIP6:
            new_var_name = self.CMIP5_to_CMIP6[var_name]
            self.var = self._var['variable_entry'][new_var_name]
        else:
            self.var = self._var['variable_entry'][var_name]
        self._load_coord_information()

    def get_frequency(self):
        return self._var['Header']['frequency']


class CMORCheck(object):

    def __init__(self, cube, table, frequency=None, fail_on_error=False):
        self.cube = cube
        self._failerr = fail_on_error
        self._errors = list()
        self._cmor = CMORTable(table, self.cube.var_name)
        if frequency is None:
            frequency = self._cmor.get_frequency()
        self.frequency = frequency

    def check(self):
        self._check_rank()
        self._check_dim_names()
        self._check_coords()
        self._check_time_coord()
        self._check_var_metadata()
        self._check_fill_value()
        self._check_data_range()

        if self.has_errors():
            msg = 'There were errors in variable {0}:\n {1}'
            msg = msg.format(self.cube.standard_name, '\n '.join(self._errors))
            raise CMORCheckError(msg)

    def _check_fill_value(self):
        # Iris removes _FillValue/missing_value information if data has none
        #  of these values. If there are values == _FillValue then it will
        #  be encoded in the numpy.ma object created.
        #
        #  => Very difficult to check!
        pass

    def _check_var_metadata(self):
        # Set generic error message
        msg = '{} for {} is {}, not {}'
        # Check standard_name
        if self._cmor.var['standard_name']:
            if self.cube.standard_name != self._cmor.var['standard_name']:
                self.report_error(msg, 'standard_name',
                                  self.cube.var_name,
                                  self.cube.standard_name,
                                  self._cmor.var['standard_name'])

        # Check units
        if self._cmor.var['units']:
            if str(self.cube.units) != self._cmor.var['units']:
                self.report_error(msg, 'units',
                                  self.cube.var_name,
                                  self.cube.units,
                                  self._cmor.var['units'])

        # Check other variable attributes that match entries in cube.attributes
        attrs = ['positive']
        for attr in attrs:
            if self._cmor.var[attr]:
                if self.cube.attributes[attr] != self._cmor.var[attr]:
                    self.report_error(msg, attr,
                                      self.cube.var_name,
                                      self.cube.attributes[attr],
                                      self._cmor.var[attr])

    def _check_data_range(self):
        # Set generic error message
        msg = 'Variable {} has values {} ({})'
        # Check data is not less than valid_min
        if self._cmor.var['valid_min']:
            valid_min = float(self._cmor.var['valid_min'])
            if np.any(self.cube.data < valid_min):
                self.report_error(msg, '< valid_min',
                                  self.cube.name(), valid_min)
        # Check data is not greater than valid_max
        if self._cmor.var['valid_max']:
            valid_max = float(self._cmor.var['valid_max'])
            if np.any(self.cube.data > valid_max):
                self.report_error(msg, '> valid_max',
                                  self.cube.name(), valid_max)

    def _check_rank(self):
        # Count rank, excluding scalar dimensions
        rank = 0
        for (axis, cmor) in self._cmor.coords.items():
            if cmor == 'generic_level' or not cmor['value']:
                rank += 1
        # Extract dimension coordinates from cube
        dim_coords = self.cube.coords(dim_coords=True)
        # Check number of dimension coords matches rank
        if len(dim_coords) != rank:
            self.report_error('Coordinate rank does not match')

    def _check_dim_names(self):
        for (axis, cmor) in self._cmor.coords.items():
            if axis == 'none':
                axis = None
            if cmor == 'generic_level':
                if self.cube.coords(axis=axis, dim_coords=True):
                    coord = self.cube.coord(axis=axis, dim_coords=True)
                    if not coord.standard_name:
                        self.report_error('generic_level coordinate does not have a standard name')
                else:
                    self.report_error('generic_level coordinate does not exist')
            else:
                var_name = cmor['out_name']
                # Check for var_name in cube coordinates
                #  Don't insist on dim_coords as could be scalar
                if self.cube.coords(var_name=var_name): #, axis=axis):
                    coord = self.cube.coord(var_name=var_name) #, axis=axis)
                    if coord.standard_name != cmor['standard_name']:
                        self.report_error('standard_name for {} is {}, not {}',
                                          var_name,
                                          coord.standard_name,
                                          cmor['standard_name'])
                else:
                    self.report_error('Coordinate {} does not exist', var_name)

    def _check_coords(self):
        for (axis, cmor) in self._cmor.coords.items():
            if axis == 'none':
                axis = None

            # Cannot check generic_level coords as no CMOR information
            #  Do we want to do a basic check for units, etc.?
            if cmor != 'generic_level':
                var_name = cmor['out_name']

                # Get coordinate var_name as it exists!
                try:
                    coord = self.cube.coord(var_name=var_name,
                                            dim_coords=True) #,
                                            # axis=axis)
                except iris.exceptions.CoordinateNotFoundError:
                    continue

                self._check_coord(cmor, coord, var_name)

    def _check_coord(self, cmor, coord, var_name):
        # Check units
        if cmor['units']:
            if str(coord.units) != cmor['units']:
                self.report_error('units for {} is {}, not {}',
                                  var_name, coord.units, cmor['units'])
        self._check_coord_monotonicity_and_direction(cmor, coord, var_name)
        self._check_coord_values(cmor, coord, var_name)

    def _check_coord_monotonicity_and_direction(self, cmor, coord, var_name):
        if not coord.is_monotonic():
            self.report_error('Coord {} is not monotonic', var_name)
        if cmor['stored_direction']:
            if cmor['stored_direction'] == 'increasing':
                if coord.points[0] > coord.points[1]:
                    msg = 'Coord {} is not increasing'
                    self.report_error(msg, var_name)
            elif cmor['stored_direction'] == 'decreasing':
                if coord.points[0] < coord.points[1]:
                    msg = 'Coord {} is not decreasing'
                    self.report_error(msg, var_name)

    def _check_coord_values(self, cmor, coord, var_name):
        # Check requested coordinate values exist in coord.points
        if cmor['requested']:
            cmor_points = [float(val) for val in cmor['requested']]
            coord_points = list(coord.points)
            for point in cmor_points:
                if point not in coord_points:
                    msg = 'Coord {} does not contain {} {}'
                    self.report_error(msg, var_name, str(point),
                                      str(coord.units))

        # Check coordinate value ranges
        if cmor['valid_min']:
            valid_min = float(cmor['valid_min'])
            if np.any(coord.points < valid_min):
                msg = 'Coord {} has values < valid_min ({})'
                self.report_error(msg, var_name, valid_min)
        if cmor['valid_max']:
            valid_max = float(cmor['valid_max'])
            if np.any(coord.points > valid_max):
                msg = 'Coord {} has values > valid_max ({})'
                self.report_error(msg, var_name, valid_max)

    def _check_time_coord(self):
        try:
            coord = self.cube.coord('time', dim_coords=True)  # , axis='T')
        except iris.exceptions.CoordinateNotFoundError:
            return

        if not coord.units.is_time_reference():
            self.report_error('Coord time does not have time reference units')

        if self.frequency:
            tolerance = 0.001
            if self.frequency == 'dec':
                target_interval = (3600 - tolerance, 3660 + tolerance)
            elif self.frequency == 'yr':
                target_interval = (360 - tolerance, 366 + tolerance)
            elif self.frequency == 'mon':
                target_interval = (28 - tolerance, 31 + tolerance)
            elif self.frequency == 'day':
                target_interval = (1 - tolerance, 1 + tolerance)
            elif self.frequency.endswith('hr'):

                frequency = self.frequency[:-2]
                if frequency == 'sub':
                    frequency = 1.0 / 24
                    target_interval = (-tolerance, frequency + tolerance)
                else:
                    frequency = float(frequency) / 24
                    target_interval = (frequency - tolerance, frequency + tolerance)
            else:
                self.report_error('Frequency {} not supported by checker', self.frequency)
                return
            for i in range(len(coord.points)-1):
                interval = coord.points[i + 1] - coord.points[i]
                if interval < target_interval[0] or interval > target_interval[1]:
                    self.report_error('Frequency {} does not mach input data', self.frequency)
                    break

    def has_errors(self):
        return len(self._errors) > 0

    def report_error(self, message, *args):
        msg = message.format(*args)
        if self._failerr:
            raise CMORCheckError(msg)
        else:
            self._errors.append(msg)


class CMORCheckError(Exception):
    pass


def main():
    #data_folder = '/Users/nube/esmval_data'
    data_folder = '/home/paul/ESMValTool/data'
    example_datas = [
        ('ETHZ_CMIP5/historical/Amon/ps/GFDL-ESM2G/r1i1p1', 'ps', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/ps/MIROC5/r1i1p1', 'ps', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/ps/MIROC-ESM/r1i1p1', 'ps', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/ta/CMCC-CESM/r1i1p1', 'ta', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/ta/GFDL-ESM2G/r1i1p1', 'ta', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/ta/bcc-csm1-1/r1i1p1', 'ta', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/tro3/GFDL-ESM2G/r1i1p1', 'tro3', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/tro3/MIROC5/r1i1p1', 'tro3', 'Amon'),
        ('ETHZ_CMIP5/historical/Amon/tro3/MIROC-ESM/r1i1p1', 'tro3', 'Amon'),
        ('ETHZ_CMIP5/historical/OImon/sic/EC-EARTH/r1i1p1', 'sic', 'SImon'),
        ('ETHZ_CMIP5/historical/OImon/sic/HadCM3/r1i1p1', 'sic', 'SImon'),
        ('ETHZ_CMIP5/historical/OImon/sic/MRI-ESM1/r1i1p1', 'sic', 'SImon'),
        ('CMIP6/1pctCO2/Amon/ua/MPI-ESM-LR/r1i1p1f1', 'ua', 'Amon'),
        ('CMIP6/1pctCO2/Amon/tas/MPI-ESM-LR/r1i1p1f1', 'tas', 'Amon'),
        ('CMIP6/1pctCO2/day/tas/MPI-ESM-LR/r1i1p1f1', 'tas', 'day'),
        ('CMIP6/1pctCO2/day/pr/MPI-ESM-LR/r1i1p1f1', 'pr', 'day'),
        ('CMIP6/1pctCO2/cfDay/hur/MPI-ESM-LR/r1i1p1f1', 'hur', 'CFday'),
        ('CMIP6/1pctCO2/LImon/snw/MPI-ESM-LR/r1i1p1f1', 'snw', 'LImon'),
        ('CMIP6/1pctCO2/Lmon/cropFrac/MPI-ESM-LR/r1i1p1f1', 'cropFrac', 'Lmon'),
        ('CMIP6/1pctCO2/Oyr/co3/MPI-ESM-LR/r1i1p1f1', 'co3', 'Oyr'),
        ]

    # Use this callback to fix anything Iris tries to break!
    # noinspection PyUnusedLocal
    def merge_protect_callback(raw_cube, field, filename):
        # Remove attributes that cause issues with merging and concatenation
        for attr in ['creation_date', 'tracking_id', 'history']:
            if attr in raw_cube.attributes:
                del raw_cube.attributes[attr]
        # Iris chooses to change longitude and latitude units to degrees
        #  regardless of value in file
        if raw_cube.coords('longitude'):
            raw_cube.coord('longitude').units = 'degrees_east'
        if raw_cube.coords('latitude'):
            raw_cube.coord('latitude').units = 'degrees_north'

    for (example_data, var_name, table) in example_datas:
        print('\n' + example_data)

        try:
            # Load cubes
            files = os.path.join(data_folder, example_data, '*.nc')
            # Suppress warnings associated with missing 'areacella' measure
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',
                                        'Missing CF-netCDF measure variable',
                                        UserWarning)
                warnings.filterwarnings('ignore',
                                        'Missing CF-netCDF boundary variable',
                                        UserWarning)
                var_lambda = lambda cube: cube.var_name==var_name
                cubes = iris.load(files, iris.Constraint(cube_func=var_lambda),
                                  callback=merge_protect_callback)
            # Concatenate data to single cube
            cube = cubes.concatenate_cube()
            # Create checker for loaded cube
            checker = CMORCheck(cube, table)  # , fail_on_error=True)
            # Run checks
            checker.check()

        except (iris.exceptions.ConstraintMismatchError,
                iris.exceptions.ConcatenateError,
                CMORCheckError) as ex:
            print(ex)

if __name__ == '__main__':
    main()
