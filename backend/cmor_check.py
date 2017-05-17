import numpy as np
import iris
import iris.coords
import iris.exceptions
import iris.coord_categorisation
import os
import warnings
import cf_units

from backend.variable_info import VariablesInfo

iris.FUTURE.cell_datetime_objects = True
iris.FUTURE.netcdf_promote = True


class CMORCheck(object):
    """
    Class used to check the CMOR-compliance of the data.
    It can also fix some minor errors
    """

    _attr_msg = '{}: {} should be {}, not {}'
    _does_msg = '{}: does not {}'
    _is_msg = '{}: is not {}'
    _vals_msg = '{}: has values {} {}'
    _contain_msg = '{}: does not contain {} {}'

    def __init__(self, cube, table, variables_info, frequency=None,
                 fail_on_error=False, automatic_fixes=False):
        var_info = variables_info.get_variable(table, cube.var_name)
        if var_info is None:
            raise CMORCheckError('Variable {0} was not '
                                 'recognized'.format(cube.var_name))

        self.cube = cube
        self._failerr = fail_on_error
        self._errors = list()
        self._warnings = list()
        self._cmor_var = var_info
        if frequency is None:
            frequency = self._cmor_var.frequency
        self.frequency = frequency
        self.automatic_fixes = automatic_fixes

    def check_metadata(self):
        """
        Checks the data. Raises CMORCheckError if an error is found.
        """
        self._check_rank()
        self._check_var_metadata()
        self._check_fill_value()
        self._check_dim_names()
        self._check_coords()
        self._check_time_coord()

        if self.has_warnings():
            msg = 'There were warnings in variable {0}:\n {1}'
            msg = msg.format(self.cube.var_name, '\n '.join(self._warnings))
            print(msg)

        if self.has_errors():
            msg = 'There were errors in variable {0}:\n {1}'
            msg = msg.format(self.cube.var_name, '\n '.join(self._errors))
            raise CMORCheckError(msg)
        
        self._add_auxiliar_time_coordinates()

    def check_data(self):
        """
        Checks the data. Raises CMORCheckError if an error is found.
        """
        self._check_data_range()

        # Check units
        if self._cmor_var.units:
            if str(self.cube.units) != self._cmor_var.units:
                self.cube.convert_units(self._cmor_var.units)

        if self.has_errors():
            msg = 'There were errors in variable {0}:\n {1}'
            msg = msg.format(self.cube.var_name, '\n '.join(self._errors))
            raise CMORCheckError(msg)

    def _check_fill_value(self):
        # Iris removes _FillValue/missing_value information if data has none
        #  of these values. If there are values == _FillValue then it will
        #  be encoded in the numpy.ma object created.
        #
        #  => Very difficult to check!
        pass

    def _check_var_metadata(self):

        # Check standard_name
        if self._cmor_var.standard_name:
            if self.cube.standard_name != self._cmor_var.standard_name:
                self.report_error(self._attr_msg, self.cube.var_name,
                                  'standard_name',
                                  self._cmor_var.standard_name,
                                  self.cube.standard_name)

        # Check units
        if self._cmor_var.units:
            if not self.cube.units.is_convertible(self._cmor_var.units):
                self.report_error('Variable {0} units () can not be '
                                  'converted to {2}',
                                  self.cube.var_name, self._cmor_var.units,
                                  self.cube.units)

        # Check other variable attributes that match entries in cube.attributes
        attrs = ('positive',)
        for attr in attrs:
            attr_value = getattr(self._cmor_var, attr)
            if attr_value:
                if self.cube.attributes[attr] != attr_value:
                    self.report_error(self._attr_msg, self.cube.var_name, attr,
                                      attr_value, self.cube.attributes[attr])

    def _check_data_range(self):
        # Check data is not less than valid_min
        if self._cmor_var.valid_max:
            valid_min = float(self._cmor_var.valid_max)
            if np.any(self.cube.data < valid_min):
                self.report_error(self._vals_msg, self.cube.var_name,
                                  '< {} ='.format('valid_min'), valid_min)
        # Check data is not greater than valid_max
        if self._cmor_var.valid_max:
            valid_max = float(self._cmor_var.valid_max)
            if np.any(self.cube.data > valid_max):
                self.report_error(self._vals_msg, self.cube.var_name,
                                  '> {} ='.format('valid_max'), valid_max)

    def _check_rank(self):
        # Count rank, excluding scalar dimensions
        rank = 0
        for (axis, coordinate) in self._cmor_var.coordinates.items():
            if coordinate.generic_level or not coordinate.value:
                rank += 1
        # Extract dimension coordinates from cube
        dim_coords = self.cube.coords(dim_coords=True)
        # Check number of dimension coords matches rank
        if len(dim_coords) != rank:
            self.report_error(self._does_msg, self.cube.var_name,
                              'match coordinate rank')

    def _check_dim_names(self):
        for (axis, coordinate) in self._cmor_var.coordinates.items():
            if axis == 'none':
                axis = None
            if coordinate.generic_level:
                var_name = 'generic_level'
                if self.cube.coords(axis=axis, dim_coords=True):
                    cube_coord = self.cube.coord(axis=axis, dim_coords=True)
                    if not cube_coord.standard_name:
                        self.report_error(self._does_msg, var_name,
                                          'have standard_name')
                    # Test for units that match standard_name?
                    # Check for attributes that must exist for a generic_level
                    attrs = []  # 'positive']
                    for attr in attrs:
                        if attr not in self.cube.attributes:
                            self.report_error(self._does_msg, var_name,
                                              'have {}'.format(attr))
                else:
                    self.report_error(self._does_msg, var_name, 'exist')
            else:
                try:
                    cube_coord = self.cube.coord(var_name=coordinate.out_name)
                    if cube_coord.standard_name != coordinate.standard_name:
                        self.report_error(self._attr_msg, coordinate.out_name,
                                          'standard_name',
                                          coordinate.standard_name,
                                          cube_coord.standard_name)
                except iris.exceptions.CoordinateNotFoundError:
                    self.report_error(self._does_msg, coordinate.name, 'exist')

    def _check_coords(self):
        for (axis, coordinate) in self._cmor_var.coordinates.items():
            # Cannot check generic_level coords as no CMOR information
            #  Do we want to do a basic check for units, etc.?
            if coordinate.generic_level:
                continue
            var_name = coordinate.out_name

            # Get coordinate var_name as it exists!
            try:
                coord = self.cube.coord(var_name=var_name,
                                        dim_coords=True)  # ,
                #                        axis=axis)
            except iris.exceptions.CoordinateNotFoundError:
                continue

            self._check_coord(coordinate, coord, var_name)

    def _check_coord(self, cmor, coord, var_name):
        if coord.var_name == 'time':
            return
        # Check units
        if cmor.units:
            if str(coord.units) != cmor.units:
                fixed = False
                if self.automatic_fixes:
                    try:
                        new_unit = cf_units.Unit(cmor.units,
                                                 coord.units.calendar)
                        coord.convert_units(new_unit)
                        fixed = True
                    except ValueError:
                        pass
                if not fixed:
                    self.report_error(self._attr_msg, var_name, 'units',
                                      cmor.units, coord.units)
        self._check_coord_monotonicity_and_direction(cmor, coord, var_name)
        self._check_coord_values(cmor, coord, var_name)

    def _check_coord_monotonicity_and_direction(self, cmor, coord, var_name):

        if not coord.is_monotonic():
            self.report_error(self._is_msg, var_name, 'monotonic')

        if cmor.stored_direction:
            if cmor.stored_direction == 'increasing':
                if coord.points[0] > coord.points[1]:
                    self.report_error(self._is_msg, var_name, 'increasing')

            elif cmor.stored_direction == 'decreasing':
                if coord.points[0] < coord.points[1]:
                    self.report_error(self._is_msg, var_name, 'decreasing')

    def _check_coord_values(self, coord_info, coord, var_name):
        # Check requested coordinate values exist in coord.points
        if coord_info.requested:
            cmor_points = [float(val) for val in coord_info.requested]
            coord_points = list(coord.points)
            for point in cmor_points:
                if point not in coord_points:
                    self.report_warning(self._contain_msg, var_name,
                                        str(point), str(coord.units))

        l_fix_coord_value = False

        # Check coordinate value ranges
        if coord_info.valid_min:
            valid_min = float(coord_info.valid_min)
            if np.any(coord.points < valid_min):
                if coord_info.standard_name == 'longitude' and \
                        self.automatic_fixes:
                    l_fix_coord_value = True
                else:
                    self.report_error(self._vals_msg, var_name,
                                      '< {} ='.format('valid_min'), valid_min)

        if coord_info.valid_max:
            valid_max = float(coord_info.valid_max)
            if np.any(coord.points > valid_max):
                if coord_info.standard_name == 'longitude' and \
                        self.automatic_fixes:
                    l_fix_coord_value = True
                else:
                    self.report_error(self._vals_msg, var_name,
                                      '> {} ='.format('valid_max'), valid_max)

        if l_fix_coord_value:
            lon_extent = iris.coords.CoordExtent(coord, 0.0, 360., True, False)
            self.cube = self.cube.intersection(lon_extent)

    def _check_time_coord(self):
        try:
            coord = self.cube.coord('time', dim_coords=True)  # , axis='T')
            var_name = coord.var_name
        except iris.exceptions.CoordinateNotFoundError:
            return

        if not coord.units.is_time_reference():
            self.report_error(self._does_msg, var_name,
                              'have time reference units')
        else:
            coord.convert_units(cf_units.Unit('days since 1950-01-01 00:00:00',
                                              calendar=coord.units.calendar))
            simplified_cal = self._simplify_calendars(coord.units.calendar)
            coord.units = cf_units.Unit(coord.units.name,
                                        simplified_cal)

        if self.frequency:
            tol = 0.001
            intervals = {'dec': (3600, 3660),
                         'yr': (360, 366),
                         'mon': (28, 31),
                         'day': (1,1)}
            if self.frequency in intervals:
                interval = intervals[self.frequency]
                target_interval = (interval[0] - tol, interval[1] + tol)
            elif self.frequency.endswith('hr'):

                frequency = self.frequency[:-2]
                if frequency == 'sub':
                    frequency = 1.0 / 24
                    target_interval = (-tol, frequency + tol)
                else:
                    frequency = float(frequency) / 24
                    target_interval = (frequency - tol, frequency + tol)
            else:
                msg = '{}: Frequency {} not supported by checker'
                self.report_error(msg, var_name, self.frequency)
                return
            for i in range(len(coord.points)-1):
                interval = coord.points[i + 1] - coord.points[i]
                if (interval < target_interval[0] or
                   interval > target_interval[1]):
                    msg = '{}: Frequency {} does not match input data'
                    self.report_error(msg, var_name, self.frequency)
                    break

    CALENDARS = [['standard', 'gregorian'],
                 ['proleptic_gregorian'],
                 ['noleap', '365_day'],
                 ['all_leap', '366_day'],
                 ['360_day'],
                 ['julian'],
                 ['none']]

    @staticmethod
    def _simplify_calendars(calendar):
        for calendar_type in CMORCheck.CALENDARS:
            if calendar in calendar_type:
                return calendar_type[0]

    def has_errors(self):
        return len(self._errors) > 0

    def has_warnings(self):
        return len(self._warnings) > 0

    def report_error(self, message, *args):
        msg = message.format(*args)
        if self._failerr:
            raise CMORCheckError(msg)
        else:
            self._errors.append(msg)

    def report_warning(self, message, *args):
        msg = message.format(*args)
        if self._failerr:
            print('WARNING: {0}'.format(msg))
        else:
            self._warnings.append(msg)

    def _add_auxiliar_time_coordinates(self):
        coords = [coord.name() for coord in self.cube.aux_coords]
        if 'day_of_month' not in coords:
            iris.coord_categorisation.add_day_of_month(self.cube, 'time')
        if 'day_of_year' not in coords:
            iris.coord_categorisation.add_day_of_year(self.cube, 'time')
        if 'month_number' not in coords:
            iris.coord_categorisation.add_month_number(self.cube, 'time')
        if 'year' not in coords:
            iris.coord_categorisation.add_year(self.cube, 'time')


class CMORCheckError(Exception):
    pass


def main():
    data_folder = '/Users/nube/esmval_data'
    # data_folder = '/home/paul/ESMValTool/data'
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
        # ('CMIP6/1pctCO2/Amon/ua/MPI-ESM-LR/r1i1p1f1', 'ua', 'Amon'),
        # ('CMIP6/1pctCO2/Amon/tas/MPI-ESM-LR/r1i1p1f1', 'tas', 'Amon'),
        # ('CMIP6/1pctCO2/day/tas/MPI-ESM-LR/r1i1p1f1', 'tas', 'day'),
        # ('CMIP6/1pctCO2/day/pr/MPI-ESM-LR/r1i1p1f1', 'pr', 'day'),
        # ('CMIP6/1pctCO2/cfDay/hur/MPI-ESM-LR/r1i1p1f1', 'hur', 'CFday'),
        # ('CMIP6/1pctCO2/LImon/snw/MPI-ESM-LR/r1i1p1f1', 'snw', 'LImon'),
        # ('CMIP6/1pctCO2/Lmon/cropFrac/MPI-ESM-LR/r1i1p1f1', 'cropFrac',
        #  'Lmon'),
        # ('CMIP6/1pctCO2/Oyr/co3/MPI-ESM-LR/r1i1p1f1', 'co3', 'Oyr'),
        ]

    def get_attr_from_field_coord(ncfield, coord_name, attr):
        if coord_name is not None:
            attrs = ncfield.cf_group[coord_name].cf_attrs()
            attr_val = [value for (key, value) in attrs if key == attr]
            if attr_val:
                return attr_val[0]
        return None

    # Use this callback to fix anything Iris tries to break!
    # noinspection PyUnusedLocal
    def merge_callback(raw_cube, field, filename):
        # Remove attributes that cause issues with merging and concatenation
        for attr in ['creation_date', 'tracking_id', 'history']:
            if attr in raw_cube.attributes:
                del raw_cube.attributes[attr]

        # Updates to cube
        # Iris removes several attributes, attempting to re-add them
        # attrs = ['_FillValue', 'missing_value']
        # for attr in attrs:
        #     if attr not in raw_cube.attributes:
        #         attrval = get_attr_from_field(field, attr)
        #         if attrval is not None:
        #             raw_cube.attributes[attr] = attrval

        # Updates to coordinates
        for coord in raw_cube.coords():
            # Iris chooses to change longitude and latitude units to degrees
            #  regardless of value in file, so reinstating file value
            if coord.standard_name in ['longitude', 'latitude']:
                units = get_attr_from_field_coord(field,
                                                  coord.var_name,
                                                  'units')
                if units is not None:
                    coord.units = units
    variables_info = VariablesInfo()
    for (example_data, var_name, table) in example_datas:
        print('\n' + example_data)

        try:
            # Load cubes for requested variable in given files
            files = os.path.join(data_folder, example_data, '*.nc')
            # Suppress warnings associated with missing cell measure and
            #  boundary variables. This may be a problem long term as these
            #  variables are expected to be in the same NetCDF files as the
            #  model data. This is not the case with CMIP data.
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',
                                        'Missing CF-netCDF measure variable',
                                        UserWarning)
                warnings.filterwarnings('ignore',
                                        'Missing CF-netCDF boundary variable',
                                        UserWarning)

                def cube_var_name(raw_cube):
                    return raw_cube.var_name == var_name
                var_cons = iris.Constraint(cube_func=cube_var_name)
                cubes = iris.load(files, var_cons, callback=merge_callback)
            # Concatenate data to single cube, i.e. merge time series
            cube = cubes.concatenate_cube()
            # Create checker for loaded cube
            checker = CMORCheck(cube, table, variables_info,
                                automatic_fixes=True)  # ,
            #                     fail_on_error=True)
            # Run checks
            checker.check_metadata()
            checker.check_data()

        except (iris.exceptions.ConstraintMismatchError,
                iris.exceptions.ConcatenateError,
                CMORCheckError) as ex:
            print(ex)

if __name__ == '__main__':
    main()
