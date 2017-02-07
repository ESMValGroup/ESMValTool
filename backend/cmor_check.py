import numpy as np
import iris
import os
import json

iris.FUTURE.cell_datetime_objects = True
iris.FUTURE.netcdf_promote = True


def merge_protect_callback(cube, field, file):
    for attr in ['creation_date', 'tracking_id', 'history']:
        if attr in cube.attributes:
            del cube.attributes[attr]


class CMORCheck(object):

    def __init__(self, cube, fail_on_error=False):
        self.cube = cube
        self._failerr = fail_on_error
        self._errors = list()
        self._cmor_tables_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                'cmip6-cmor-tables',
                                                'Tables')
        self._load_variable_information()
        self._load_coord_information()

    def _load_coord_information(self):
        json_data = open(os.path.join(self._cmor_tables_folder,
                                      'CMIP6_coordinate.json')).read()
        self.coord_json_data = json.loads(json_data)
        self.coord_json_data = self.coord_json_data['axis_entry']
        # Don't look here!
        self.coord_json_data['time']['units'] = u'days since 1850-1-1 00:00:00'

    def _load_variable_information(self):
        json_data = open(os.path.join(self._cmor_tables_folder,
                                      'CMIP6_Amon.json')).read()
        self.var_json_data = json.loads(json_data)
        self.var_json_data = self.var_json_data['variable_entry'][self.cube.var_name]

    def check(self):
        self._check_rank()
        self._check_dim_names()
        self._check_coords()
        self._check_time_coord()
        #self._check_var_metadata()
        #self._check_fill_value()
        #self._check_data_range()

        if len(self._errors) > 0:
            raise CMORCheckError('There were errors in variable {0}:\n  {1}'.format(self.cube.standard_name,
                                                                                    '\n  '.join(self._errors)))

    def _check_rank(self):
        if self.var_json_data['dimensions']:
            # Check number of dim_coords matches rank required
            dim_names = self.var_json_data['dimensions'].split()
            rank = len(dim_names)
            dim_coords = self.cube.coords(dim_coords=True)
            if len(dim_coords) != rank:
                self.report_error('Coordinate rank does not match')

    def _check_dim_names(self):
        if self.var_json_data['dimensions']:
            # Get names of dimension variables from variable CMOR table
            var_names = self.var_json_data['dimensions'].split()
            for var_name in var_names:
                # Get CMOR metadata for coordinate
                cmor = self.coord_json_data[var_name]
                # Check for variable name in coordinate CMOR table,
                #  get out_name as that is coordinate variable name
                #  in variable file.
                out_var_name = cmor['out_name']
                # Check for out_var_name in variable coordinates
                if not self.cube.coords(var_name=out_var_name, dim_coords=True):
                    self.report_error('Coordinate {} does not exist', var_name)
                else:
                    coord = self.cube.coord(var_name=out_var_name, dim_coords=True)
                    if coord.standard_name != cmor['standard_name']:
                        self.report_error('standard_name for {0} is {1}, not {2}',
                                          var_name, coord.standard_name,
                                          cmor['standard_name'])

    def _check_time_coord(self):
        coord = self.cube.coord('time', dim_coords=True)
        if not coord.units.is_time():
            self.report_error('Coord time does not have time units')
        if not coord.units.calendar:
            self.report_error('Coord time units does not contain a calendar')

    def _check_coords(self):
        if self.var_json_data['dimensions']:
            # Get names of dimension variables from variable CMOR table
            var_names = self.var_json_data['dimensions'].split()
            for var_name in var_names:
                # Get CMOR metadata for coordinate
                cmor = self.coord_json_data[var_name]
                # Check for variable name in coordinate CMOR table,
                #  get out_name as that is coordinate variable name
                #  in variable file.
                out_var_name = cmor['out_name']
                # Get coordinate out_var_name as it exists!
                coord = self.cube.coord(var_name=out_var_name, dim_coords=True)

                # Check units
                if cmor['units']:
                    if str(coord.units) != cmor['units']:
                        self.report_error('Units for {0} is {1}, not {2}',
                                          var_name, coord.units,
                                          cmor['units'])

                # Check monotonicity and direction
                if not coord.is_monotonic():
                    self.report_error('Coord {} is not monotonic', var_name)
                if cmor['stored_direction']:
                    if cmor['stored_direction'] == 'increasing':
                        if coord.points[0] > coord.points[1]:
                            self.report_error('Coord {} is not increasing', var_name)
                    elif cmor['stored_direction'] == 'decreasing':
                        if coord.points[0] < coord.points[1]:
                            self.report_error('Coord {} is not decreasing', var_name)

                # Check requested coordinate values exist in coord.points
                if cmor['requested']:
                    cmor_points = [float(val) for val in cmor['requested']]
                    coord_points = list(coord.points)
                    for point in cmor_points:
                        if point not in coord_points:
                            self.report_error('Coord {0} does not contain {1} {2}',
                                              var_name, str(point),
                                              str(coord.units))

                # Check coordinate value ranges
                if cmor['valid_min']:
                    valid_min = float(cmor['valid_min'])
                    if np.any(coord.points < valid_min):
                        self.report_error('Coord {} has values < valid_min', var_name)
                if cmor['valid_max']:
                    valid_max = float(cmor['valid_max'])
                    if np.any(coord.points > valid_max):
                        self.report_error('Coord {} has values > valid_max', var_name)

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
    #example_data = 'ETHZ_CMIP5/historical/Amon/ta/CMCC-CESM/r1i1p1'
    example_data = 'CMIP6/1pctCO2/Amon/ua/MPI-ESM-LR/r1i1p1f1'

    try:
        # Load cubes
        files = os.path.join(data_folder, example_data, '*.nc')
        cubes = iris.load(files, callback=merge_protect_callback)
        # Concatenate data to single cube to deal with multiple netcdf files
        cube = cubes.concatenate_cube()
        # Create checker for loaded cube
        checker = CMORCheck(cube) #, fail_on_error=True)
        # Run checks
        checker.check()

    except iris.exceptions.ConstraintMismatchError, ex:
        print ex

    except iris.exceptions.ConcatenateError, ex:
        print ex

    except CMORCheckError, ex:
        print(ex)

if __name__ == '__main__':
    main()
