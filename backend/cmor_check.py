import numpy as np
import iris
import os
import json

iris.FUTURE.cell_datetime_objects = True
iris.FUTURE.netcdf_promote = True


class CMORCheck(object):

    def __init__(self, cube):
        self.cube = cube
        self._errors = list()
        self._cmor_tables_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cmip6-cmor-tables',
                                                'Tables')
        self._load_variable_information()
        self._load_coord_information()

    def _load_coord_information(self):
        json_data = open(os.path.join(self._cmor_tables_folder, 'CMIP6_coordinate.json')).read()
        self.coord_json_data = json.loads(json_data)
        self.coord_json_data = self.coord_json_data['axis_entry']

    def _load_variable_information(self):
        json_data = open(os.path.join(self._cmor_tables_folder, 'CMIP6_Amon.json')).read()
        self.var_json_data = json.loads(json_data)
        self.var_json_data = self.var_json_data['variable_entry'][self.cube.var_name]

    def check(self):
        self._check_rank()
        self._check_dim_names()
        self._check_coords()
        #self._check_time_coord()
        #self._check_var_metadata()
        #self._check_fill_value()
        #self._check_data_range()

        if len(self._errors) > 0:
            for error in self._errors:
                print(error)
            raise CMORCheckError('There were errors in variable {0}: \n{1}'.format(self.cube.standard_name,
                                                                                   '\n'.join(self._errors)))

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
            dim_names = self.var_json_data['dimensions'].split()
            for coord_name in dim_names:
                # Check for variable name in coordinate CMOR table,
                #  get out_name as that is coordinate variable name
                #  in variable file
                if coord_name in self.coord_json_data:
                    out_var_name = self.coord_json_data[coord_name]['out_name']
                # Check for out_var_name in variable coordinates
                if not self.cube.coords(var_name=out_var_name, dim_coords=True):
                    self.report_error('Coordinate {} does not exist', coord_name)

    def _check_time_coord(self):
        coord = self.cube.coord('time')
        if not coord.units.is_time():
            self.report_error('Coord time does not have time units')
        if not coord.units.calendar:
            self.report_error('Coord time units does not contain a calendar')
        if not coord.is_monotonic():
             self.report_error('Coord time is not monotonic')

    def _check_coords(self):
        # Check metadata for dim_coords
        for coord in self.cube.coords(dim_coords=True):
            if coord.name() != 'time':
                # Get CMOR metadata for coordinate
                cmor_table = self.coord_json_data[coord.var_name]

                # Check units
                if cmor_table['units']:
                    if str(coord.units) != cmor_table['units']:
                        self.report_error('Units for {0} are {1}, not {2}',
                                          coord.name(), coord.units,
                                          cmor_table['units'])

                # Check monotonicity and direction
                if not coord.is_monotonic():
                    self.report_error('Coord {} is not monotonic', coord.name())
                if cmor['stored_direction']:
                    if cmor['stored_direction'] == 'increasing':
                        if coord.points[0] > coord.points[1]:
                            self.report_error('Coord {} is not increasing', coord.name())
                    elif cmor['stored_direction'] == 'decreasing':
                        if coord.points[0] < coord.points[1]:
                            self.report_error('Coord {} is not decreasing', coord.name())

                # Check coordinate value ranges
                if cmor_table['valid_min']:
                    valid_min = float(cmor_table['valid_min'])
                    if np.any(coord.points < valid_min):
                        self.report_error('Coord {} has values < valid_min', coord.name())
                if cmor_table['valid_max']:
                    valid_max = float(cmor_table['valid_max'])
                    if np.any(coord.points > valid_max):
                        self.report_error('Coord {} has values > valid_max', coord.name())

    def report_error(self, message, *args):
        self._errors.append(message.format(*args))


class CMORCheckError(Exception):
    pass


def main():
    #data_folder = '/Users/nube/esmval_data'
    data_folder = '/home/paul/ESMValTool/data/BACKEND-DATA'
    cube = iris.load_cube(os.path.join(data_folder, 'ETHZ_CMIP5/historical/Amon/ta/CMCC-CESM/r1i1p1',
                                       'ta_Amon_CMCC-CESM_historical_r1i1p1_200001-200212.nc'))
    checker = CMORCheck(cube)
    checker.check()

if __name__ == '__main__':
    main()
