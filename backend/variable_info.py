"""
Classes to read variable information from CMIP6 tables and make it easily
available for the other components of ESMValTool
"""
import os
import json
import glob


class CMIP6Info(object):

    _CMIP_5to6_varname = {
        'sic': 'siconc',
        'tro3': 'o3',
    }

    def __init__(self, cmor_tables_path=None):
        """
        Class to read CMIP6-like data request

        Parameters
        ----------
        cmor_tables_path: basestring
            Path to the folder containing the Tables folder with the json files
        """

        cmor_tables_path = self._get_cmor_path(cmor_tables_path)

        self._cmor_folder = os.path.join(cmor_tables_path, 'Tables')

        self.tables = {}

        self._load_coordinates()

        for json_file in glob.glob(os.path.join(self._cmor_folder, '*.json')):
            print('json_file: {}'.format(json_file))
            if 'CV_test' in json_file or 'grids' in json_file:
                continue
            self._load_table(json_file)

    @staticmethod
    def _get_cmor_path(cmor_tables_path):
        if not cmor_tables_path:
            cwd = os.path.dirname(os.path.realpath(__file__))
            cmor_tables_path = os.path.join(cwd, 'cmip6-cmor-tables')
        return cmor_tables_path

    def _load_table(self, json_file):
        with open(json_file) as inf:
            raw_data = json.loads(inf.read())
            if not self._is_table(raw_data):
                return
            header = raw_data['Header']
            name = header['table_id'][6:]
            self.tables[name] = {}

            generic_levels = header['generic_levels'].split()
            frequency = header['frequency']

            for var_name, var_data in raw_data['variable_entry'].items():
                var = VariableInfo(var_name)
                var.frequency = frequency
                var.read_json(var_data)
                self._assign_dimensions(var, generic_levels)
                self.tables[name][var_name] = var

    def _assign_dimensions(self, var, generic_levels):
        for dimension in var.dimensions:
            if dimension in generic_levels:
                coord = CoordinateInfo(dimension)
                coord.generic_level = True
                coord.axis = 'Z'
            else:
                coord = self.coords[dimension]

            axis = coord.axis
            if not axis:
                axis = 'none'

            var.coordinates[axis] = coord

    def _load_coordinates(self):
        self.coords = {}
        for json_file in glob.glob(os.path.join(self._cmor_folder,
                                                '*coordinate*.json')):
            with open(json_file) as inf:
                table_data = json.loads(inf.read())
                for coord_name in table_data['axis_entry'].keys():
                    coord = CoordinateInfo(coord_name)
                    coord.read_json(table_data['axis_entry'][coord_name])
                    self.coords[coord_name] = coord

    def get_variable(self, table, short_name):
        """
        Search and return the variable info

        Parameters
        ----------
        table: basestring
            Table name
        short_name: basestring
            Variable's short name

        Returns
        -------
        VariableInfo
            Returns the VariableInfo object for the requested variable if
            found, returns None if not
        """
        try:
            return self.tables[table][short_name]
        except KeyError:
            if short_name in CMIP6Info._CMIP_5to6_varname:
                new_short_name = CMIP6Info._CMIP_5to6_varname[short_name]
                return self.get_variable(table, new_short_name)
            return None

    @staticmethod
    def _is_table(table_data):
        if 'variable_entry' not in table_data:
            return False
        if 'Header' not in table_data:
            return False
        return True


class JsonInfo(object):
    def __init__(self):
        """
        Base class fo the info classes to provide common utility methods to
        read json variables
        """
        self._json_data = None

    def _read_json_variable(self, parameter):
        """
        Reads a json parameter in json_data

        Parameters
        ----------
        parameter: str
            parameter to read

        Returns
        -------
        str
            Option's value or empty string if parameter is not present
        """
        if parameter not in self._json_data:
            return ''
        return str(self._json_data[parameter])

    def _read_json_list_variable(self, parameter):
        """
        Reads a json list parameter in json_data

        Parameters
        ----------
        parameter: str
            parameter to read

        Returns
        -------
        str
            Option's value or empty list if parameter is not present
        """
        if parameter not in self._json_data:
            return []
        return self._json_data[parameter]


class VariableInfo(JsonInfo):

    def __init__(self, short_name):
        """
        Class to read and store variable information

        Parameters
        ----------
        short_name: str
            variable's short name
        """
        super(VariableInfo, self).__init__()
        self.short_name = short_name
        self.standard_name = ''
        self.long_name = ''
        self.units = ''
        self.valid_min = ''
        self.valid_max = ''
        self.frequency = ''
        self.positive = ''

        self.dimensions = []
        self.coordinates = {}

        self.derived = False
        self.required_vars = []

        self._json_data = None

    def read_json(self, json_data):
        """
        Reads variable information from json.

        Non-present options will be set to empty

        Parameters
        ----------
        json_data: dict
            dictionary created by the json reader containing
            variable information
        """
        self._json_data = json_data

        self.standard_name = self._read_json_variable('standard_name')
        self.long_name = self._read_json_variable('long_name')
        self.units = self._read_json_variable('units')
        self.valid_min = self._read_json_variable('valid_min')
        self.valid_max = self._read_json_variable('valid_max')
        self.positive = self._read_json_variable('positive')

        self.dimensions = self._read_json_variable('dimensions').split()


class CoordinateInfo(JsonInfo):

    def __init__(self, name):
        """
        Class to read and store coordinate information

        Parameters
        ----------
        name: str
            coordinate's name
        """
        super(CoordinateInfo, self).__init__()
        self.name = name
        self.generic_level = False

        self.axis = ""
        self.value = ""
        self.standard_name = ""
        self.long_name = ""
        self.out_name = ""
        self.var_name = ""
        self.units = ""
        self.stored_direction = ""
        self.requested = []
        self.valid_min = ""
        self.valid_max = ""

    def read_json(self, json_data):
        """
        Reads coordinate information from json.

        Non-present options will be set to empty

        Parameters
        ----------
        json_data: dict
            dictionary created by the json reader containing
            coordinate information
        """
        self._json_data = json_data

        self.axis = self._read_json_variable('axis')
        self.value = self._read_json_variable('value')
        self.out_name = self._read_json_variable('out_name')
        self.var_name = self._read_json_variable('var_name')
        self.standard_name = self._read_json_variable('standard_name')
        self.long_name = self._read_json_variable('long_name')
        self.units = self._read_json_variable('units')
        self.stored_direction = self._read_json_variable('stored_direction')
        self.valid_min = self._read_json_variable('valid_min')
        self.valid_max = self._read_json_variable('valid_max')
        self.requested = self._read_json_list_variable('requested')
