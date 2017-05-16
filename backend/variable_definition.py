import os
import json
import glob


class VariablesInfo(object):
    _CMIP_5to6_varname = {
        'sic': 'siconc',
        'tro3': 'o3',
    }

    def __init__(self, cmor_tables_path=None):

        if not cmor_tables_path:
            cwd = os.path.dirname(os.path.realpath(__file__))
            cmor_tables_path = os.path.join(cwd, 'cmip6-cmor-tables')
        self._cmor_folder = os.path.join(cmor_tables_path, 'Tables')

        self.table = {}

        for json_file in glob.glob(os.path.join(self._cmor_folder, '*coordinate*.json')):
            with open(json_file) as inf:
                table_data = json.loads(inf.read())
                if 'coordinate' in json_file:
                    self.coords = table_data
                    continue

        for json_file in glob.glob(os.path.join(self._cmor_folder, '*.json')):
            print('json_file: {}'.format(json_file))
            if 'CV_test' in json_file or 'grids' in json_file:
                continue
            with open(json_file) as inf:
                table_data = json.loads(inf.read())
                if 'variable_entry' not in table_data or 'Header' not in table_data or\
                                'table_id' not in table_data['Header']:
                    continue
                table_name = table_data['Header']['table_id'][6:]
                if 'generic_levels' in table_data['Header']:
                    generic_levels = table_data['Header']['generic_levels'].split()
                else:
                    generic_levels = ()
                self.table[table_name] = {}
                if 'frequency' in table_data['Header']:
                    frequency = table_data['Header']['frequency']
                else:
                    frequency = ''

                for var_name, var_data in table_data['variable_entry'].items():
                    var = Variable(var_name)
                    var.read_json(var_data)
                    var.frequency = frequency
                    for dimension in var.dimensions:
                        if dimension in generic_levels:
                            coord = 'generic_level'
                            axis = 'Z'
                        else:
                            coord = self.coords['axis_entry'][dimension]
                            axis = coord['axis']
                            if not axis:
                                axis = 'none'

                        if axis not in self.coords:
                            var.coordinates[axis] = coord

                    self.table[table_name][var_name] = var

    def get_variable(self, table, short_name):
        try:
            return self.table[table][short_name]
        except KeyError:
            if short_name in VariablesInfo._CMIP_5to6_varname:
                return self.get_variable(table,
                                         VariablesInfo._CMIP_5to6_varname[short_name])
            return None


class Variable(object):

    def __init__(self, short_name):
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
        self._json_data = json_data

        self.standard_name = self._read_json_variable('standard_name')
        self.long_name = self._read_json_variable('long_name')
        self.units = self._read_json_variable('units')
        self.valid_min = self._read_json_variable('valid_min')
        self.valid_max = self._read_json_variable('valid_max')
        self.positive = self._read_json_variable('positive')

        self.dimensions = self._read_json_variable('dimensions').split()


    def _read_json_variable(self, var_name):
        if var_name not in self._json_data:
            return ''
        return str(self._json_data[var_name])

    def compute(self):
        pass

if __name__ == '__main__':
    var_info = VariablesInfo()
    pass
