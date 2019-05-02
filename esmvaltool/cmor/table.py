"""
CMOR information reader for ESMValTool.

Read variable information from CMOR 2 and CMOR 3 tables and make it easily
available for the other components of ESMValTool
"""
import errno
import glob
import json
import logging
import os

logger = logging.getLogger(__name__)

CMOR_TABLES = {}
"""dict of str, obj: CMOR info objects."""


def read_cmor_tables(cfg_developer):
    """Read cmor tables required in the configuration.

    Parameters
    ----------
    cfg_developer : dict of str
        Parsed config-developer file

    """
    custom = CustomInfo()
    CMOR_TABLES['custom'] = custom

    for table in cfg_developer:
        project = cfg_developer[table]

        cmor_type = project.get('cmor_type', 'CMIP5')
        table_path = project.get('cmor_tables', cmor_type.lower())
        table_path = os.path.expandvars(os.path.expanduser(table_path))

        cmor_strict = project.get('cmor_strict', True)
        if cmor_strict:
            default = None
        else:
            default = custom
        if cmor_type == 'CMIP5':
            CMOR_TABLES[table] = CMIP5Info(
                table_path, default=default,
            )
        elif cmor_type == 'CMIP6':
            CMOR_TABLES[table] = CMIP6Info(
                table_path, default=default,
            )


class CMIP6Info(object):
    """
    Class to read CMIP6-like data request.

    This uses CMOR 3 json format

    Parameters
    ----------
    cmor_tables_path: basestring
        Path to the folder containing the Tables folder with the json files

    """

    _CMIP_5to6_varname = {
        'sic': 'siconc',
        'sit': 'sithick',
        'tro3': 'o3',
    }

    def __init__(self, cmor_tables_path, default=None):
        cmor_tables_path = self._get_cmor_path(cmor_tables_path)

        self._cmor_folder = os.path.join(cmor_tables_path, 'Tables')
        self.default = default

        self.tables = {}

        self._load_coordinates()
        for json_file in glob.glob(os.path.join(self._cmor_folder, '*.json')):
            if 'CV_test' in json_file or 'grids' in json_file:
                continue
            self._load_table(json_file)

    @staticmethod
    def _get_cmor_path(cmor_tables_path):
        if os.path.isdir(cmor_tables_path):
            return cmor_tables_path
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, 'tables', cmor_tables_path)
        return cmor_tables_path

    def _load_table(self, json_file):
        with open(json_file) as inf:
            raw_data = json.loads(inf.read())
            if not self._is_table(raw_data):
                return
            table = TableInfo()
            header = raw_data['Header']
            table.name = header['table_id'][6:].split('_')[-1]
            self.tables[table.name] = table

            generic_levels = header['generic_levels'].split()
            table.frequency = header.get('frequency', '')
            table.realm = header.get('realm', '')

            for var_name, var_data in raw_data['variable_entry'].items():
                var = VariableInfo('CMIP6', var_name)
                if 'frequency' in var_data:
                    var.frequency = var_data['frequency']
                else:
                    var.frequency = table.frequency
                var.read_json(var_data)
                self._assign_dimensions(var, generic_levels)
                table[var_name] = var

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
        for json_file in glob.glob(
                os.path.join(self._cmor_folder, '*coordinate*.json')):
            with open(json_file) as inf:
                table_data = json.loads(inf.read())
                for coord_name in table_data['axis_entry'].keys():
                    coord = CoordinateInfo(coord_name)
                    coord.read_json(table_data['axis_entry'][coord_name])
                    self.coords[coord_name] = coord

    def get_table(self, table):
        """
        Search and return the table info.

        Parameters
        ----------
        table: basestring
            Table name

        Returns
        -------
        TableInfo
            Return the TableInfo object for the requested table if
            found, returns None if not

        """
        return self.tables.get(table)

    def get_variable(self, table, short_name):
        """
        Search and return the variable info.

        Parameters
        ----------
        table: basestring
            Table name
        short_name: basestring
            Variable's short name

        Returns
        -------
        VariableInfo
            Return the VariableInfo object for the requested variable if
            found, returns None if not

        """
        try:
            return self.tables[table][short_name]
        except KeyError:
            if short_name in CMIP6Info._CMIP_5to6_varname:
                new_short_name = CMIP6Info._CMIP_5to6_varname[short_name]
                return self.get_variable(table, new_short_name)
            if self.default:
                return self.default.get_variable(table, short_name)
            return None

    @staticmethod
    def _is_table(table_data):
        if 'variable_entry' not in table_data:
            return False
        if 'Header' not in table_data:
            return False
        return True


class TableInfo(dict):
    """Container class for storing a CMOR table."""

    def __init__(self, *args, **kwargs):
        """Create a new TableInfo object for storing VariableInfo objects."""
        super(TableInfo, self).__init__(*args, **kwargs)
        self.name = ''
        self.frequency = ''
        self.realm = ''


class JsonInfo(object):
    """
    Base class for the info classes.

    Provides common utility methods to read json variables
    """

    def __init__(self):
        self._json_data = {}

    def _read_json_variable(self, parameter):
        """
        Read a json parameter in json_data.

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
        Read a json list parameter in json_data.

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
    """Class to read and store variable information."""

    def __init__(self, table_type, short_name):
        """
        Class to read and store variable information.

        Parameters
        ----------
        short_name: str
            variable's short name

        """
        super(VariableInfo, self).__init__()
        self.table_type = table_type
        self.modeling_realm = []
        """Modeling realm"""
        self.short_name = short_name
        """Short name"""
        self.standard_name = ''
        """Standard name"""
        self.long_name = ''
        """Long name"""
        self.units = ''
        """Data units"""
        self.valid_min = ''
        """Minimum admitted value"""
        self.valid_max = ''
        """Maximum admitted value"""
        self.frequency = ''
        """Data frequency"""
        self.positive = ''
        """Increasing direction"""

        self.dimensions = []
        """List of dimensions"""
        self.coordinates = {}
        """Coordinates"""

        self._json_data = None

    def read_json(self, json_data):
        """
        Read variable information from json.

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
        self.modeling_realm = \
            self._read_json_variable('modeling_realm').split()

        self.dimensions = self._read_json_variable('dimensions').split()


class CoordinateInfo(JsonInfo):
    """Class to read and store coordinate information."""

    def __init__(self, name):
        """
        Class to read and store coordinate information.

        Parameters
        ----------
        name: str
            coordinate's name

        """
        super(CoordinateInfo, self).__init__()
        self.name = name
        self.generic_level = False

        self.axis = ""
        """Axis"""
        self.value = ""
        """Coordinate value"""
        self.standard_name = ""
        """Standard name"""
        self.long_name = ""
        """Long name"""
        self.out_name = ""
        """
        Out name

        This is the name of the variable in the file
        """
        self.var_name = ""
        """Short name"""
        self.units = ""
        """Units"""
        self.stored_direction = ""
        """Direction in which the coordinate increases"""
        self.requested = []
        """Values requested"""
        self.valid_min = ""
        """Minimum allowed value"""
        self.valid_max = ""
        """Maximum allowed value"""

    def read_json(self, json_data):
        """
        Read coordinate information from json.

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


class CMIP5Info(object):
    """
    Class to read CMIP5-like data request.

    Parameters
    ----------
    cmor_tables_path: basestring
       Path to the folder containing the Tables folder with the json files

    """

    def __init__(self, cmor_tables_path, default=None):
        cmor_tables_path = self._get_cmor_path(cmor_tables_path)

        self._cmor_folder = os.path.join(cmor_tables_path, 'Tables')
        if not os.path.isdir(self._cmor_folder):
            raise OSError(errno.ENOTDIR, "CMOR tables path is not a directory",
                          self._cmor_folder)

        self.tables = {}
        self.coords = {}
        self.default = default
        self._current_table = None
        self._last_line_read = None

        for table_file in glob.glob(os.path.join(self._cmor_folder, '*')):
            if '_grids' in table_file:
                continue
            self._load_table(table_file)

    @staticmethod
    def _get_cmor_path(cmor_tables_path):
        if os.path.isdir(cmor_tables_path):
            return cmor_tables_path
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, 'tables', cmor_tables_path)
        return cmor_tables_path

    def _load_table(self, table_file, table_name=''):
        if table_name and table_name in self.tables:
            # special case used for updating a table with custom variable file
            table = self.tables[table_name]
        else:
            # default case: table name is first line of table file
            table = None

        self._read_table_file(table_file, table)

    def _read_table_file(self, table_file, table=None):
        with open(table_file) as self._current_table:
            self._read_line()
            while True:
                key, value = self._last_line_read
                if key == 'table_id':
                    table = TableInfo()
                    table.name = value[len('Table '):]
                    self.tables[table.name] = table
                elif key == 'frequency':
                    table.frequency = value
                elif key == 'modeling_realm':
                    table.realm = value
                elif key == 'generic_levels':
                    for dim in value.split(' '):
                        coord = CoordinateInfo(dim)
                        coord.generic_level = True
                        coord.axis = 'Z'
                        self.coords[dim] = coord
                elif key == 'axis_entry':
                    self.coords[value] = self._read_coordinate(value)
                    continue
                elif key == 'variable_entry':
                    table[value] = self._read_variable(value, table.frequency)
                    continue
                if not self._read_line():
                    return

    def _read_line(self):
        line = self._current_table.readline()
        if line == '':
            return False
        if line.startswith('!'):
            return self._read_line()
        line = line.replace('\n', '')
        if '!' in line:
            line = line[:line.index('!')]
        line = line.strip()
        if not line:
            self._last_line_read = ('', '')
        else:
            index = line.index(':')
            self._last_line_read = (line[:index].strip(),
                                    line[index + 1:].strip())
        return True

    def _read_coordinate(self, value):
        coord = CoordinateInfo(value)
        while self._read_line():
            key, value = self._last_line_read
            if key in ('variable_entry', 'axis_entry'):
                return coord
            if key == 'requested':
                coord.requested = value.split(' ')
                continue
            if hasattr(coord, key):
                setattr(coord, key, value)
        return coord

    def _read_variable(self, short_name, frequency):
        var = VariableInfo('CMIP5', short_name)
        var.frequency = frequency
        while self._read_line():
            key, value = self._last_line_read
            if key in ('variable_entry', 'axis_entry'):
                break
            if key in ('dimensions', 'modeling_realm'):
                setattr(var, key, value.split())
            elif hasattr(var, key):
                setattr(var, key, value)
        for dim in var.dimensions:
            var.coordinates[dim] = self.coords[dim]
        return var

    def get_table(self, table):
        """
        Search and return the table info.

        Parameters
        ----------
        table: basestring
            Table name

        Returns
        -------
        TableInfo
            Return the TableInfo object for the requested table if
            found, returns None if not

        """
        return self.tables.get(table)

    def get_variable(self, table, short_name):
        """
        Search and return the variable info.

        Parameters
        ----------
        table: basestring
            Table name
        short_name: basestring
            Variable's short name

        Returns
        -------
        VariableInfo
            Return the VariableInfo object for the requested variable if
            found, returns None if not

        """
        var_info = self.tables.get(table, {}).get(short_name, None)
        if not var_info and self.default:
            return self.default.get_variable(table, short_name)
        return var_info


class CustomInfo(CMIP5Info):
    """
    Class to read custom var info for ESMVal.

    Parameters
    ----------
    cmor_tables_path: basestring or None
        Full path to the table or name for the table if it is present in
        ESMValTool repository

    """

    def __init__(self, cmor_tables_path=None):
        cwd = os.path.dirname(os.path.realpath(__file__))
        self._cmor_folder = os.path.join(cwd, 'tables', 'custom')
        self.tables = {}
        table = TableInfo()
        table.name = 'custom'
        self.tables[table.name] = table
        self._coordinates_file = os.path.join(
            self._cmor_folder,
            'CMOR_coordinates.dat',
        )
        self.coords = {}
        self._read_table_file(self._coordinates_file, self.tables['custom'])
        for dat_file in glob.glob(os.path.join(self._cmor_folder, '*.dat')):
            if dat_file == self._coordinates_file:
                continue
            self._read_table_file(dat_file, self.tables['custom'])

    def get_table(self, table):
        """
        Search and return the table info.

        Parameters
        ----------
        table: basestring
            Table name

        Returns
        -------
        TableInfo
            Return the TableInfo object for the requested table if
            found, returns None if not

        """
        return self.tables.get(table)

    def get_variable(self, table, short_name):
        """
        Search and return the variable info.

        Parameters
        ----------
        table: basestring
            Table name
        short_name: basestring
            Variable's short name

        Returns
        -------
        VariableInfo
            Return the VariableInfo object for the requested variable if
            found, returns None if not

        """
        return self.tables['custom'].get(short_name, None)

    def _read_table_file(self, table_file, table=None):
        with open(table_file) as self._current_table:
            self._read_line()
            while True:
                key, value = self._last_line_read
                if key == 'generic_levels':
                    for dim in value.split(' '):
                        coord = CoordinateInfo(dim)
                        coord.generic_level = True
                        coord.axis = 'Z'
                        self.coords[dim] = coord
                elif key == 'axis_entry':
                    self.coords[value] = self._read_coordinate(value)
                    continue
                elif key == 'variable_entry':
                    table[value] = self._read_variable(value, None)
                    continue
                if not self._read_line():
                    return
