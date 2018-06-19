"""Fixes for CESM1-BGC model"""
import shutil
import six
from cf_units import Unit
from netCDF4 import Dataset

from ..fix import Fix


class nbp(Fix):
    """Fixes for nbp variable"""

    def fix_file(self, filepath, output_dir):
        """
        Apply fixes to the files prior to creating the cube.

        Should be used only to fix errors that prevent loading or can
        not be fixed in the cube (i.e. those related with missing_value
        and _FillValue)

        Parameters
        ----------
        filepath: basestring
            file to fix
        output_dir: basestring
            path to the folder to store the fixe files, if required

        Returns
        -------
        basestring
            Path to the corrected file. It can be different from the original
            filepath if a fix has been applied, but if not it should be the
            original filepath

        """
        temp = Fix.get_fixed_filepath(output_dir, filepath)

        original_dataset = Dataset(filepath)
        new_dataset = Dataset(temp, mode='w')

        for dim_name, dimension in six.iteritems(original_dataset.dimensions):
            new_dataset.createDimension(dim_name, dimension.size)

        for var_name, variable in six.iteritems(original_dataset.variables):
            fill_value = variable._FillValue
            if var_name == 'nbp':
                fill_value = 1e+33
            new_var = new_dataset.createVariable(var_name, variable.datatype,
                                                 variable.dimensions,
                                                 zlib=True,
                                                 fill_value=fill_value)
            attr = {k: variable.getncattr(k) for k in variable.ncattrs()}
            del attr['_FillValue']
            attr['missing_value'] = 1e+33
            new_var.setncatts(attr)
            new_var[...] = variable[...]
        original_dataset.close
        new_dataset.close()
        return temp


class co2(Fix):
    """Fixes for co2 variable"""

    def fix_data(self, cube):
        """
        Fix data

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        return cube * 28.966 / 44.0


class allvars(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes time units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        time = cube.coord('time')
        if time.units.name == 'day since 1-01-01 00:00:00.000000 UTC':
            time.units = Unit('days since 1850-01-01 00:00:00',
                              time.units.calendar)
        return cube
