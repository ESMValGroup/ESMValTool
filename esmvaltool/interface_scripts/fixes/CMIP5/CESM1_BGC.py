"""Fixes for CESM1-BGC model"""
import shutil
from cf_units import Unit
from netCDF4 import Dataset


from esmvaltool.interface_scripts.fixes.fix import Fix


class nbp(Fix):
    """Fixes for nbp variable"""

    def fix_file(self, filepath, var_path):
        """
        Fix file

        Fixes missing_value and _FillValue of nbp

        Parameters
        ----------
        filepath: str
            path to the file to fix
        var_path: str
            path to the preprocessor folder

        Returns
        -------
        str
            path to fixed file

        """
        temp = Fix.get_fixed_filepath(filepath, preproc_dir)
        shutil.copy(filepath, temp)
        original_dataset = Dataset(temp, mode='a')
        original_var = original_dataset.variables['nbp']
        attr = {k: original_var.getncattr(k) for k in original_var.ncattrs()}
        attr['missing_value'] = 1e+33
        attr['_FillValue'] = 1e+33
        original_var.setncatts(attr)
        original_dataset.close()
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
        time.units = Unit('days since 1850-01-01 00:00:00',
                          time.units.calendar)
        return cube
