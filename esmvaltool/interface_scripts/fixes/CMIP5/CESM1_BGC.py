from esmvaltool.interface_scripts.fixes.fix import Fix
from netCDF4 import Dataset
from cf_units import Unit


class nbp(Fix):
    def fix_file(self, filepath):
        original_dataset = Dataset(filepath, mode='a')
        original_var = original_dataset.variables['nbp']
        attr = {k: original_var.getncattr(k) for k in original_var.ncattrs()}
        attr['missing_value'] = 1e+33
        attr['_FillValue'] = 1e+33
        original_var.setncatts(attr)
        original_dataset.close()
        return filepath


class co2(Fix):
    def fix_data(self, cube):
        return cube * 28.966 / 44.0


class allvars(Fix):
    def fix_metadata(self, cube):
        time = cube.coord('time')
        time.units = Unit('days since 1850-01-01 00:00:00', time.units.calendar)
        return cube
