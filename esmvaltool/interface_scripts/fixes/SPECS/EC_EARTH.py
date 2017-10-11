from esmvaltool.interface_scripts.fixes.fix import Fix
from netCDF4 import Dataset


class ta(Fix):

    def fix_file(self, filepath):
        handler = Dataset(filepath, 'r+')
        handler.variables['ta'].coordinates = 'air_pressure forecast_period'
        handler.close()
        return filepath


class sic(Fix):

    def fix_file(self, filepath):
        handler = Dataset(filepath, 'a')
        handler.contact = 'Pierre-Antoine Bretonniere, pierre-antoine.bretonniere@bsc.es, ' \
                          'Javier Vegas Regidor, javier.vegas@bsc.es'
        handler['sic'].interval_operation = 3600.0
        handler.variables['sic'].coordinates = 'lat lon'
        handler.close()
        return filepath

    def fix_metadata(self, cube):
        cube.units = '1.0'
        print(cube)
        return cube

