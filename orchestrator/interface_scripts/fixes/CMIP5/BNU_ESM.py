from orchestrator.interface_scripts.fixes.fix import Fix
from cf_units import Unit


class fgco2(Fix):

    def fix_metadata(self, cube):
        cube.units = Unit('kg m-2 s-1')
        return cube

    def fix_data(self, cube):
        return cube * 12.0 / 44.0


class ch4(Fix):

    def fix_metadata(self, cube):
        cube.units = Unit('1e-9')
        return cube

    def fix_data(self, cube):
        return cube * 29.0/16.0 * 1.e9


class co2(Fix):

    def fix_metadata(self, cube):
        cube.units = Unit('1e-6')
        return cube

    def fix_data(self, cube):
        return cube * 29.0/44.0 * 1.e6


class spco2(Fix):

    def fix_data(self, cube):
        return cube * 1.e6

# No clear way to apply this fix now that we are working with cubes, not files

# class sftlf(Fix):
#
#     def fix_metadata(self):
#         self.cube = self.cube * 1.e6

#   if (name.eq."sftlf") then
#       files = systemfunc("ls " + INFILE)
#       f=addfile(files(0), "r")
#       tmp=f->lat
#       var&lat = tmp
#       delete(tmp)
#       delete(f)
#       ret = 0
#   end if
#

