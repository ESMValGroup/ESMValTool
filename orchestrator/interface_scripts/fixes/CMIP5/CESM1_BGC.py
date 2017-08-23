from orchestrator.interface_scripts.fixes.fix import Fix

# Can not be done in the cube

# if (name.eq."nbp") then
#     var@_FillValue = 1.e+33
#     var@missing_value = var@_FillValue
#     ret = 0
# end if


class co2(Fix):

    def fix_data(self, cube):
        return cube * 28.966 / 44.0
