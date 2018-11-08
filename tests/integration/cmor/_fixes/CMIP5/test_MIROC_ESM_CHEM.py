import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.MIROC_ESM_CHEM import tro3


class TestTro3(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='tro3', units='J')
        self.fix = tro3()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1000)
        self.assertEqual(cube.units, Unit('J'))


# if (name .eq. "tro3") then
#     ; files say unit of ozone is "1e-9" ut unit is actually "1e-6"
#     var = var * 1.0e3
#     if (iscoord(var, "time")) then
#         do it = 1, dimsizes(var&time) - 1
#             if (var&time(it).eq.0) then
#                 tt = tointeger(cd_calendar(var&time(it-1), 0))
#                 tt(0, 1) = tt(0, 1) + 1  ; month
#                 if (tt(0, 1).gt.12) then
#                     tt(0, 1) = 1
#                     tt(0, 0) = tt(0, 0) + 1  ; year
#                 end if
#                 var&time(it) = cd_inv_calendar(\
#                     tt(0, 0), tt(0, 1), tt(0, 2), tt(0, 3), \
#                     tt(0, 4), tt(0, 5), var&time@units, 0)
#             end if
#         end do
#     ret = 0
#     end if
# end if
