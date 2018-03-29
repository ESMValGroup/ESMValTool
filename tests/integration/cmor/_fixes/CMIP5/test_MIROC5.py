import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.MIROC5 import sftof


class TestGpp(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='sftof', units='J')
        self.fix = sftof()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))

    # dayspermonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    #
    # if ((name.eq."snc".or.name.eq."snw").and.FIELD.eq."T2Ds".and. \
    #     ENSEMBLE.eq."r1i1p1") then
    #     opt = 0
    #     opt@calendar = var&time@calendar
    #     ; get start date from time attribute "days_since_xxxx"
    #     t = 0.0
    #     t@calendar = var&time@calendar
    #     t@units = var&time@units
    #     res = cd_calendar(t, -5)
    #     yy = res(0, 0)
    #     mm = res(0, 1)
    #     dd = res(0, 2)
    #     do ii = 0, dimsizes(var&time) - 1
    #         var&time(ii) = tofloat(cd_inv_calendar(yy, mm, dd, 12, 0, 0, \
    #                                var&time@units, opt))
    #         dd = dd + 1
    #         if (dd.gt.dayspermonth(mm-1)) then
    #             mm = mm + 1
    #             dd = 1
    #         end if
    #         if (mm.gt.12) then
    #             mm = 1
    #             yy = yy + 1
    #         end if
    #     end do
    #     ret = 0
    # end if
