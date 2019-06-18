# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for MIROC5 model."""
import numpy as np
from ..fix import Fix


class sftof(Fix):
    """Fixes for sftof."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 100
        cube.metadata = metadata
        return cube


class snw(Fix):
    """Fixes for snw."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 100
        cube.metadata = metadata
        return cube


class snc(snw):
    """Fixes for snc."""

    # dayspermonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    # if ((name.eq."snc".or.name.eq."snw").and.FIELD.eq."T2Ds".and. \
    #     ENSEMBLE.eq."r1i1p1") then
    #     opt = 0
    #     opt@calendar = var&time@calendar
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


class msftmyz(Fix):
    """Fixes for msftmyz."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes mask

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube.data = np.ma.array(cube.data)
        cube.data = np.ma.masked_where(cube.data.mask + (cube.data == 0.),
                                       cube.data)

        return cube
