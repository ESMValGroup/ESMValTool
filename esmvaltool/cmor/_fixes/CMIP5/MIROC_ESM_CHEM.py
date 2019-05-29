# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for MIROC ESM CHEM"""
from ..fix import Fix


class tro3(Fix):
    """Fixes for tro3"""

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
        cube *= 1000
        cube.metadata = metadata
        return cube


# if (name .eq. "tro3") then
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
