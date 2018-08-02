"""Fixes for ISCCP"""

from ..fix import Fix


class clisccp(Fix):
    """Class to fix clisccp"""

    def fix_metadata(self, cube):
        """
        Fix metadata for clisccp

        Fix cube name to:
        isccp_cloud_area_fraction
        """
        cube.rename('isccp_cloud_area_fraction')
        return cube
