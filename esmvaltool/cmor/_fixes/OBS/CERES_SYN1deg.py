"""Fixes for CERES-SYN1deg model"""
from ..fix import Fix


class rlutcs(Fix):

    """Fixes for rlutcs"""
    def fix_metadata(self, cubes):
        self.get_cube_from_list(cubes).standard_name = 'toa_outgoing_longwave_flux_assuming_clear_sky'
        return cubes

class rsdt(Fix):

    """Fixes for rsdt"""
    def fix_metadata(self, cubes):
        self.get_cube_from_list(cubes).standard_name = 'toa_incoming_shortwave_flux'
        return cubes



