"""Test derivation of `fgco2_grid`."""
import mock

import esmvaltool.preprocessor._derive.fgco2_grid as fgco2_grid

CUBES = 'mocked cubes'
STD_NAME = 'surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon'


@mock.patch.object(fgco2_grid, 'grid_area_correction', autospec=True)
def test_fgco2_grid_calculation(mock_grid_area_correction):
    """Test calculation of `fgco2_grid."""
    derived_var = fgco2_grid.DerivedVariable()
    derived_var.calculate(CUBES)
    mock_grid_area_correction.assert_called_once_with(
        CUBES, STD_NAME, ocean_var=True)
