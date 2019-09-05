"""Test derivation of `nbp_grid`."""
import mock

import esmvaltool.preprocessor._derive.nbp_grid as nbp_grid

CUBES = 'mocked cubes'
STD_NAME = ('surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_'
            'carbon_due_to_all_land_processes')


@mock.patch.object(nbp_grid, 'grid_area_correction', autospec=True)
def test_nbp_grid_calculation(mock_grid_area_correction):
    """Test calculation of `nbp_grid."""
    derived_var = nbp_grid.DerivedVariable()
    derived_var.calculate(CUBES)
    mock_grid_area_correction.assert_called_once_with(CUBES, STD_NAME)
