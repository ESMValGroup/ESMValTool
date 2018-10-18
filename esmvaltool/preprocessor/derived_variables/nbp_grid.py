"""Derivation of variable `nbp_grid`."""


import iris
from iris import Constraint

from esmvaltool._config import get_config_user_file
from esmvaltool._data_finder import get_input_fx_filelist
from ._derived_variable import DerivedVariable


class nbp_grid(DerivedVariable):  # noqa
    """Derivation of variable `nbp_grid`."""

    def get_required(self, frequency):
        """Get variable `short_name` and `field` pairs required for derivation.

        Parameters
        ----------
        frequency : str
            Frequency of the desired derived variable.

        Returns
        -------
        list of tuples
            List of tuples (`short_name`, `field`) of all variables required
            for derivation.

        """
        return [('nbp', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute net biome production relative to grid cell area.

        By default, `nbp` is defined relative to land area. For easy spatial
        integration, the original quantity is multiplied by the land area
        fraction (`sftlf`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `nbp` (`surface_net_downward_mass_flux_of_
            carbon_dioxide_expressed_as_carbon_due_to_all_land_processes`)

        Returns
        -------
        iris.cube.Cube
            `Cube` containing net biome production relative to grid cell area.

        """
        nbp_cube = cubes.extract_strict(
            Constraint(name='surface_net_downward_mass_flux_of_carbon_dioxide_'
                            'expressed_as_carbon_due_to_all_land_processes'))

        # Get sftlf
        cfg = get_config_user_file()
        self.variable['fx_files'] = ['sftlf']
        fx_files_dict = get_input_fx_filelist(variable=self.variable,
                                              rootpath=cfg['rootpath'],
                                              drs=cfg['drs'])

        # Load sftlf if possible and correct nbp
        if fx_files_dict['sftlf']:
            sftlf_cube = iris.load_cube(fx_files_dict['sftlf'])
            nbp_cube.data = nbp_cube.data * sftlf_cube.data / 100.0
        return nbp_cube
