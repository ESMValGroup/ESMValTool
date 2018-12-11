"""Fixes for GCP."""
import iris
import numpy as np

from ..fix import Fix


class nbp_global(Fix):
    """Class to fix nbp_global."""

    def fix_file(self, filepath, output_dir):
        """Apply fixes to the files prior to creating the cube.

        Fix `var_name`, `standard_name` and `long_name`.

        """
        standard_name = ('surface_net_downward_mass_flux_of_carbon_dioxide_'
                         'expressed_as_carbon_due_to_all_land_processes_'
                         'globally_averaged')
        cube = iris.load_cube(filepath)
        iris.std_names.STD_NAMES[standard_name] = {
            'canonical_units': cube.units,
        }
        cube.var_name = 'nbp_global'
        cube.standard_name = standard_name
        cube.long_name = ('Carbon Mass Flux out of Atmosphere due to Net '
                          'Biospheric Production on Land (globally averaged)')
        new_path = Fix.get_fixed_filepath(output_dir, filepath)
        iris.save(cube, new_path)
        return new_path

    # def fix_metadata(self, cube):
    #     """
    #     Fix metadata for npb.

    #     Add latitude and longitude coordinates as GCP data is actually T0M.

    #     FIXME:
    #         Better introduce a new variable (e.g. 'nbp_t0m'), but at the
    #         moment NCL is not able to process multiple variables (#531).
    #         Best solution right now: rename "T0M" -> "T2Ms" in the filename.
    #     """
    #     time_coord = cube.coord('time')
    #     lat_coord = iris.coords.DimCoord([0],
    #                                      standard_name='latitude',
    #                                      long_name='latitude coordinate',
    #                                      var_name='lat',
    #                                      units='degrees_north',
    #                                      bounds=[-90, 90])
    #     lon_coord = iris.coords.DimCoord([180],
    #                                      standard_name='longitude',
    #                                      long_name='longitude coordinate',
    #                                      var_name='lon',
    #                                      units='degrees_east',
    #                                      bounds=[0, 360])
    #     metadata = cube.metadata

    #     # Create new cube with latitude and longitude coordinate
    #     new_data = cube.data[:, np.newaxis, np.newaxis]
    #     cube = iris.cube.Cube(new_data,
    #                           dim_coords_and_dims=[(time_coord, 0),
    #                                                (lat_coord, 1),
    #                                                (lon_coord, 2)],
    #                           **metadata._asdict())
    #     cube.attributes['field'] = "T2M"

    #     return cube
