"""Derivation of variable `lwp`."""


import logging

from iris import Constraint

from .derived_variable import DerivedVariable

logger = logging.getLogger(__name__)


class lwp(DerivedVariable):  # noqa
    """Derivation of variable `lwp`."""

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
        return [('clwvi', 'T2' + frequency + 's'),
                ('clivi', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute liquid water path.

        Liquid water path is calculated by subtracting `clivi` (ice water) from
        `clwvi` (condensed water path).
        Note: Some datasets output the variable `clwvi` which only contains
        lwp. In these cases, the input `clwvi` cube is just returned.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `clwvi`
            (`atmosphere_cloud_condensed_water_content`) and `clivi`
            (`atmosphere_cloud_ice_content`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing liquid water path.

        """
        clwvi_cube = cubes.extract_strict(
            Constraint(name='atmosphere_cloud_condensed_water_content'))
        clivi_cube = cubes.extract_strict(
            Constraint(name='atmosphere_cloud_ice_content'))

        dataset = clwvi_cube.attributes.get('model_id')
        project = clwvi_cube.attributes.get('project_id')
        # Should we check that the model_id/project_id are the same on both
        # cubes?

        bad_datasets = [
            'CESM1-CAM5-1-FV2', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM',
            'CMCC-CMS', 'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR',
            'CCSM4', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC-ESM-CHEM',
            'MIROC-ESM', 'CSIRO-Mk3-6-0', 'MPI-ESM-MR', 'MPI-ESM-LR',
            'MPI-ESM-P',
        ]
        if ((project in ["CMIP5", "CMIP5_ETHZ"] and dataset in bad_datasets)
                or (project == 'OBS' and dataset == 'UWisc')):
            logger.info(
                "Assuming that variable clwvi from %s dataset %s "
                "contains only liquid water", project, dataset)
            lwp_cube = clwvi_cube
        else:
            lwp_cube = clwvi_cube - clivi_cube

        return lwp_cube
