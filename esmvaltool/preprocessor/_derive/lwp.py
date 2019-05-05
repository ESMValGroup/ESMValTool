"""Derivation of variable `lwp`."""

import logging

from iris import Constraint

from ._baseclass import DerivedVariableBase

logger = logging.getLogger(__name__)


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `lwp`."""

    # Required variables
    required = [
        {
            'short_name': 'clwvi'
        },
        {
            'short_name': 'clivi'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute liquid water path.

        Note
        ----
        Some datasets output the variable `clwvi` which only contains `lwp`. In
        these cases, the input `clwvi` cube is just returned.

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
            'CESM1-CAM5-1-FV2',
            'CESM1-CAM5',
            'CMCC-CESM',
            'CMCC-CM',
            'CMCC-CMS',
            'IPSL-CM5A-MR',
            'IPSL-CM5A-LR',
            'IPSL-CM5B-LR',
            'CCSM4',
            'IPSL-CM5A-MR',
            'MIROC-ESM',
            'MIROC-ESM-CHEM',
            'MIROC-ESM',
            'CSIRO-Mk3-6-0',
            'MPI-ESM-MR',
            'MPI-ESM-LR',
            'MPI-ESM-P',
        ]
        if (project in ["CMIP5", "CMIP5_ETHZ"] and dataset in bad_datasets):
            logger.info(
                "Assuming that variable clwvi from %s dataset %s "
                "contains only liquid water", project, dataset)
            lwp_cube = clwvi_cube
        else:
            lwp_cube = clwvi_cube - clivi_cube

        return lwp_cube
