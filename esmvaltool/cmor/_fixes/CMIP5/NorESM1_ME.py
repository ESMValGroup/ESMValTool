# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for NorESM1-ME model."""

import numpy as np

from ..fix import Fix


class tas(Fix):
    """Fixes for tas."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Some coordinate points vary for different files of this dataset (for
        different time range). This fix removes these inaccuracies by rounding
        the coordinates.

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        for cube in cubes:
            for coord in cube.coords(dim_coords=True):
                for attr in ('points', 'bounds'):
                    setattr(coord, attr, np.round(getattr(coord, attr), 12))
        return cubes
