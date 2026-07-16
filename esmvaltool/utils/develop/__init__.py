"""Interface to utility commands for develop command group."""

import sys
from pathlib import Path

from esmvaltool.utils.develop.compare import compare


class DevelopCommand:
    """Development utilities."""

    def compare(self, reference_dir, current_dir, verbose=False):
        """Compare a recipe run to a reference run.

        Returns True if the runs were identical, False otherwise.

        Parameters
        ----------
        reference_dir : str
            Results directory from reference run
        current_dir : str
            Results directory from run to be tested
        verbose : bool
            Display more information on differences
        """
        same = compare(Path(reference_dir), Path(current_dir), verbose)

        sys.exit(int(not same))
