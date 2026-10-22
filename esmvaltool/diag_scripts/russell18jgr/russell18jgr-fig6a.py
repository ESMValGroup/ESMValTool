"""Density-layer volume transport across 30S (russell18jgr figure 6a).

Python port of russell18jgr-fig6a.ncl.

Plots the volume transport (Sv) across 30S in the density layers of
Talley (2008): dark blue bars are the layer totals (compare with the
magenta observed values), narrow red bars are equal subdivisions of
each layer.  See russell_fig6_shared.py for the implementation.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import russell_fig6_shared

from esmvaltool.diag_scripts.shared import run_diagnostic

if __name__ == "__main__":
    with run_diagnostic() as config:
        russell_fig6_shared.run_fig6(config, mode="volume")
