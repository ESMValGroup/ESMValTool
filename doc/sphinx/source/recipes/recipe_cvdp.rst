.. _recipes_cvdp:

The Climate Variability Diagnostics Package (CVDP)
==================================================

Overview
--------
The Climate Variability Diagnostics Package (CVDP) developed by NCAR's Climate Analysis Section is an analysis tool that documents the major modes of climate variability in models and observations, including ENSO, Pacific Decadal Oscillation, Atlantic Multi-decadal Oscillation, Northern and Southern Annular Modes, North Atlantic Oscillation, Pacific North and South American teleconnection patterns. For details please refer to the [1] and [2].

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_cvdp.yml

Diagnostics are stored in diag_scripts/cvdp/

    * cvdp_wrapper.py

User settings in recipe
-----------------------

Currently, the recipe must be used with a single dataset entry.

Variables
---------

* ts (atmos, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4mips data for download
instructions.*


References
----------
[1] http://www.cesm.ucar.edu/working_groups/CVC/cvdp/

[2] https://github.com/NCAR/CVDP-ncl

Example plots
-------------

.. figure::  /recipes/figures/cvdp/nam.prreg.ann.png
   :align:   center

   Atmospheric Modes of Variability; pr (annual)
