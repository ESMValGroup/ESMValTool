.. ocean_diagnostics:


========================
Ocean Diagnostics README
========================


Overview
--------

This package contains several diagnostics which produce some simple figures.

While these are built to evaluate the ocean, they almost


Diagnostics
===========
Diagnostics are stored in: esmvaltool/diag_scripts/ocean

The following python modules are included in the ocean diagnostics package.
* diagnostic_maps.py
* diagnostic_maps_quad.py
* diagnostic_omz.py
* diagnostic_profiles.py
* diagnostic_seaice.py
* diagnostic_timeseries.py
* diagnostic_tools.py
* diagnostic_transects.py
* ocean_diagnostics.rst


diagnostic_maps.py
-----------------------
This diagnostics produces a spatial map from a NetCDF.
It requires the input netCDF to have the following dimensions. Either:
* A two dimensional file: latitude, longitude
* A three dimensional file: depth, latitude, longitude
In the case of a 3D netCDF file, this diagnostic produces a map for EVERY layer.
For this reason, we recommend extracting a small number of specific layers in
the preprocessor, using the `extract_layer` preprocessor.

This diagnostic also includes the optional arguments, `threshold` and
`thresholds`.
* threshold: a single float.
* thresholds: a list of floats.
Only one of these arguments should be provided. These two arguments produce a
second kind of diagnostic map plot: a contour map showing the spatial
distribution of the threshold value,  for each dataset. Alternatively, if the
thresholds argument is used instead of the threshold, the single-dataset contour
map shows the contours of all thresholds.
In addition to the single dataset contour, a multi-dataset contour map is also
produced.
In the case of the threshold value, a single mult-dataset contour figure is
produced. In the case of the thresholds list, a mult-dataset contour figure is
produced for each value in the thresholds list.


diagnostic_maps_quad.py
-----------------------


diagnostic_omz.py
-----------------------


diagnostic_profiles.py
-----------------------


diagnostic_timeseries.py
------------------------


diagnostic_tools.py
-----------------------


diagnostic_seaice.py
-----------------------


diagnostic_transects.py
-----------------------
Produces plots of transect. This is a


ocean_diagnostics.rst
-----------------------
This is the README file that you are currently reading!








Associated recipes
========================


Recipes are stored in: esmvaltool/recipes

The following recipes are known to use these diagnostics:
* recipe_OceanBGC.yml
.. * recipe_OxygenMinimumZones.yml
* recipe_OceanPhysics.yml
* recipe_OceanQuadMap.yml
* recipe_SeaIceExtent.yml




recipe_OceanPhysics.yml
-----------------------
