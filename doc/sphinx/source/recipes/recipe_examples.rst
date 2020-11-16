.. _recipe_examples:

Example recipes
===============

Overview
--------

These are example recipes calling example diagnostic scripts.

The recipe examples/recipe_python.yml produces time series plots of global mean
temperature and for the temperature in Amsterdam.
It also produces a map of global temperature in January 2020.

The recipe examples/recipe_extract_shape.yml produces a map of the mean
temperature in the Elbe catchment over the years 2000 to 2002.
Some example shapefiles for use with this recipe are available
`here <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/diag_scripts/shapeselect/testdata>`__,
make sure to download all files with the same name but different extensions.

For detailed instructions on obtaining input data, please refer to
:ref:`inputdata`. However, in case you just quickly want to run through the
example, you can use the following links to obtain the data from ESGF:

  * `BCC-ESM1 <http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/BCC/BCC-ESM1/historical/r1i1p1f1/Amon/tas/gn/v20181214/tas_Amon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc>`_
  * `CanESM2 <http://esgf2.dkrz.de/thredds/fileServer/lta_dataroot/cmip5/output1/CCCma/CanESM2/historical/mon/atmos/Amon/r1i1p1/v20120718/tas/tas_Amon_CanESM2_historical_r1i1p1_185001-200512.nc>`_

Please refer to the terms of use for `CMIP5
<https://pcmdi.llnl.gov/mips/cmip5/terms-of-use.html>`_ and `CMIP6
<https://pcmdi.llnl.gov/CMIP6/TermsOfUse/TermsOfUse6-1.html>`_ data.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * examples/recipe_python.yml
    * examples/recipe_extract_shape.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * examples/diagnostic.py: visualize results and store provenance information


User settings in recipe
-----------------------

#. Script ``examples/diagnostic.py``

   *Required settings for script*

   * ``quickplot: plot_type``: which of `Iris' quickplot <https://scitools.org.uk/iris/docs/latest/iris/iris/quickplot.html>`_ functions to use.
     Arguments that are accepted by these functions can also be specified here, e.g. ``cmap``.
     Preprocessors need to be configured such that the resulting data matches the plot type, e.g. a timeseries or a map.

Variables
---------

* tas (atmos, monthly, longitude, latitude, time)

Example plots
-------------

.. _global_map:
.. figure::  /recipes/figures/examples/map.png
   :align:   center

   Air temperature in January 2000 (BCC-ESM1 CMIP6).

.. _timeseries:
.. figure::  /recipes/figures/examples/timeseries.png
   :align:   center

   Amsterdam air temperature (multimodel mean of CMIP5 CanESM2 and CMIP6 BCC-ESM1).

.. _elbe:
.. figure::  /recipes/figures/examples/elbe.png
   :align:   center

   Mean air temperature over the Elbe catchment during 2000-2002 according to CMIP5 CanESM2.
