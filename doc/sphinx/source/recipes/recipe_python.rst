.. _recipe_python:

Example recipe Python
=====================

Overview
--------

This is an example recipe calling a simple diagnostic script written in Python.
The recipe produces time series plots of global mean temperature and for the
temperature in Amsterdam. It also produces a map of global temperature in
January 2020.

For detailed instructions on obtaining input data, please refer to
:ref:`_inputdata`. However, in case you just quickly want to run through the
example, you can use the following links to obtain the data from ESGF:

  * `BCC-ESM1 <http://cmip.bcc.cma.cn/thredds/fileServer/cmip6_data/CMIP/BCC/BCC-ESM1/historical/r1i1p1f1/Amon/tas/gn/v20181214/tas_Amon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc>`_

  * `CanESM2 <http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esg_dataroot/AR5/CMIP5/output/CCCma/CanESM2/historical/mon/atmos/tas/r1i1p1/tas_Amon_CanESM2_historical_r1i1p1_185001-200512.nc>`_

Please refer to the terms of use for `CMIP5
<https://pcmdi.llnl.gov/mips/cmip5/terms-of-use.html>`_ and `CMIP6
<https://pcmdi.llnl.gov/CMIP6/TermsOfUse/TermsOfUse6-1.html>`_ data.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * examples/recipe_python.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * examples/diagnostic.py: visualize cubes and store provenance information


User settings in recipe
-----------------------

#. Script ``examples/diagnostic.py``

   *Required settings for script*

   * ``quickplot: plot_type``: which of `Iris' quickplot <https://scitools.org.uk/iris/docs/latest/iris/iris/quickplot.html>`_ functions to use.


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
