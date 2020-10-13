.. _recipe_python:

Example recipe Python
=====================

Overview
--------

This is an example recipe calling a simple diagnostic script written in Python.
The recipe produces time series plots of global mean temperature and for the
temperature in Amsterdam. It also produces a map of global temperature in
January 2020.

Data are available e.g. via the `climate4impact portal <https://climate4impact.eu/impactportal/data/esgfsearch.jsp>`_. The links below will help you find them:

  * `BCC-ESM1 <https://climate4impact.eu/impactportal/data/esgfsearch.jsp#project=CMIP6&variable=tas&frequency=mon&experiment_id=historical&member_id=r1i1p1f1&source_id=BCC-ESM1&>`_

  * `CanESM2 <https://climate4impact.eu/impactportal/data/esgfsearch.jsp#project=CMIP5&variable=tas&experiment=historical&model=CanESM2&ensemble=r1i1p1&time_frequency=mon&>`_

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
