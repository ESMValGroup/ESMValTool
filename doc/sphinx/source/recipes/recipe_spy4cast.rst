.. _recipe_monitor:

Monitor
=======

Overview
--------

These recipes and diagnostics allow plotting arbitrary preprocessor output,
i.e., arbitrary variables from arbitrary datasets.
In addition, a :ref:`base class
<api.esmvaltool.diag_scripts.monitor.monitor_base>` is provided that allows a
convenient interface for all monitoring diagnostics.


Available recipes and diagnostics
---------------------------------

Recipes are stored in `recipes/monitor`

* recipe_mca.yml
* recipe_nn4cast.yml

Diagnostics are stored in `diag_scripts/monitor/`

* :ref:`mca.py <api.esmvaltool.diag_scripts.mca.py>`:
  Monitoring diagnostic to plot arbitrary preprocessor output.
* :ref:`diagnostic_nn4cast.py <api.esmvaltool.diag_scripts.diagnostic_nn4cast.py>`:
  Monitoring diagnostic to plot EOF maps and associated PC timeseries.



User settings
-------------

It is recommended to use a vector graphic file type (e.g., SVG) for the output
format when running this recipe, i.e., run the recipe with the
:ref:`configuration options <esmvalcore:config_options>` ``output_file_type:
svg``.
Note that map and profile plots are rasterized by default.
Use ``rasterize_maps: false`` or ``rasterize: false`` (see `Recipe settings`_)
in the recipe to disable this.

Recipe settings
~~~~~~~~~~~~~~~

A list of all possible configuration options that can be specified in the
recipe is given for each diagnostic individually (see previous section).



Variables
---------

Any combination of two variables can be used, but they must be on regular grids

Example plots
-------------

.. _fig_climglobal:
.. figure::  /recipes/figures/monitor/clim.png
   :align:   center
   :width:   14cm

   Global climatology of tas.

