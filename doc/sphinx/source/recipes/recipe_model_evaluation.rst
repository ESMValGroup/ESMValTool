.. _recipe_model_evaluation:

General model evaluation
========================

Overview
--------

These recipes and diagnostics provide a basic climate model evaluation with
observational data.
This is especially useful to get an overview of the performance of a
simulation.
The diagnostics used here allow plotting arbitrary preprocessor output, i.e.,
arbitrary variables from arbitrary datasets.


Available recipes and diagnostics
---------------------------------

Recipes are stored in `recipes/model_evaluation`

* recipe_model_evaluation_basics.yml
* recipe_model_evaluation_clouds_clim.yml
* recipe_model_evaluation_clouds_cycles.yml
* recipe_model_evaluation_precip_zonal.yml

Diagnostics are stored in `diag_scripts/monitor/`

* :ref:`multi_datasets.py
  <api.esmvaltool.diag_scripts.monitor.multi_datasets>`:
  Monitoring diagnostic to show multiple datasets in one plot (incl. biases).


User settings
-------------

It is recommended to use a vector graphic file type (e.g., SVG) for the output
format when running this recipe, i.e., run the recipe with the
:ref:`configuration options <esmvalcore:config_options>` ``output_file_type:
svg``.
Note that map and profile plots are rasterized by default.
Use ``rasterize: false`` in the recipe to disable
this.


Recipe settings
~~~~~~~~~~~~~~~

A list of all possible configuration options that can be specified in the
recipe is given for each diagnostic individually (see links given for the
available diagnostics in the previous section).


Variables
---------

Any, but the variables' number of dimensions should match the ones expected by
each diagnostic (see links given for the available diagnostics in the previous
section).


Example plots
-------------

.. figure::  /recipes/figures/model_evaluation/map_tas_MPI-ESM1-2-HR_Amon.jpg
   :align:   center
   :width:   14cm

   Global climatology of 2m near-surface air temperature.

.. figure::  /recipes/figures/model_evaluation/map_swcre_MPI-ESM1-2-HR_Amon.jpg
   :align:   center
   :width:   14cm

   Global climatology of the shortwave cloud radiative effect (SWCRE).

.. figure::  /recipes/figures/model_evaluation/timeseries_rtnt_ambiguous_dataset_Amon.jpg
   :align:   center
   :width:   14cm

   Time series of the global mean top-of-the-atmosphere net radiative flux.

.. figure::  /recipes/figures/model_evaluation/variable_vs_lat_pr_Amon.jpg
   :align:   center
   :width:   14cm

   Zonal mean precipitation.

.. figure::  /recipes/figures/model_evaluation/annual_cycle_clt_southerocean_Amon.jpg
   :align:   center
   :width:   14cm

   Annual cycle of Southern Ocean total cloud cover.

.. figure::  /recipes/figures/model_evaluation/diurnal_cycle_clt_sepacific_3hr.png
   :align:   center
   :width:   14cm

   Diurnal cycle of Southeast Pacific total cloud cover.
