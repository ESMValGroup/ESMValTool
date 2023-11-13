.. _recipe_model_evaluation:

General model evaluation
========================

Overview
--------

These recipes and diagnostics provide a basic climate model evaluation with
observational data.
This is especially useful to get a rough idea about the performance of a
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
files when running this recipe, i.e., run the recipe with the command line
option ``--output_file_type=svg`` or use ``output_file_type: svg`` in your
:ref:`esmvalcore:user configuration file`.
Note that map and profile plots are rasterized by default.
Use ``rasterize: false`` in the recipe to disable
this.


Recipe settings
~~~~~~~~~~~~~~~

A list of all possible configuration options that can be specified in the
recipe is given for each diagnostic individually (see previous section).


Variables
---------

Any, but the variables' number of dimensions should match the ones expected by each plot.


Example plots
-------------

.. _fig_climglobal:
.. figure::  /recipes/figures/monitor/clim.png
   :align:   center
   :width:   14cm

Global climatology of tas.
