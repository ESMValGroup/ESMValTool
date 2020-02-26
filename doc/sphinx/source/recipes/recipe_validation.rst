.. _recipes_validation:

Zonal and Meridional Means
==========================

Overview
--------

This functional diagnostic takes two models designated by CONTROL and EXPERIMENT and compares them via a number of
analyses. Optionally a number of observational datasets can be added for processing. There are three types of standard analysis:
lat_lon, meridional_mean and zonal_mean. Each of these diagnostics can be run on a separate basis (each an entry to diagnostics/scripts).
The lat_lon analysis produces the following plots: a simple global plot for each variable for each dataset, a global plot for the
difference between CONTROL and EXPERIMENT, a global plot for the difference between CONTROL and each of the observational datasets.
The meridional_mean and zonal_mean produce variable vs coordinate (``latitude`` or ``longitude``) with both ``CONTROL`` and ``EXPERIMENT`` curves
in each plot, for the entire duration of time specified and also, if the user wishes, for each season (seasonal means): winter DJF, spring MAM, summer JJA, autumn SON (by setting ``seasonal_analysis: true`` in the recipe).

At least regridding on a common grid for all model and observational datasets should be performed in preprocessing (if datasets
are on different grids). Also note that currently it is not allowed to use the same dataset (with varying parameters like experiment
or ensemble) for both CONTROL and EXPERIMENT (the use case for comparison between different experiments or ensembles for the same model
will be implemented in a future release).

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_validation.yml (CMIP5)
* recipe_validation_CMIP6.yml (CMIP6)

Diagnostics are stored in diag_scripts/

* validation.py
* shared/_validation.py

User settings
-------------

#. validation.py

   *Required settings for script*

   * title: title of the analysis, user defined;
   * control_model: control dataset name e.g. UKESM1-0-LL;
   * exper_model: experiment dataset name e.g. IPSL-CM6A-LR;
   * observational_datasets: list of at least one element; if no OBS wanted comment out; e.g. ['ERA-Interim'];
   * analysis_type: use any of: lat_lon, meridional_mean, zonal_mean;
   * seasonal_analysis: boolean, if seasonal means are needed e.g. ``true``;

Variables
---------

* any variable

Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips or OBS or ana4mips can be used.*

* any observations
* it is important to note that all observational data should go through the same preprocessing as model data

References
----------

* none, basic technical analysis

Example plots
-------------

.. figure:: /recipes/figures/validation/Merid_Mean_DJF_longitude_tas_UKESM1-0-LL_vs_IPSL-CM6A-LR.png
   :width: 70 %

   Meridional seasonal mean for winter (DJF) comparison beween CMIP6 UKESM1 and IPSL models.

.. figure:: /recipes/figures/validation/Zonal_Mean_DJF_latitude_tas_UKESM1-0-LL_vs_IPSL-CM6A-LR.png
   :width: 70 %

   Zonal seasonal mean for winter (DJF) comparison beween CMIP6 UKESM1 and IPSL models.
