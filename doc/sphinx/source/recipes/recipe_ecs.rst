.. _recipes_ecs:

Equilibrium climate sensitivity
===============================

Overview
--------


Equilibrium climate sensitivity is defined as the change in global mean
temperature as a result of a doubling of the atmospheric CO\ :sub:`2`
concentration compared to pre-industrial times after the climate system has
reached a new equilibrium. This recipe uses a regression method based on
`Gregory et al. (2004)`_ to calculate it for several CMIP models.

.. _`Gregory et al. (2004)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003GL018747


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_ecs.yml


Diagnostics are stored in diag_scripts/

   * climate_metrics/ecs.py


User settings in recipe
-----------------------

#. Preprocessor

   * ``area_statistics`` (*operation: mean*): Calculate global mean.

#. Script climate_metrics/ecs.py

   * ``calculate_mmm``, *bool*, optional (default: ``True``): Calculate
     multi-model mean ECS.
   * ``read_external_file``, *str*, optional: Read ECS and net climate feedback
     parameter from external file. All other input data is ignored.
   * ``seaborn_settings``, *dict*, optional: Options for seaborn's ``set()``
     method (affects all plots), see
     https://seaborn.pydata.org/generated/seaborn.set.html.


Variables
---------

* *rlut* (atmos, monthly, longitude, latitude, time)
* *rsdt* (atmos, monthly, longitude, latitude, time)
* *rsut* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Gregory, Jonathan M., et al. "A new method for diagnosing radiative forcing
  and climate sensitivity." Geophysical research letters 31.3 (2004).


Example plots
-------------

.. _fig_ecs_1:
.. figure:: /recipes/figures/ecs/CanESM2.png
   :align: center
   :width: 50%

   Scatterplot between TOA radiance and global mean surface temperature anomaly
   for 150 years of the abrupt 4x CO2 experiment including linear regression to
   calculate ECS for CanESM2 (CMIP5).
