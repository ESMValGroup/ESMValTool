.. _recipes_ecs:

Constraining future Indian Summer Monsoon projections with  the present-day precipitation over the tropical western Pacific
===========================================================================================================================

Overview
--------


Following `Li et al. (2017)`_ the change between present-day and future Indian Summer Monsoon (ISM) precipitation is constrained
using the precipitation over the tropical western Pacific compared to
a fixed, observed amound of 6 mm d-1 from Global Precipitation Climatology Project (GPCP) for 1980-2009 `(Adler et al., 2003)`_.
For CMIP6 historical data for 1980-2009 should be used. For CMIP5 historical data only from 1980-2005 should be used, due to the length of the data sets.
At the moment it is not possiible to use a combined [historical, rcp] data set, because the diagnostic requires that a historical data set exists.

.. _`(Adler et al., 2003)`: https://journals.ametsoc.org/doi/abs/10.1175/1525-7541%282003%29004%3C1147%3ATVGPCP%3E2.0.CO%3B2
.. _`Li et al. (2017)`: https://www.nature.com/articles/nclimate3387


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_li17natcc.yml
   * recipe_li17natcc_cmip6.yml


Diagnostics are stored in diag_scripts/

   * emergent_constraints/lif1f2.py


User settings in recipe
-----------------------

#. Preprocessor

   * ``area_statistics`` (*operation: mean*): Calculate global mean.

#. Script climate_metrics/ecs.py

   * ``calculate_mmm``, *bool*, optional (default: ``True``): Calculate
     multi-model mean ECS.
   * ``read_external_file``, *str*, optional: Read ECS and net climate feedback
     parameter from external file. Can be given relative to the diagnostic
     script or as absolute path.
   * ``seaborn_settings``, *dict*, optional: Options for seaborn's ``set()``
     method (affects all plots), see
     https://seaborn.pydata.org/generated/seaborn.set.html.

#. Script climate_metrics/create_barplot.py

   * ``label_attribute``, *str*, optional: Attribute of the cube which is used
     as label for the different input files in the barplot.
   * ``patterns``, *list of str*, optional: Patterns to filter list of input
     files.
   * ``seaborn_settings``, *dict*, optional: Options for seaborn's ``set()``
     method (affects all plots), see
     https://seaborn.pydata.org/generated/seaborn.set.html.
   * ``sort_ascending``, *bool*, optional (default: ``False``): Sort bars in
     ascending order.
   * ``sort_descending``, *bool*, optional (default: ``False``): Sort bars in
     descending order.
   * ``value_labels``, *bool*, optional (default: ``False``): Label bars with
     value of that bar.
   * ``y_range``, *list of float*, optional: Range for the Y axis of the plot.

#. Script climate_metrics/create_scatterplot.py

   * ``dataset_style``, *str*, optional: Name of the style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
   * ``pattern``, *str*, optional: Pattern to filter list of input files.
   * ``seaborn_settings``, *dict*, optional: Options for seaborn's ``set()``
     method (affects all plots), see
     https://seaborn.pydata.org/generated/seaborn.set.html.
   * ``y_range``, *list of float*, optional: Range for the Y axis of the plot.


Variables
---------

* *pr* (atmos, monthly, longitude, latitude, time)
* *ua* (atmos, monthly, longitude, latitude, pressure level, time)
* *va* (atmos, monthly, longitude, latitude, pressure level, time)
* *ts* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Li, G., Xie, S. P., He, C., and Chen, Z. S.: Western Pacific emergent constraint lowers projected increase in Indian summer monsoon rainfall, Nat Clim Change, 7, 708-+, 2017


Example plots
-------------

.. _fig_ecs_1:
.. figure:: /recipes/figures/ecs/CanESM2.png
   :align: center
   :width: 50%

   Scatterplot between TOA radiance and global mean surface temperature anomaly
   for 150 years of the abrupt 4x CO2 experiment including linear regression to
   calculate ECS for CanESM2 (CMIP5).
