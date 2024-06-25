.. _recipes_seaborn_diag:

Seaborn Diagnostics
===================

Overview
--------

These recipes showcase the use of the Seaborn diagnostic that provides a
high-level interface to `Seaborn <https://seaborn.pydata.org>`__ for ESMValTool
recipes.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_seaborn.yml

Diagnostics are stored in diag_scripts/

* :ref:`seaborn_diag.py <api.esmvaltool.diag_scripts.seaborn_diag>`


Variables
---------

Arbitrary variables are supported.


Observations and reformat scripts
---------------------------------

Arbitrary datasets are supported.


References
----------

* Waskom, M. L. (2021), seaborn: statistical data visualization, Journal of
  Open Source Software, 6(60), 3021, doi:10.21105/joss.03021.


Example plots
-------------

.. _fig_seaborn_1:
.. figure:: /recipes/figures/seaborn/ta_vs_lat.jpg
   :align: center
   :width: 50%

   Monthly and zonal mean temperatures vs. latitude in the period 1991-2014 for
   two Earth system models (CESM2-WACCM and GFDL-ESM4).
   Colors visualize the corresponding pressure levels.

.. _fig_seaborn_2:
.. figure:: /recipes/figures/seaborn/regional_pr_hists.jpg
   :align: center
   :width: 50%

   Spatiotemporal distribution of daily precipitation in the period 2005-2014
   for six IPCC AR6 regions simulated by two Earth system models (CESM2-WACCM
   and GFDL-ESM4).
   Each day in each grid cell in the corresponding regions is considered with
   equal weight.
