
.. _recipes_martin18grl:

Agricultural Droughts in CMIP6 Future Projections
=================================================

Overview
--------

The two recipes presented here evaulate historical simulations of 18 CMIP6
models and analyse their projections for three different future pathways.
The results are published in Lindenlaub (2025).


Available recipes and diagnostics
---------------------------------

Recipes are stored in ``recipes/droughts/``

   * recipe_lindenlaub25_historical.yml
   * recipe_lindenlaub25_scenarios.yml

Diagnostics used by this recipes:

   * droughts/pet.R
   * droughts/spei.R
   * :ref:`droughts/diffmap.py <api.esmvaltool.diag_scripts.droughts.diffmap>`
   * :ref:`droughts/distribution.py <api.esmvaltool.diag_scripts.droughts.distribution>`
   * :ref:`droughts/event_area_timeseries.py <api.esmvaltool.diag_scripts.droughts.event_area_timeseries>`
   * :ref:`droughts/pattern_correlation.py <api.esmvaltool.diag_scripts.droughts.pattern_correlation>`
   * :ref:`droughts/timeseries_scenarios.py <api.esmvaltool.diag_scripts.droughts.timeseries_scenarios>`
   * :ref:`droughts/regional_hexagons.py <api.esmvaltool.diag_scripts.droughts.regional_hexagons>`

Data
----

Soil moisture is evaluated and discussed, but not required for PET and SPEI
calculation. 
``tasmin``, ``tasmax``, ``sfcWind``, ``ps``, ``rsds`` are used to approximate 
``evspsblpot`` for ERA5 and 18 CMIP6 datasets.
``SPEI`` is calculated from ``evspsblpot`` and ``pr``. 

References
----------

* Lindenlaub, L. (2025). Agricultural Droughts in CMIP6 Future Projections.
  Journal of Climate, 38(1), 1-15. https://doi.org/10.1029/2025JC012345
