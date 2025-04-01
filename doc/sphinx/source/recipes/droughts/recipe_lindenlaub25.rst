
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


TODO: ERA5 native on the fly or cmorizer?

Reference Data (ERA5, CDS-SM, CRU) can be downloaded and cmorized by the
esmvaltool executing the `data` command:

```
esmvaltool data download ERA5 pr
esmvaltool data format ERA5 pr
```

``tasmin`` and ``tasmax`` are not directly available in the monthly averaged data
product. The download script `diag_scripts/droughts/download_era5_tasminmax.py`
can be used to download and preprocess the data based on the
``minimum_2m_temperature_since_previous_post_processing`` and
``maximum_2m_temperature_since_previous_post_processing`` variables from ERA5
daily data on single levels. The output files are compatible with the esmvaltool
and can be copied into the ``native6/Tier3/ERA5/v1/mon/``
directory. Using the esmvaltool environment ensures that all required libraries
are available. The script can be run with the following command:

```
python diag_scripts/droughts/download_era5_tasminmax.py
```

For more options use the ``--help`` flag.

The CMIP6 data can be downloaded automatically by the ESMValTool. Just ensure
that ``esgf_download`` is set to ``True`` or ``when_missing`` in the 
user configuration.

Figures
-------

References
----------

* Lindenlaub, L. (2025). Agricultural Droughts in CMIP6 Future Projections.
  Journal of Climate, 38(1), 1-15. https://doi.org/10.1029/2025JC012345
