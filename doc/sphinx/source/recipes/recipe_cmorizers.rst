.. _recipe_cmorizers:

CMORizer recipes
=================

Overview
--------

These are CMORizer recipes calling CMORizer diagnostic scripts.

ESMValCore supports ERA5 hourly and monthly datasets in their native
format, see :ref:`inputdata_native_datasets`.
and `ERA5 data documentation <https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation>`_.
It may be useful in some cases to create ERA5 daily CMORized data. This can be
achieved by using a CMORizer *recipe*,
see `recipe_daily_era5.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/cmorizers/recipe_daily_era5.yml>`_.
This recipe reads native, hourly ERA5 data, performs a daily aggregation
preprocessor, and then calls a diagnostic that operates on the data. In this
example, the diagnostic renames the files to the standard OBS6 file names. The output
are thus daily, CMORized ERA5 data, that can be used through the OBS6 project.
As such, this example recipe creates a local pool of CMORized data. The advantage, in this
case, is that the daily aggregation is performed only once, which can save a lot
of time and compute if it is used often.

The example CMORizer recipe can be run like any other ESMValTool recipe:

.. code-block:: bash

    esmvaltool run cmorizers/recipe_daily_era5.yml

Note that the ``recipe_daily_era5.yml`` adds the next day of the new year to
the input data. This is because one of the fixes needed for the ERA5 data is to
shift the time axis of non-instantaneous variables half an hour back in time, resulting in a missing
record on the last day of the year. ERA5 data can be downloaded using `era5cli <https://era5cli.readthedocs.io>`_.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * cmorizers/recipe_daily_era5.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * cmorizers/era5.py: generates output filename


User settings in recipe
-----------------------

#. cmorizers/recipe_daily_era5.yml

   *Required add_one_day preprocessor settings:*

   * start_year: 1990
   * start_month: 1
   * start_day: 1
   * end_year: 1991
   * end_month: 1
   * end_day: 1

These settings should not be changed
   * daily_mean:
         operator: mean
   * daily_min:
         operator: min
   * daily_max:
         operator: max

Variables
---------

#. cmorizers/recipe_daily_era5.yml

   * clt
   * evspsbl
   * evspsblpot
   * mrro
   * pr
   * prsn
   * ps
   * psl
   * rlds
   * rls
   * rsds
   * rsdt
   * rss
   * tas
   * tasmax
   * tasmin
   * tdps
   * ts
   * tsn
   * uas
   * vas

References
----------

* Hersbach, H., et al., Quarterly Journal of the Royal Meteorological Society, 730, 1999-2049, doi:10.1002/qj.3803, 2020.
