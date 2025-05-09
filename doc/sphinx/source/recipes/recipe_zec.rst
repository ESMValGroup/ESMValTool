.. _recipes_zec:

Zero Emissions Commitment (ZEC)
===============================

Overview
--------

The Zero Emissions Commitment (ZEC) quantifies the change in global mean
temperature expected to occur after net carbon dioxide (CO\ :sub:`2`)
emissions cease. ZEC is therefore important to consider when estimating
the remaining carbon budget. Calculation of ZEC requires dedicated simulations
with emissions set to zero, branching off a base simulation with emissions.
In CMIP6 the simulations were part of ZECMIP, with the simulations called
``esm-1pct-brch-xPgC`` branching off the ``1pctCO2`` simulation when emissions
reach x PgC. The default x was 1000PgC, with additional simulations for 750PgC
and 2000PgC. In CMIP7, ZEC simulations (``esm-flat10-zec``) are part of the
fast track and branch off (``esm-flat10``) with constant emissions of 10GtC/yr
at year 100 (Sanderson 2024).


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_zec.yml


Diagnostics are stored in diag_scripts/

   * :ref:`climate_metrics/zec.py <api.esmvaltool.diag_scripts.climate_metrics.zec>`


User settings in recipe
-----------------------

*Note*: The preprocessor settings should not be changed. When changing the
recipe, ensure that the reference simulation period (generally ``1pctCO2``
or ``esm-flat10``, using the ``anomaly_base`` preprocessor block) is the
20-year period around the time the ZEC simulation branches off of.
This point may vary in CMIP6 for 1pctCO2, while in CMIP7 it should be at
year 100 of ``esm-flat10``.

* Preprocessor

   * ``area_statistics`` (*operation: mean*): Calculate global mean.
   * ``annual_statistics`` (*operation: mean*): For the ZEC simulation,
     compute annual statistics.
   * ``climate_statistics`` (*operation: mean, period: full*): For the
     reference simulation, compute the 20-year time average centered around
     the time when the ZEC simulation starts.

.. _tcr.py:

* Script climate_metrics/zec.py

   * ``zec_year``, *list*, optional (default: ``[50]``): Calculate ZEC for
     the 20-year average centered around year x, multiple values are possible.
     Barplots are generated for all ZEC\ :sub:`x`.
   * ``experiments``, *dict*, optional: When using non-default experiments
     to calculate ZEC, the experiment setting is required with values for
     ``reference`` and ``simulation``. Accepted by default are:
     ``reference: esm-flat10, 1pctCO2``
     ``simulation: esm-flat10-zec, esm-1pct-brch-1000PgC``


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* MacDougall, A. H., Frölicher, T. L., Jones, C. D., Rogelj, J., Matthews,
  H. D., Zickfeld, K., Arora, V. K., Barrett, N. J., Brovkin, V., Burger,
  F. A., Eby, M., Eliseev, A. V., Hajima, T., Holden, P. B., Jeltsch-Thömmes,
  A., Koven, C., Mengis, N., Menviel, L., Michou, M., Mokhov, I. I., Oka, A.,
  Schwinger, J., Séférian, R., Shaffer, G., Sokolov, A., Tachiiri, K.,
  Tjiputra, J., Wiltshire, A., and Ziehn, T.: Is there warming in the
  pipeline? A multi-model analysis of the Zero Emissions Commitment from
  CO\ :sub:`2`, Biogeosciences, 17, 2987–3016,
  https://doi.org/10.5194/bg-17-2987-2020, 2020.
* Sanderson, B. M., Brovkin, V., Fisher, R., Hohn, D., Ilyina, T., Jones, C.,
  Koenigk, T., Koven, C., Li, H., Lawrence, D., Lawrence, P., Liddicoat, S.,
  Macdougall, A., Mengis, N., Nicholls, Z., O'Rourke, E., Romanou, A.,
  Sandstad, M., Schwinger, J., Seferian, R., Sentman, L., Simpson, I., Smith,
  C., Steinert, N., Swann, A., Tjiputra, J., and Ziehn, T.: flat10MIP: An
  emissions-driven experiment to diagnose the climate response to positive,
  zero, and negative CO\ :sub:`2` emissions, EGUsphere [preprint],
  https://doi.org/10.5194/egusphere-2024-3356, 2024.

Example plots
-------------

.. _fig_zec_1:
.. figure:: /recipes/figures/zec/zec_ts.png
   :align: center
   :width: 50%

   Time series of ZEC - temperature change after cessation of emissions.

.. _fig_zec_2:
.. figure:: /recipes/figures/zec/zec_bar.png
   :align: center
   :width: 50%

   Barplot for ascending values of Zec\ :sub:`x`, with x = 50 in this case,
   symbolizing the Zero Emissions Commitment 50 years after emissions cease.
