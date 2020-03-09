.. _recipes_deangelis15nat:

Evaluate water vapor short wave radiance absorption schemes of ESMs with the observations.
==========================================================================================================================

Overview
--------


The recipe reproduces figures from `DeAngelis et al. (2015)`_:
Figure 1b to 4 from the main part as well as extended data figure 1 and 2.
This paper compares models with different schemes for water vapor short wave radiance absorption with the observations.
Schemes using pseudo-k-distributions with more than 20 exponential terms show the best results.

.. _`DeAngelis et al. (2015)`: https://www.nature.com/articles/nature15770


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_deangelis15nat.yml

Diagnostics are stored in diag_scripts/

   * deangelis15nat/deangelisf1b.py
   * deangelis15nat/deangelisf2ext.py
   * deangelis15nat/deangelisf3f4.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models.
deangelisf1b.py:

deangelisf2ext.py:

deangelisf3f4.py:
For each model, two experiments must be given:
a pre industrial control run, and a scenario with 4 times CO\ :sub:`2`\.
Possibly, 150 years should be given, but shorter time series work as well.


Variables
---------

deangelisf1b.py:

deangelisf2ext.py:

deangelisf3f4.py:
* *rsnstcs* (atmos, monthly, longitude, latitude, time)
* *rsnstcsnorm* (atmos, monthly, longitude, latitude, time)
* *prw* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

deangelisf3f4.py:
* *rsnstcs*
CERES-EBAF

* *prw*
ERA-Interim, SSMI

References
----------

* DeAngelis, A. M., Qu, X., Zelinka, M. D., and Hall, A.: An observational radiative constraint on hydrologic cycle intensification, Nature, 528, 249, 2015.


Example plots
-------------
