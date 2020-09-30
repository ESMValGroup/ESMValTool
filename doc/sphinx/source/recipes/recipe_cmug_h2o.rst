.. _recipes_deangelis15nat_cmug:

Evaluate water vapor short wave radiance absorption schemes of ESMs with the observations, including ESACCI data.
==========================================================================================================================

Overview
--------


The recipe reproduces figures 3 and 4 from `DeAngelis et al. (2015)`_:
See also doc/sphinx/source/recipes/recipe_deangelis15nat.rst 
This paper compares models with different schemes for water vapor short wave radiance absorption with the observations.
Schemes using pseudo-k-distributions with more than 20 exponential terms show the best results.

.. _`DeAngelis et al. (2015)`: https://www.nature.com/articles/nature15770


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_deangelis15nat_cmug.yml

Diagnostics are stored in diag_scripts/

   * deangelis15nat/deangelisf3f4.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models.

deangelisf3f4.py:
For each model, two experiments must be given:
a pre industrial control run, and a scenario with 4 times CO\ :sub:`2`\.
Possibly, 150 years should be given, but shorter time series work as well.
Currently, HOAPS data are incuded as place holder for expected ESACCI-WV data, type CDR-2:
Gridded monthly time series of TCWV in units of kg/m2 (corresponds to prw)
that cover the global land and ocean areas with a spatial resolution of 0.05° / 0.5° 
for the period July 2002 to December 2017.


Variables
---------

deangelisf3f4.py:
* *rsnstcs* (atmos, monthly, longitude, latitude, time)
* *rsnstcsnorm* (atmos, monthly, longitude, latitude, time)
* *prw* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

deangelisf1b.py:
* None

deangelisf2ext.py:
* None

deangelisf3f4.py:

* *rsnstcs*:
   CERES-EBAF

* *prw*
   HOAPS, planed for ESACCI-WV data, type CDR-2


References
----------

* DeAngelis, A. M., Qu, X., Zelinka, M. D., and Hall, A.: An observational radiative constraint on hydrologic cycle intensification, Nature, 528, 249, 2015.


Example plots
-------------



.. _fig3b:
.. figure:: /recipes/figures/deangelis15nat/fig3b_cmug.png
   :align: center
   :width: 50%

   Scatter plot and regression line the between the ratio of the change of net short wave radiation (rsnst) and the change of the Water Vapor Path (prw) against the ratio of the change of netshort wave radiation for clear skye (rsnstcs) and the the change of surface temperature (tas). The width of horizontal shading for models and the vertical dashed lines for observations (Obs.) represent statistical uncertainties of the ratio, as the 95% confidence interval (CI) of the regression slope to the rsnst versus prw curve. For the observations the minimum of the lower bounds of all CIs to the maximum of the upper bounds of all CIs is shown.
