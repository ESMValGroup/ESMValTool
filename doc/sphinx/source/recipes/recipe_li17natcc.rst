.. _recipes_li17natcc:

Constraining future Indian Summer Monsoon projections with the present-day precipitation over the tropical western Pacific
==========================================================================================================================

Overview
--------


Following `Li et al. (2017)`_ the change between present-day and future Indian Summer Monsoon (ISM) precipitation is constrained
using the precipitation over the tropical western Pacific compared to
a fixed, observed amount of 6 mm d\ :sup:`-1` from Global Precipitation Climatology Project (GPCP) `(Adler et al., 2003)`_ for 1980-2009.
For CMIP6, historical data for 1980-2009 should be used. For CMIP5 historical data from 1980-2005 should be used, due to the length of the data sets.
At the moment it is not possible to use a combined ``['historical', 'rcp']`` data set, because the diagnostic requires that a historical data set is given.

.. _`(Adler et al., 2003)`: https://journals.ametsoc.org/doi/abs/10.1175/1525-7541%282003%29004%3C1147%3ATVGPCP%3E2.0.CO%3B2
.. _`Li et al. (2017)`: https://www.nature.com/articles/nclimate3387


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_li17natcc.yml


Diagnostics are stored in diag_scripts/

   * emergent_constraints/lif1f2.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models. For each model, two experiments must be given: 
one historical run, possibly between 1980-2009 and one other model experiment. The user can choose the other model experiment, 
but it needs to be the same for all given models. 
The start and end year for the second data set can be choosen by the user, but should be consistent for all models 
(the same for future scenarios, the same length for other experiments). Different ensemble members are not possible, yet.


Variables
---------

* *pr* (atmos, monthly, longitude, latitude, time)
* *ua* (atmos, monthly, longitude, latitude, plev, time)
* *va* (atmos, monthly, longitude, latitude, plev, time)
* *ts* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Li, G., Xie, S. P., He, C., and Chen, Z. S.: Western Pacific emergent constraint lowers projected increase in Indian summer monsoon rainfall, Nat Clim Change, 7, 708-+, 2017


Example plots
-------------

.. _li17natcc_fig2a:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2a.png
   :align: center
   :width: 50%

   Scatter plot of the simulated tropical western Pacific precipitation (mm d\ :sup:`-1`\ ) versus projected average ISM (Indian Summer Monsoon) rainfall changes under the ssp585 scenario. The red line denotes the observed present-day western Pacific precipitation and the inter-model correlation (r) is shown. (CMIP6).

.. _li17natcc_fig2b:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2b.png
   :align: center
   :width: 50%

   Scatter plot of the uncorrected versus corrected average ISM (Indian Summer Monsoon) rainfall change ratios (% per degree Celsius of global SST warming). The error bars for the Multi-model mean indicate the standard deviation spread among models and the 2:1 line (y = 0.5x) is used to illustrate the Multi-model mean reduction in projected rainfall increase. (CMIP6).

.. _li17natcc_fig2c:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2c.png
   :align: center
   :width: 50%

   Multi-model mean rainfall change due to model error. Box displays the area used to define the average ISM (Indian Summer Monsoon) rainfall. Precipitation changes are normalized by the corresponding global mean SST increase for each model. (CMIP6).

.. _li17natcc_fig2d:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2d.png
   :align: center
   :width: 50%

   Corrected multi-model mean rainfall change. Box displays the area used to define the average ISM (Indian Summer Monsoon) rainfall. Precipitation changes are normalized by the corresponding global mean SST increase for each model. (CMIP6).
