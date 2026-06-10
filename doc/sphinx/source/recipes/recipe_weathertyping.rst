.. _recipes_weathertyping:

Lamb Weathertypes
===================

Overview
--------

A diagnostic to calculate Lamb weathertypes over a given region. Furthermore,
correlations between weathertypes and precipitation patterns over a given area can be calculated
and 'combined' or 'simplified' weathertypes can be derived. Additionally, mean fields, as well as
anomalies and standard deviations can be plotted.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_weathertyping.yml

Diagnostics are stored in esmvaltool/diag_scripts/weathertyping/

    * weathertyping.py: calculate lamb and simplified WT, plot mean, anomalies and std for each WT for psl, tas, pr


User settings in recipe
-----------------------

#. weathertyping.py

   *Required settings for script*

   *Optional settings for script*

   * correlation_threshold: correlation threshold for selecting similar WT pairs, only needed if automatic_slwt==True and predefined_slwt==False. default=0.9
   * rmse_threshold: rmse threshold for selecting similar WT pairs, only needed if automatic_slwt==True and predefined_slwt==False. default=0.002
   * plotting: if true, create plots of means, anomalies and std for psl, tas, prcp
   * automatic_slwt: if true, automatically combine WT with similar precipitation patterns over specified area (via thresholds of correlation and rmse OR via predefined_slwt)
   * predefined_slwt: dictionary of mappings between weathertypes

.. note::

  predefined_slwt can be a dictionary where keys are slwt and the values are arrays of lwt OR where keys are lwt and values are slwt

   *Required settings for variables*

   *Optional settings for variables*

   *Required settings for preprocessor*

   *Optional settings for preprocessor*

   *Color tables*


Variables
---------

* psl (atmos, day, time longitude latitude)
* tas (atmos, day, time longitude latitude)
* tp (atmos, day, time longitude latitude)
* pr (atmos, day, time longitude latitude)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

This recipe currently only works with the following reanalysis and observation datasets:

* E-OBS: European Climate Assessment & Dataset gridded daily precipitation sum
* ERA5: ECMWF reanalysis

References
----------

* Maraun, D., Truhetz, H., & Schaffer, A. (2021). Regional climate model biases, their dependence on synoptic circulation biases and the potential for bias adjustment: A process-oriented evaluation of the Austrian regional climate projections. Journal of Geophysical Research: Atmospheres, 126, e2020JD032824. https://doi.org/10.1029/2020JD032824

* Jones, P.D., Hulme, M. and Briffa, K.R. (1993), A comparison of Lamb circulation types with an objective classification scheme. Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

Example plots
-------------

.. _fig_weathertyping_1:
.. figure::  /recipes/figures/weathertyping/lwt_1_ERA5__psl_mean_1958-2014.png
   :align:   center

   PSL mean map of Lamb WT 1 for ERA5.

.. _fig_weathertyping_2:
.. figure::  /recipes/figures/weathertyping/lwt_1_TaiESM1_r1i1p1f1_psl_mean_1950-2014.png
   :align:   center

   PSL mean map of Lamb WT 1 for TaiESM1.

.. _fig_weathertyping_3:
.. figure::  /recipes/figures/weathertyping/slwt_EOBS_4_ERA5__psl_mean_1958-2014.png
   :align:   center

   PSL mean map of slwt_EOBS 4 for ERA5 (in this case combined Lamb WT 24 and 23).

.. _fig_weathertyping_4:
.. figure::  /recipes/figures/weathertyping/correlation_matrix_E-OBS_1958-2014.png
   :align:   center

   Heatmap of correlation values for Lamb WTs 1-27.

.. _fig_weathertyping_5:
.. figure::  /recipes/figures/weathertyping/ERA5__lwt_rel_occurrence_1958-2014.png
   :align:   center

   Stackplot of seasonal relative occurrences of each WT for ERA5.
