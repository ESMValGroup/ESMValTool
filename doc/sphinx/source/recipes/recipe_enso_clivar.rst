.. _recipes_enso_diagnostics:

ENSO CLIVAR metrics - reproducing ENSO characteristics, lifecycle, seasonality, teleconnections
=================================================================================================

Overview
--------

Reproducing ENSO metrics from the
[CLIVAR package](https://pcmdi.llnl.gov/pmp-preliminary-results/interactive_plot/portrait_plot/enso_metric/enso_metrics_interactive_portrait_plots_v20231121.html)
with dive down plots of ENSO diagnostics including characteristics, lifecycle, seasonality, feedbacks and teleconnections.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_enso_climatology1.yml
* recipe_enso_climatology_level2-3.yml
* recipe_enso_characteristics.yml
* recipe_enso_divedown9-10_14-15.yml
* recipe_enso_divedown11-13.yml
* recipe_enso_feedback.yml
* recipe_enso_feedbacknhf.yml
* recipe_enso_teleconnection.yml

Diagnostics are stored in esmvaltool/diag_scripts/enso_metrics/

* **matrix.py**: reads metrics in work_dir from csv file written out in climatology_diagnostic1, use for other groups of metrics
* **climatology_diagnostic1.py**: computes basic climatology for level 1 metrics
* **climatology_diagnosticlevel2.py**: creates basic climatology level 2 plots
* **climatology_diagnosticlevel3.py**: creates basic climatology level 3 plots
* **enso_diag1metrics.py**: computes basic ENSO characteristics metrics
* **enso_divedowns/ampseasskew11-13.py**: dive down scripts for ENSO amplitude, seasonality and skewness
* **enso_divedowns/durdiver_14-15.py**: dive down scripts for ENSO duration and diversity
* **enso_divedowns/lifecycle_10.py**: dive down scripts for ENSO lifecycle
* **enso_divedowns/pattern_9.py**: dive down scripts for ENSO pattern
* **feedback/feedback_metrics_lvl4.py**: feedback level4 plots
* **feedback/feedback_metrics.py**: feedback level3 plots
* **feedback/feedback_nhf_lvl3_4.py**: feedback NHF level3-4 plots
* **feedback/feedback_nhf_metric.py**: feedback NHF metrics
* **teleconnections_metrics.py**: computes and plots teleconnections metrics


User settings in recipe
-----------------------

#. Script: **matrix.py**

   *Required settings for script*

   * **diag_metrics**: diagnostic name and script name in *yml* of the diagnostic that computes all the metrics so it can find the *csv* in the `work_dir` - eg. diagnostic_metrics/plot_script


Variables
---------

* tos (Omon, monthly)
* areacello (Ofx)
* pr (Amon, monthly)
* ts (Amon, monthly)
* areacella (fx)
* tauu (Amon, monthly)
* zos (Omon, monthly)
* hfns (hfls+hfss, Amon, monthly)

Observations and reformat scripts
---------------------------------

* HadISST
* TropFLUX
* GPCP-SG
* ERA-Interim

References
----------

* [1] https://pcmdi.llnl.gov/pmp-preliminary-results/interactive_plot/portrait_plot/enso_metric/enso_metrics_interactive_portrait_plots_v20231121.html
* [2] https://github.com/CLIVAR-PRP/ENSO_metrics/

Example plots
-------------

.. _fig_teleconnections:
.. figure:: /recipes/figures/enso_metrics/ACCESS-ESM1-5_DJF_ts_telecon.png
   :align: center

   PR or SST anomalies on Earth (between 60°S-60°N), showing the location associated with ENSO.

.. _fig_seasonality_level4:
.. figure:: /recipes/figures/enso_metrics/ACCESS-ESM1-5_12seasonality_level_4.png
   :align: center

   Seasonality level 4 plot of boreal summer/winter SSTA standard deviation.

.. _fig_feedback_level2:
.. figure:: /recipes/figures/enso_metrics/ACCESS-ESM1-5_SSH_SST_lvl2.png
   :align: center

   Feedback level 2 plot of SSH and SST anomalies.

.. _fig_metrics:
.. figure:: /recipes/figures/enso_metrics/plot_matrix.png
   :align: center

   Portrait plot of normalized metrics for models against reference observation.
