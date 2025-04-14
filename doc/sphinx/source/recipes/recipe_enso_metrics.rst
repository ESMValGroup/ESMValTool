.. _recipes_enso_metrics:

ENSO CLIVAR metrics - reproducing ENSO characteristics lifecycle, seasonality and teleconnections
=================================================================================================

Overview
--------

Reproducing some ENSO metrics from the 
[CLIVAR package](https://pcmdi.llnl.gov/pmp-preliminary-results/interactive_plot/portrait_plot/enso_metric/enso_metrics_interactive_portrait_plots_v20231121.html)
specifically ones used for REF

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_enso_metrics.yml

Diagnostics are stored in esmvaltool/diag_scripts/enso_metrics/

* **enso_diag1metrics.py**: metrics for basic ENSO characteristics
* **matrix.py**: reads metrics in work_dir from csv file written out in climatology_diagnostic1, use for other groups of metrics


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

.. _fig_seasonality:
.. figure:: /recipes/figures/enso_metrics/seasonality.png
   :align: center

   Ratio of winter over spring's standard deviation of sea surface temperature anomalies (SSTA) in the central equatorial Pacific (Niño3.4 averaged), showing the seasonal timing of SSTA. 
   All models are on the same line.

.. _fig_lifecycle:
.. figure:: /recipes/figures/enso_metrics/lifecycle.png
   :align: center

   Temporal structure of sea surface temperature anomalies (SSTA) in the central equatorial Pacific (Niño3.4 averaged), showing the temporal evolution of SSTA associated with ENSO. 
   Observation is dashed black line.

.. _fig_teleconnections:
.. figure:: /recipes/figures/enso_metrics/ACCESS-CM2_DJF_ts_telecon.png
   :align: center

   PR or SST anomalies on Earth (between 60°S-60°N), showing the location associated with ENSO.

.. _fig_metrics:
.. figure:: /recipes/figures/enso_metrics/plot_matrix.png
   :align: center

   Portrait plot of normalized metrics for models against reference observation.