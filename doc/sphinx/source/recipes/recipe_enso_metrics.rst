.. _recipes_enso_metrics:

ENSO CLIVAR metrics - reproducing ENSO characteristics lifecycle and seasonality
================================================================================

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

   * **diag_metrics**: diagnostic name and script name in *yml* of the diagnostic that computes all the metrics 
   so it can find the *csv* in the `work_dir` - eg. diagnostic_metrics/plot_script


Variables
---------

* tos (Omon, monthly)
* areacello (Ofx)



Observations and reformat scripts
---------------------------------


* HadISST
* TropFLUX


References
----------

* [1] https://pcmdi.llnl.gov/pmp-preliminary-results/interactive_plot/portrait_plot/enso_metric/enso_metrics_interactive_portrait_plots_v20231121.html
* [2] https://github.com/CLIVAR-PRP/ENSO_metrics/

Example plots
-------------

.. _fig_seasonality:
.. figure::  /recipes/figures/enso_metrics/seasonality.png
   :align:   center

   Add figure caption here.

.. _fig_lifecycle:
.. figure::  /recipes/figures/enso_metrics/lifecycle.png
   :align:   center

   Add figure caption here.
