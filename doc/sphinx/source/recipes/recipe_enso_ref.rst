.. _recipes_enso_metrics:

ENSO CLIVAR metrics - reproducing background climatology and ENSO characteristics
=================================================================================

Overview
--------

Reproducing some ENSO metrics from the 
[CLIVAR package](https://pcmdi.llnl.gov/pmp-preliminary-results/interactive_plot/portrait_plot/enso_metric/enso_metrics_interactive_portrait_plots_v20231121.html)


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_enso_characteristics.yml
* recipe_enso_basicclimatology.yml

Diagnostics are stored in esmvaltool/diag_scripts/enso_metrics/

* **enso_diag1metrics.py**: metrics for basic ENSO characteristics.
* **climatology_diagnostic1.py**: metrics for background climatology.
* **climatology_diagnosticlevel2.py**: creating map plots for background climatology.


User settings in recipe
-----------------------



Variables
---------

* tos (Omon, monthly)
* areacello (Ofx)
* pr (Amon, monthly)
* ts (Amon, monthly)
* tauu (Amon, monthly)


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
.. figure:: /recipes/figures/enso_metrics/ACCESS-CM2_12seasonality.png
   :align: center

   Ratio of winter over spring's standard deviation of sea surface temperature anomalies (SSTA) in the central equatorial Pacific (Niño3.4 averaged), showing the seasonal timing of SSTA. 
   All models are on the same line.

.. _fig_lifecycle:
.. figure:: /recipes/figures/enso_metrics/ACCESS-CM2_10lifecycle.png
   :align: center

   Temporal structure of sea surface temperature anomalies (SSTA) in the central equatorial Pacific (Niño3.4 averaged), showing the temporal evolution of SSTA associated with ENSO. 
   Observation is dashed black line.
