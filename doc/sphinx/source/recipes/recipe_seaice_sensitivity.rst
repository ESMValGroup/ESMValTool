.. _recipe_seaice_sensitivity:

Sea Ice Sensitivity
===================

Overview
--------

This recipe and diagnostic calculates the rate of sea ice area loss per degree of warming, as in plot 1d of `Notz et al.`_ and figure 3e of `Roach et al.`_. The figures used in the plots are output to csv files in the work directory.

.. _`Notz et al.`: https://doi.org/10.1029/2019GL086749
.. _`Roach et al.`: https://doi.org/10.1029/2019GL086729

Available recipe and diagnostic
-------------------------------

Recipe is in `recipes`

* recipe_seaice_sensitivity.yml

Diagnostic is in `diag_scripts/seaice/`

* seaice_sensitivity.py (Plotting the sensitivity of sea ice to mean global temperature).

Recipe settings
~~~~~~~~~~~~~~~

Years to be evaluated are specified in the ``extract_test_period`` preprocessor. Arctic sea ice is evaluated using data from September, and Antarctic sea ice is evaluated using annually meaned data. This can be amended by changing the preprocessor for each variable (shown below).

.. code-block:: yaml

    pp_arctic_sept_sea_ice:
      extract_time:
        start_day: 1
        start_month: 1
        start_year: 1979
        end_day: 31
        end_month: 12
        end_year: 2014
      extract_month:
        month: 9
      extract_region:
        start_longitude: 0
        end_longitude: 360
        start_latitude: 0
        end_latitude: 90
      area_statistics:
        operator: sum
      convert_units:
        units: 1e6 km2

Datasets
--------

The recipe tries to use as many datasets as possible from the list in Table S3 of https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2019GL086749&file=grl60504-sup-0002-Table_SI-S01.pdf
Only one ensemble member is used for each model.

All, some or no datasets may be labelled in the plots, using ``label_dataset: True`` in the recipe settings.

References
----------

* SIMIP_2020: Scripts for the processing, analysis and plots of the SIMIP Community Paper (Notz et al. 2020), https://github.com/jakobdoerr/SIMIP_2020/tree/master

Example plots
-------------

.. _fig_seaice_sensitivity_1:
.. figure::  /recipes/figures/seaice/September_Arctic_sea_ice_sensitivity.png
   :align:   center
   :width:   8cm

   Plot of sensitivity of northern hemisphere sea ice area loss (millions of square kilometres) in the month of September to the annual mean global temperature change (K).

   The dashed black line shows the observational mean, the shaded area denotes one one standard deviation of observational uncertainty, as calculated by Notz et al (2020).
   The dotted grey lines reflect Notz et al estimate of a plausible range incorporating both internal variability and observational uncertainty.
   These values are configurable in the recipe, with the default values taken from Notz et al (2020):

=========  =====================
mean       -4.01 million km2 K-1
std_dev    0.32 million km2 K-1
plausible  1.28 million km2 K-1
=========  =====================

.. _fig_seaice_sensitivity_2:
.. figure::  /recipes/figures/seaice/Annual_Antarctic_sea_ice_trends.png
   :align:   center
   :width:   18cm

   Plot of the trend of annually averaged southern hemisphere sea ice area (millions of square kilometres) over time against the trend of annually and globally averaged air temperature near the surface (degrees Kelvin) over time. The values plotted are 10 times the annual trend, which was calculated using :func:`scipy.stats.linregress`, for consistency with the decadal values used in the published plot.

   The colour of each point is determined by the Pearson correlation coefficient between the two variables, and the hatching indicates a ``p_value`` greater than 0.05, both calculated using :func:`scipy.stats.linregress`.
