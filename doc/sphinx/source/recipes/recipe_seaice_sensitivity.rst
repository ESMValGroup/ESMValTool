.. _recipe_seaice_sensitivity:

Sea Ice Sensitivity
===================

Overview
--------

This recipe and diagnostic calculates the rate of sea ice area loss per degree of warming, as in plot 1d of `Notz et al.`_ and figure 3e of `Roach et al.`_.

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
      <<: *extract_test_period
      <<: *extract_sept
      <<: *nh_total_area

    pp_antarctic_avg_ann_sea_ice:
      <<: *extract_test_period
      <<: *annual_mean
      <<: *sh_total_area

    pp_avg_ann_global_temp:
      <<: *extract_test_period
      <<: *global_mean
      <<: *annual_mean

Datasets
--------

The recipe tries to use as many datasets as possible from the list in Table S3 of https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2019GL086749&file=grl60504-sup-0002-Table_SI-S01.pdf

References
----------

* SIMIP_2020: Scripts for the processing, analysis and plots of the SIMIP Community Paper (Notz et al. 2020), https://github.com/jakobdoerr/SIMIP_2020/tree/master

Example plots
-------------

.. _fig_seaice_sensitivity_1:
.. figure::  /recipes/figures/seaice/September_Arctic_sea_ice_sensitivity.png
   :align:   center
   :width:   8cm

   Plot of northern hemisphere sea ice area loss (millions of square kilometres) in the month of September per degree Kelvin of average global warming.

   The dashed black line shows the observational mean, the shaded area is within one standard deviation of the mean, and the dashed grey lines are within "plausibility" of the mean, with all values taken from Ed Blockley's code for the period 1979-2014 and defined as follows:

   * mean:      -4.01,
   * std_dev:  0.32,
   * plausible: 1.28,

.. _fig_seaice_sensitivity_2:
.. figure::  /recipes/figures/seaice/Annual_Antarctic_sea_ice_trends.png
   :align:   center
   :width:   18cm

   Plot of the trend of annually averaged southern hemisphere sea ice area (millions of square kilometres) over time against the trend of annually and globally averaged air temperature near the surface (degrees Kelvin) over time.

   The colour of each point is determined by the ``r value`` of the correlation between the two variables, and the hatching indicates a ``p value`` greater than ``0.05``, both calculated using ``scipy.stats.linregress``.
