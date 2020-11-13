.. _recipes_consecdrydays:

Consecutive dry days
====================

Overview
--------
Meteorological drought can in its simplest form be described by a lack of precipitation. First, a wet day threshold is set, which can be either a limit related to measurement accuracy, or more directly process related to an amount that would break the drought. The diagnostic calculates the longest period of consecutive dry days, which is an indicator of the worst drought in the time series. Further, the diagnostic calculates the frequency of dry periods longer than a user defined number of days.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_consecdrydays.yml

Diagnostics are stored in diag_scripts/droughtindex/

    * diag_cdd.py: calculates the longest period of consecutive dry days, and
      the frequency of dry day periods longer than a user defined length


User settings in recipe
-----------------------

#. Script diag_cdd.py

   *Required settings (script)*

   * plim: limit for a day to be considered dry [mm/day]

   * frlim: the shortest number of consecutive dry days for entering statistic on frequency of dry periods.


Variables
---------

* pr      (atmos, daily mean, time latitude longitude)


Example plots
-------------

.. _fig_consecdrydays:
.. figure::  /recipes/figures/consecdrydays/consec_example_freq.png
   :align:   center
   :width:   14cm

   Example of the number of occurrences with consecutive dry days of more than five days in the period 2001 to 2002 for the CMIP5 model bcc-csm1-1-m. 
