.. _recipes_<mynewrecipe>:

ESA CCI LST comparison to Historical Models
=====

Overview
--------

This diagnostic compares ESA CCI LST to multiple historical emsemble members of CMIP models.
It does this over a defined region for monthly values of the land surface temperature.
The result is a plot showing the mean differnce of CCI LST to model average LST, with a region on +/- one standard deviation of the model mean LST given as a measure of model variability.

The recipe and diagnostic need the all time average monthly LST from the CCI data.
We use the L3C single sensor monthy data.
A CMORizing script <path here> calculates the mean of the day time, and night time overpasses to give the all time average LST.
This is so that the Amon output from CMIP models can be used.
We created such a dataset from the Aqua MODIS data from CCI.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_esacci_lst.yml

Diagnostics are stored in esmvaltool/diag_scripts/lst/

    * <lst.py>: ***


User settings in recipe
-----------------------

#. Script <recipe_esacci_lst.yml>

   *Required settings for script*

   * No user defined inputs to the diagnostic

   *Required settings for variables*
   *** setting from recipe for model type....


   *Required settings for preprocessor*
   **** region
   

Variables
---------

* ts (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

*** CMOR stuff here ***

*Note: (1) obs4mips data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4mips data for download
instructions.*

* xxx

  *Reformat script:* <myreformatscript.py>

References
----------

* xxx

Example plots
-------------

.. _fig_mynewdiag_1:
.. figure::  /recipes/figures/<mynewdiagnostic>/awesome1.png
   :align:   center

   Add figure caption here.
