.. _recipes_wetland_ch4:

Wetland methane emissions
=====

Overview
--------

Wetland methane emissions in the Northern boreal region.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_wetland_ch4.yml

Diagnostics are stored in esmvaltool/diag_scripts/<mynewdiag>/

    * wetland_ch4/: one line scription
    * wetland_ch4/: one line scription
    * wetland_ch4/: one line scription


User settings in recipe
-----------------------

#. Script <mynewdiag.py/.ncl/.r>

   *Required settings for script*

   * xxx: zzz

   *Optional settings for script*

   *Required settings for variables*

   *Optional settings for variables*

   *Required settings for preprocessor*

   *Optional settings for preprocessor*

   *Color tables*

   * list required color tables (if any) here


Variables
---------

* wetlandCH4 (land, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
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