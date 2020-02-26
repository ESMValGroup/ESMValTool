Emergent constraints for equilibrium climate sensitivity
========================================================

Overview
--------

Calculates equilibrium climate sensitivity (ECS) versus

1) southern ITCZ index, similar to fig. 2 from Tian (2015)
2) lower tropospheric mixing index (LTMI), similar to fig. 5 from Sherwood et al. (2014)
3) tropical mid-tropospheric humidity asymmetry index, similar to fig. 4 from Tian (2015)
4) covariance of shortwave cloud reflection (Brient and Schneider, 2016)
5) climatological Hadley cell extent (Lipat et al., 2017)

The results are displayed as scatterplots.

.. note:: This recipe requires pre-calulation of the equilibrium climate
  sensitivites (ECS) for all models. The ECS values are calculated
  with recipe_ecs.yml. The netcdf file containing the ECS values
  (path and filename) is specified by diag_script_info@ecs_file.
  Alternatively, the netcdf file containing the ECS values can be
  generated with the cdl-script
  $diag_scripts/emergent_constraints/ecs_cmip.cdl (recommended method):
          
  1) save script given at the end of this recipe as ecs_cmip.cdl
  2) run command: ncgen -o ecs_cmip.nc ecs_cmip.cdl
  3) copy ecs_cmip.nc to directory given by diag_script_info@ecs_file
     (e.g. $diag_scripts/emergent_constraints/ecs_cmip.nc)


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_ecs_scatter.yml

Diagnostics are stored in diag_scripts/emergent_constraints/

    * ecs_scatter.ncl: calculate emergent constraints for ECS


User settings in recipe
-----------------------

#. Script ecs_scatter.ncl

   *Required settings (scripts)*

   * diag: emergent constraint to calculate ("itczidx", "humidx", "ltmi",
     "covrefl", "shhc")
   * ecs_file: path and filename of netCDF containing precalculated
     ECS values (see note above)

   *Optional settings (scripts)*

   * calcmm: calculate multi-model mean (True, False)
   * legend_outside: plot legend outside of scatterplots (True, False)
   * predef_minmax: use predefined internal min/max values for axes
     (True, False)
   * styleset: "CMIP5" (if not set, diagnostic will create a color table
     and symbols for plotting)
   * suffix: string to add to output filenames (e.g."cmip3")

   *Required settings (variables)*

   * reference_dataset: name of reference data set

   *Optional settings (variables)*

   none

   *Color tables*

   none


Variables
---------

* pr (atmos, monthly mean, longitude latitude time)
* hur (atmos, monthly mean, longitude latitude level time)
* hus (atmos, monthly mean, longitude latitude level time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsut (atmos, monthly mean, longitude latitude time)
* rsutcs (atmos, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude level time)
* ts (atmos, monthly mean, longitude latitude time)
* va (atmos, monthly mean, longitude latitude level time)
* wap (atmos, monthly mean, longitude latitude level time)


Observations and reformat scripts
---------------------------------

.. note:: (1) Obs4mips data can be used directly without any preprocessing.
          (2) See headers of reformat scripts for non-obs4mips data for download instructions.

* AIRS (obs4mips): hus, husStderr
* CERES-EBAF (obs4mips): rsdt, rsut, rsutcs
* ERA-Interim (OBS): hur, ta, va, wap
* HadISST (OBS): ts
* TRMM-L3 (obs4mips): pr, prStderr


References
----------

* Brient, F., and T. Schneider, J. Climate, 29, 5821-5835, doi:10.1175/JCLIM-D-15-0897.1, 2016.
* Lipat et al., Geophys. Res. Lett., 44, 5739-5748, doi:10.1002/2017GL73151, 2017.
* Sherwood et al., nature, 505, 37-42, doi:10.1038/nature12829, 2014.
* Tian, Geophys. Res. Lett., 42, 4133-4141, doi:10.1002/2015GL064119, 2015.

Example plots
-------------

.. _fig_ec_ecs_1:
.. figure::  /recipes/figures/emergent_constraints/ltmi.png
   :align:   center

   Lower tropospheric mixing index (LTMI; Sherwood et al., 2014) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_2:
.. figure::  /recipes/figures/emergent_constraints/shhc.png
   :align:   center

   Climatological Hadley cell extent (Lipat et al., 2017) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_3:
.. figure::  /recipes/figures/emergent_constraints/humidx.png
   :align:   center

   Tropical mid-tropospheric humidity asymmetry index (Tian, 2015) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_4:
.. figure::  /recipes/figures/emergent_constraints/itczidx.png
   :align:   center

   Southern ITCZ index (Tian, 2015) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_5:
.. figure::  /recipes/figures/emergent_constraints/covrefl.png
   :align:   center

   Covariance of shortwave cloud reflection (Brient and Schneider, 2016) vs.
   equilibrium climate sensitivity from CMIP5 models.
