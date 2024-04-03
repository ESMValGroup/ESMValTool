.. _recipe_perfmetrics_python:

Performance metrics for essential climate parameters in python
==============================================================

.. note::

      This recipe uses python diagnostics to reproduce parts of the evaluation
      done in the
      :ref:`original recipe based on NCL diagnostics <nml_perfmetrics>`:
      It aims for a complete replacement of all involved NCL diagnostics. So
      far, only portrait plots (including performance metrics) are supported.

Overview
--------

The goal is to create a standard recipe for the calculation of performance metrics to quantify the ability of the models to reproduce the climatological mean annual cycle for selected "Essential Climate Variables" (ECVs) plus some additional corresponding diagnostics and plots to better understand and interpret the results.

The recipe can be used to calculate performance metrics at different vertical levels (e.g., 5, 30, 200, 850 hPa as in `Gleckler et al. (2008) <http://dx.doi.org/10.1029/2007JD008972>`_) and in different regions. As an additional reference, we consider `Righi et al. (2015) <https://doi.org/10.5194/gmd-8-733-2015>`_.
Brief description of the diagnostic.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_perfmetrics_python.yml
    * recipe_perfmetrics_CMIP5_python.yml

Diagnostics are stored in esmvaltool/diag_scripts/perfmetrics/

    * portrait_plot.py: Plot metrics for any variable for multiple datasets and
      up to four references.


User settings in recipe
-----------------------

#. Script perfmetrics/portrait_plot.py

   This plot expects a scalar value in each input file and at most one input
   file for each subset of metadata that belongs to a cell or part of cell in
   the figure.
   By default cells are plotted for combinations of `short_name`,
   `dataset`, `project` and `split`.
   Where `split` is an optional extra_facet for variables.
   However, all this can be customized using the `x_by`,
   `y_by`, `group_by` and `split_by` script settings.
   For a complete and detailed list of settings see the
   :ref:`API documentation <api.esmvaltool.diag_scripts.perfmetrics.portrait_plot>`.
   While this allows very flexible use for any kind of data, there are some
   limitations as well: The grouping (separated
   plots in the figure) and normalization is always applied along the x-axis.
   With default settings this means normalizing all metrics for each variable
   and grouping all datasets by project.

   To plot distance metrics like RMSE, pearson R, bias etc. the
   :func:`distance_metrics <esmvalcore.preprocessor.derive>` preprocessor or
   custom diagnostics can be used.



Variables
---------

.. note::

   The recipe generally works for any variable that is preprocessed correctly.
   To use different preprocessors or reference datasets it could be useful
   to create different variable groups and link them with the same extra_facet
   like `variable_name` See recipe for examples. Below listed are the variables
   needed to produce the example figures.


#.  recipe_perfmetrics_CMIP5.yml

   * clt (atmos, monthly mean, longitude latitude time)
   * hus (atmos, monthly mean, longitude latitude lev time)
   * od550aer, od870aer, od550abs, od550lt1aer (aero, monthly mean, longitude latitude time)
   * pr (atmos, monthly mean, longitude latitude time)
   * rlut, rlutcs, rsut, rsutcs (atmos, monthly mean, longitude latitude time)
   * sm (land, monthly mean, longitude latitude time)
   * ta (atmos, monthly mean, longitude latitude lev time)
   * tas (atmos, monthly mean, longitude latitude time)
   * toz (atmos, monthly mean, longitude latitude time)
   * ts (atmos, monthly mean, longitude latitude time)
   * ua (atmos, monthly mean, longitude latitude lev time)
   * va (atmos, monthly mean, longitude latitude lev time)
   * zg (atmos, monthly mean, longitude latitude lev time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

The following list shows the currently used observational data sets for this recipe with their variable names and the reference to their respective reformat scripts in parentheses. Please note that obs4MIPs data can be used directly without any reformatting. For non-obs4MIPs data use `esmvaltool data info DATASET` or see headers of cmorization scripts (in `/esmvaltool/cmorizers/data/formatters/datasets/
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/data/formatters/datasets/>`_) for downloading and processing instructions.

#.  recipe_perfmetrics_CMIP5.yml

    * AIRS (hus - obs4MIPs)
    * CERES-EBAF (rlut, rlutcs, rsut, rsutcs - obs4MIPs)
    * ERA-Interim (tas, ta, ua, va, zg, hus - esmvaltool/cmorizers/data/formatters/datasets/era-interim.py)
    * ESACCI-AEROSOL (od550aer, od870aer, od550abs, od550lt1aer - esmvaltool/cmorizers/data/formatters/datasets/esacci-aerosol.ncl)
    * ESACCI-CLOUD (clt - esmvaltool/cmorizers/data/formatters/datasets/esacci-cloud.ncl)
    * ESACCI-OZONE (toz - esmvaltool/cmorizers/data/formatters/datasets/esacci-ozone.ncl)
    * ESACCI-SOILMOISTURE (sm - esmvaltool/cmorizers/data/formatters/datasets/esacci_soilmoisture.ncl)
    * ESACCI-SST (ts - esmvaltool/ucmorizers/data/formatters/datasets/esacci-sst.py)
    * GPCP-SG (pr - obs4MIPs)
    * HadISST (ts - esmvaltool/cmorizers/data/formatters/datasets/hadisst.ncl)
    * MODIS (od550aer - esmvaltool/cmorizers/data/formatters/datasets/modis.ncl)
    * NCEP-NCAR-R1 (tas, ta, ua, va, zg - esmvaltool/cmorizers/data/formatters/datasets/ncep_ncar_r1.py)
    * NIWA-BS (toz - esmvaltool/cmorizers/data/formatters/datasets/niwa_bs.ncl)
    * PATMOS-x (clt - esmvaltool/cmorizers/data/formatters/datasets/patmos_x.ncl)


References
----------


* Gleckler, P. J., K. E. Taylor, and C. Doutriaux, Performance metrics for climate models, J.
  Geophys. Res., 113, D06104, doi: 10.1029/2007JD008972 (2008).

* Righi, M., Eyring, V., Klinger, C., Frank, F., Gottschaldt, K.-D., Jöckel, P.,
  and Cionni, I.: Quantitative evaluation of oone and selected climate parameters in a set of EMAC simulations,
  Geosci. Model Dev., 8, 733, doi: 10.5194/gmd-8-733-2015 (2015).


Example plots
-------------

.. _fig_perfmetrics_python_portrait_plot:

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_5_python.png
   :width: 90%
   :align: center


   Relative space-time root-mean-square deviation (RMSD) calculated from the climatological
   seasonal cycle of CMIP5 simulations. A relative performance is displayed, with blue shading
   indicating better and red shading indicating worse performance than the median of all model results.
   A diagonal split of a grid square shows the relative error with respect to the reference data set
