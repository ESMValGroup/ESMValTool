.. _recipe_autoassess_landsurface_permafrost.rst:

Land-surface Permafrost - Autoassess diagnostics
================================================

Overview
--------


Prior and current contributors
------------------------------
Met Office:

* Prior to April 2018: Eleanor Burke, Paul Earnshaw

ESMValTool:

* Since April 2018: Porting into ESMValTool by Alistair Sellar and Valeriu Predoi


Developers
----------
Met Office:

* Eleanor Burke


ESMValTool:

* Since April 2018: Alistair Sellar and Valeriu Predoi

Review of current port in ESMValTool
------------------------------------
The code and results review of the port from native Autoassess to ESMValTool
was conducted by Alistair Sellar (`<alistair.sellar@matoffice.gov.uk>`_) and
Valeriu Predoi (`<valeriu.predoi@ncas.ac.uk>`_) in May 2021. Review consisted in
comparing results from runs using ESMValTool's port and native Autoassess using
the same models and data stretches.

Metrics and Diagnostics
-----------------------

Performance metrics (with observation-based estimates in brackets):

* permafrost area (17.46 million square km)
* fractional area of permafrost northwards of zero degree isotherm (0.47)
* soil temperature at 1m minus soil temperature at surface (-0.53 degrees C)
* soil temperature at surface minus air temperature (6.15 degrees C)
* annual amplitude at 1m / annual amplitude at the surface (0.40 unitless)
* annual amplitude at the surface / annual air temperature (0.57 unitless)


Diagnostics:

* Maps of permafrost extent and zero degC isotherm


Model Data
----------

========================================================================= ================== ============== ==============================================
Variable/Field name                                                       realm              frequency      Comment
========================================================================= ================== ============== ==============================================
Near-Surface Air Temperature (tas)                                        Atmosphere         monthly mean
Temperature of Soil (tsl)                                                 Land               monthly mean
Total Water Content of Soil Layer (mrsol)                                 Land               monthly mean   CMIP5: mrlsl
Percentage of the Grid Cell Occupied by Land (Including Lakes) (sftlf)    mask               fixed
========================================================================= ================== ============== ==============================================

The recipe takes as input a control model and experimental model, comparisons being made
with these two CMIP models.

Inputs and usage
----------------
The ``landsurface_permafrost`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
and, as any other ``autoassess`` metric, it uses the ``autoassess_area_base.py`` as general purpose
wrapper. This wrapper accepts a number of input arguments that are read through from the recipe.

This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other Autoassess metrics.

.. note::

   **Time gating for autoassess metrics.**

   To preserve the native Autoassess functionalities,
   data loading and selection on time is done somewhat
   differently for ESMValTool's autoassess metrics: the
   time selection is done in the preprocessor as per usual but
   a further time selection is performed as part of the diagnostic.
   For this purpose the user will specify a ``start:`` and ``end:``
   pair of arguments of ``scripts: autoassess_script`` (see below
   for example). These are formatted as ``YYYY/MM/DD``; this is
   necessary since the Autoassess metrics are computed from 1-Dec
   through 1-Dec rather than 1-Jan through 1-Jan. This is a temporary
   implementation to fully replicate the native Autoassess functionality
   and a minor user inconvenience since they need to set an extra set of
   ``start`` and ``end`` arguments in the diagnostic; this will be phased
   when all the native Autoassess metrics hanve been ported to ESMValTool
   review has completed.


An example of standard inputs as read by ``autoassess_area_base.py`` and passed
over to the diagnostic/metric is listed below.


.. code-block:: yaml

    scripts:
      plot_landsurf_permafrost: &plot_landsurf_permafrost_settings
        <<: *autoassess_landsurf_permafrost_settings
        control_model: MPI-ESM-LR
        exp_model: MPI-ESM-MR
        script: autoassess/plot_autoassess_metrics.py
        ancestors: ['*/autoassess_landsurf_permafrost']
        title: "Plot Land-Surface Permafrost Metrics"
        plot_name: "Permafrost_Metrics"
        diag_tag: aa_landsurf_permafrost
        diag_name: autoassess_landsurf_permafrost

References
----------
* Observed permafrost extent is from http://nsidc.org/data/ggd318.html: Brown, J.,
O. Ferrians, J. A. Heginbottom, and E. Melnikov. 2002. Circum-Arctic Map of
Permafrost and Ground-Ice Conditions, Version 2. Boulder, Colorado USA. NSIDC:
National Snow and Ice Data Center.  When calculating the global area of
permafrost the grid cells are weighted by the proportion of permafrost within
them.

* Annual mean air temperature is from: Legates, D. R., and C. J. Willmott, 1990:
Mean seasonal and spatial variability in global surface air temperature. Theor.
Appl. Climatol., 41, 11-21.  The annual mean is calculated from the seasonal
mean data available at the Met Office.

* The soil temperature metrics are calcuated following: Charles D. Koven, William
J. Riley, and Alex Stern, 2013: Analysis of Permafrost Thermal Dynamics and
Response to Climate Change in the CMIP5 Earth System Models. J. Climate, 26. 
(Table 3) http://dx.doi.org/10.1175/JCLI-D-12-00228.1 The
locations used for Table 3 were extracted from the model and the modelled
metrics calculated.


Observations Data sets
----------------------

None used in this diagnostic.

Sample Plots and metrics
------------------------
Below is a set of metrics for  UKESM1-0-LL (historical data); the table
shows a comparison made between running ESMValTool on CMIP6 CMORized
netCDF data freely available on ESGF nodes and the run made using native
Autoassess performed at the Met Office using the pp output of the model.
Comparison period was 1997/12/01 to 2002/12/01.

===============================================     ================     ====================
Metric name                                         UKESM1-0-LL;         UKESM1-0-LL;
                                                    CMIP6: AERmonZ;      pp files;
                                                    piControl, ESGF      piControl, u-aw310
===============================================     ================     ====================
attenuation 1m over surface                         0.496                0.496
attenuation surface over air                        0.492                0.493
fraction area permafrost over zerodeg               0.290                0.289
offset 1m minus surface                             0.947                0.947
offset surface minus air                            7.67                 7.71
permafrost area                                     13.5                 13.7
===============================================     ================     ====================


.. figure:: /recipes/figures/autoassess_landsurface/pf_extent_north_america_UKESM1-0-LL.png
   :scale: 50 %
   :alt: pf_extent_north_america_UKESM1-0-LL.png

   Permafrost extent and zero degC isotherm, showing North America

.. figure:: /recipes/figures/autoassess_landsurface/pf_extent_asia_UKESM1-0-LL.png
   :scale: 50 %
   :alt: pf_extent_asia_UKESM1-0-LL.png

   Permafrost extent and zero degC isotherm, showing Asia and Europe

