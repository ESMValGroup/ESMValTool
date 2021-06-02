.. _recipe_autoassess_landsurface_soilmoisture.rst:

Land-surface Soil Moisture - Autoassess diagnostics
===================================================

Overview
--------


Prior and current contributors
------------------------------
Met Office:

* Prior to April 2018: Heather Rumbold, Paul Earnshaw

ESMValTool:

* Since April 2018: Porting into ESMValTool by Alistair Sellar and Valeriu Predoi


Developers
----------
Met Office:

* Heather Rumbold


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

Performance metrics:

* median absolute error (model minus observations)

Metrics are calculated using model and observation multi-year climatologies (seasonal means) 
for meteorological seasons:
* December-January-February (djf)
* March-April-May (mam)
* June-July-August (jja)
* September-October-November (son)

Diagnostics:

* Metrics plot comparing control and experiment



Model Data
----------

=========================================== ================== ============== ==============================================
Variable/Field name                         realm              frequency      Comment
=========================================== ================== ============== ==============================================
Total Water Content of Soil Layer (mrsos)   Land               monthly mean
=========================================== ================== ============== ==============================================

The recipe takes as input a control model and experimental model, comparisons being made
with these two CMIP models.

Inputs and usage
----------------
The ``landsurface_soilmoisture`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
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
      autoassess_landsurf_soilmoisture: &autoassess_landsurf_soilmoisture_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Land-Surface Soilmoisture Diagnostic"
        area: land_surface_soilmoisture
        control_model: IPSL-CM5A-LR
        exp_model: inmcm4
        obs_models: []
        start: 1997/12/01
        end: 2002/12/01
        climfiles_root: '/gws/nopw/j04/esmeval/autoassess_specific_files/files'  # on JASMIN

References
----------
* Dorigo, W.A., Wagner, W., Albergel, C., Albrecht, F.,  Balsamo, G., Brocca, L., Chung, D., Ertl, M., Forkel, M., Gruber, A., Haas, E., Hamer, D. P. Hirschi, M., Ikonen, J., De Jeu, R. Kidd, R.  Lahoz, W., Liu, Y.Y., Miralles, D., Lecomte, P. (2017).  ESA CCI Soil Moisture for improved Earth system understanding: State-of-the art and future directions. In Remote Sensing of Environment, 2017,  ISSN 0034-4257, https://doi.org/10.1016/j.rse.2017.07.001.

* Gruber, A., Scanlon, T., van der Schalie, R., Wagner, W., Dorigo, W. (2019). Evolution of the ESA CCI Soil Moisture Climate Data Records and their underlying merging methodology. Earth System Science Data 11, 717-739, https://doi.org/10.5194/essd-11-717-2019


Observations Data sets
----------------------

1999-2008 climatologies (seasonal means) from ESA ECV Soil Moisture Dataset v1.
Produced by the ESA CCI soil moisture project: https://www.esa-soilmoisture-cci.org/node/93


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
soil moisture median absolute error djf             0.0708               0.0708
soil moisture median absolute error mam             0.0665               0.0671
soil moisture median absolute error jja             0.0571               0.0564
soil moisture median absolute error son             0.0656               0.0661
===============================================     ================     ====================

.. figure:: /recipes/figures/autoassess_landsurface/Soilmoisture_Metrics.png
   :scale: 50 %
   :alt: Soilmoisture_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation
