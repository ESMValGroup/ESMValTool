.. _recipe_autoassess_landsurface_surfrad.rst:

Land-surface Surface Radiation - Autoassess diagnostics
=======================================================

Overview
--------


Prior and current contributors
------------------------------
Met Office:

* Prior to April 2018: John Edwards, Paul Earnshaw

ESMValTool:

* Since April 2018: Porting into ESMValTool by Alistair Sellar and Valeriu Predoi


Developers
----------
Met Office:

* John Edwards


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

* median absolute error (model minus observations) net surface shortwave (SW) radiation
* median absolute error (model minus observations) net surface longwave (LW) radiation

Metrics are calculated using model and observation multi-year climatologies (seasonal means) 
for meteorological seasons:
* December-January-February (djf)
* March-April-May (mam)
* June-July-August (jja)
* September-October-November (son)
* Annual mean (ann)


Diagnostics:

* Metrics plot comparing control and experiment


Model Data
----------

========================================================================= ================== ============== ==============================================
Variable/Field name                                                       realm              frequency      Comment
========================================================================= ================== ============== ==============================================
Surf SW net all sky (rsns)                                                Atmosphere         monthly mean
Surf LW net all sky (rlns)                                                Atmosphere         monthly mean
Percentage of the Grid Cell Occupied by Land (Including Lakes) (sftlf)    mask               fixed
========================================================================= ================== ============== ==============================================

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
      autoassess_landsurf_surfrad: &autoassess_landsurf_surfrad_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Land-Surface Diagnostic Surfrad Metric"
        area: land_surface_surfrad
        control_model: UKESM1-0-LL
        exp_model: UKESM1-0-LL
        obs_models: [CERES-EBAF]
        obs_type: obs4mips
        start: 1997/12/01
        end: 2002/12/01


References
----------
Loeb, N. G., D. R. Doelling, H. Wang, W. Su, C. Nguyen, J. G. Corbett, L. Liang, C. Mitrescu, F. G. Rose, and S. Kato, 2018: Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) Top-of-Atmosphere (TOA) Edition-4.0 Data Product. J. Climate, 31, 895-918, doi: 10.1175/JCLI-D-17-0208.1.

Kato, S., F. G. Rose, D. A. Rutan, T. E. Thorsen, N. G. Loeb, D. R. Doelling, X. Huang, W. L. Smith, W. Su, and S.-H. Ham, 2018: Surface irradiances of Edition 4.0 Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) data product, J. Climate, 31, 4501-4527, doi: 10.1175/JCLI-D-17-0523.1


Observations Data sets
----------------------

2000-2009 climatologies (seasonal means) from CERES-EBAF Ed2.7.


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
Net SW median absolute error ann                    4.88                 4.93
Net LW median absolute error ann                    3.98                 3.81
Net SW median absolute error djf                    6.51                 6.69
Net LW median absolute error djf                    5.27                 5.23
Net SW median absolute error mam                    4.31                 4.68
Net LW median absolute error mam                    4.51                 4.46
Net SW median absolute error jja                    6.47                 6.11
Net LW median absolute error jja                    5.37                 5.70
Net SW median absolute error son                    5.60                 5.50
Net LW median absolute error son                    4.77                 4.69
===============================================     ================     ====================

.. figure:: /recipes/figures/autoassess_landsurface/Surfrad_Metrics.png
   :scale: 50 %
   :alt: Surfrad_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation
