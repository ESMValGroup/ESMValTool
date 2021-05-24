.. _recipe_autoassess_landsurface_permafrost.rst:

Land-surface Permafrost - Autoassess diagnostics
================================================

Overview
--------


Prior and current contributors
------------------------------
Met Office:

+ Prior to May 2016:
* Since May 2016: Alistair Sellar and Paul Earnshaw

ESMValTool:

* Since April 2018: Porting into ESMValTool by Alistair Sellar and Valeriu Predoi


Developers
----------
Met Office:


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


Diagnostics:

* 


Model Data
----------

===========================   ================== ============== ==============================================
Variable/Field name           realm              frequency      Comment
===========================   ================== ============== ==============================================
Near-Surface Air Temperature (tas)            Atmosphere         monthly mean
Temperature of Soil (tsl)      Land         monthly mean
Percentage of the Grid Cell Occupied by Land (Including Lakes) (sftlf)       mask         fixed
Total Water Content of Soil Layer (mrsol) Emon    monthly mean          CMIP5: ??
===========================   ================== ============== ==============================================

The recipe takes as input a control model and experimental model, comparisons being made
with these two CMIP models; additionally it can take observational datas input.

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
* Andrews, A. E., and Coauthors, 2001: Mean ages of stratospheric air derived from in situ observations of CO2, CH4, and N2O. J. Geophys. Res.,   106 (D23), 32295-32314.
* Dee, D. P., and Coauthors, 2011: The ERA-Interim reanalysis: configuration and performance of the data assimilation system. Q. J. R. Meteorol.  Soc, 137, 553-597, doi:10.1002/qj.828.
* Engel, A., and Coauthors, 2009: Age of stratospheric air unchanged within uncertainties over the past 30 years. Nat. Geosci., 2, 28-31, doi:10  .1038/NGEO388.

Observations Data sets
----------------------

None used in this diagnostic.

Sample Plots and metrics
------------------------
Below is a set of metrics for  UKESM1-0-LL (historical data); the table
shows a comparison made between running ESMValTool on CMIP6 CMORized
netCDF data freely available on ESGF nodes and the run made using native
Autoassess performed at the Met Office using the pp output of the model.

===============================================     ================     ====================
Metric name                                         UKESM1-0-LL;         UKESM1-0-LL;
                                                    CMIP6: AERmonZ;      pp files;
                                                    historical, ESGF     historical, u-bc179
===============================================     ================     ====================
metrics here
===============================================     ================     ====================

Some notes on the comparison runs here (location of runs, ideally, path to results too)

Add figures here:
.. figure:: /recipes/figures/autoassess_stratosphere/metrics.png
   :scale: 50 %
   :alt: metrics.png

   Standard metrics plot comparing standard metrics from UKESM1-0-LL and HadGEM3-GC31.


