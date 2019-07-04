.. _recipe_autoassess_stratosphere.rst:

Autoassess Stratosphere
=======================

Overview
--------

Polar night jet / easterly jet strengths are defined as the maximum / minimum wind
speed of the climatological zonal mean jet, and measure how realistic the zonal
wind climatology is in the stratosphere.

Extratropical temperature at 50hPa (area averaged poleward of 60 degrees) is important
for polar stratospheric cloud formation (in winter/spring), determining the amount of
heterogeneous ozone depletion simulated by models with interactive chemistry schemes.

The Quasi-Biennial Oscillation (QBO) is a good measure of tropical variability in the
stratosphere.  Zonal mean zonal wind at 30hPa is used to define the period and amplitude
of the QBO.

The tropical tropopause cold point (100hPa, 10S-10N) temperature is an important factor in
determining the stratospheric water vapour concentrations at entry point (70hPa, 10S-10N),
and this in turn is important for the accurate simulation of stratospheric chemistry and
radiative balance.

Age of air diagnoses how well the stratospheric meridional mass circulation is simulated,
important for the simulated thermal structure of the stratosphere and the distribution of
chemical species throughout the stratosphere.


Prior and current contributors
------------------------------
Met Office:

* Prior to May 2008: Neal Butchart
* May 2008 - May 2016: Steven C Hardiman
* Since May 2016: Alistair Sellar and Paul Earnshaw

ESMValTool:

* Since April 2018: Porting into ESMValTool by Valeriu Predoi


Developers
----------
Met Office:

* Prior to May 2008: Neal Butchart
* May 2008 - May 2016: Steven C Hardiman

ESMValTool:

* Since April 2018: Valeriu Predoi

Metrics and Diagnostics
-----------------------

Performance metrics:

* Polar night jet: northern hem (January) vs. ERA Interim (observational uncertainty from UKMO analysis)
* Polar night jet: southern hem (July) vs. ERA Interim (observational uncertainty from UKMO analysis)
* Easterly jet: southern hem (January) vs. ERA Interim (observational uncertainty from UKMO analysis)
* Easterly jet: northern hem (July) vs. ERA Interim (observational uncertainty from UKMO analysis)
* 50 hPa temperature: 60N-90N (DJF) vs. ERA Interim (observational uncertainty from MERRA)
* 50 hPa temperature: 60N-90N (MAM) vs. ERA Interim (observational uncertainty from MERRA)
* 50 hPa temperature: 90S-60S (JJA) vs. ERA Interim (observational uncertainty from MERRA)
* 50 hPa temperature: 90S-60S (SON) vs. ERA Interim (observational uncertainty from MERRA)
* QBO period at 30 hPa vs. ERA Interim (observational uncertainty from UKMO analysis)
* QBO amplitude at 30 hPa (westward) vs. ERA Interim (observational uncertainty from UKMO analysis)
* QBO amplitude at 30 hPa (eastward) vs. ERA Interim (observational uncertainty from UKMO analysis)
* 100 hPa equatorial temp (annual mean) vs. ERA Interim (observational uncertainty from MERRA)
* 100 hPa equatorial temp (annual cycle strength) vs. ERA Interim (observational uncertainty from MERRA)
* 70 hPa 10S-10N water vapour (annual mean) vs. MERRA (observational uncertainty from ERA Interim)

Diagnostics:

* Age of stratospheric air vs. observations from Andrews et al. (2001) and Engel et al. (2009)


Model Data
----------

===========================   ================== ============== ==============================================
Variable/Field name           realm              frequency      Comment
===========================   ================== ============== ==============================================
Eastward wind (ua)            Atmosphere         monthly mean   original stash: x-wind, no stash
Air temperature (ta)          Atmosphere         monthly mean   original stash: m01s30i204
Specific humidity (hus)       Atmosphere         monthly mean   original stash: m01s30i205
===========================   ================== ============== ==============================================

Inputs and usage
----------------
This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other AUtoassess metrics. An example below:

.. code-block:: yaml

    scripts:
      autoassess_strato_test_1: &autoassess_strato_test_1_settings
        script: autoassess/autoassess_area_base.py  # the base wrapper
        title: "Autoassess Stratosphere Diagnostic Metric"  # title
        area: stratosphere  # assesment area
        control_model: UKESM1-0-LL-hist  # control dataset name
        exp_model: UKESM1-0-LL-piCont  # experiment dataset name
        obs_models: [ERA-Interim]  # list to hold models that are NOT for metrics but for obs operations
        additional_metrics: [ERA-Interim]  # list to hold additional datasets for metrics
        start: 2004/12/01  # start date in native Autoassess format
        end: 2014/12/01  # end date in native Autoassess format


References
----------
* Andrews, A. E., and Coauthors, 2001: Mean ages of stratospheric air derived from in situ observations of CO2, CH4, and N2O. J. Geophys. Res.,   106 (D23), 32295-32314.
* Dee, D. P., and Coauthors, 2011: The ERA-Interim reanalysis: configuration and performance of the data assimilation system. Q. J. R. Meteorol.  Soc, 137, 553-597, doi:10.1002/qj.828.
* Engel, A., and Coauthors, 2009: Age of stratospheric air unchanged within uncertainties over the past 30 years. Nat. Geosci., 2, 28-31, doi:10  .1038/NGEO388.
* Rienecker, M. M., and Coauthors, 2011: MERRA: NASA’s Modern-Era Retrospective Analysis for Research and Applications. J. Climate, 24, 3624-364  8, doi:http://dx.doi.org/10.1175/JCLI-D-11-00015.1.


Observations Data sets
----------------------

ERA-Interim data (Dee et al., 2011) and MERRA data (Rienecker et al., 2011) can be obtained online from ECMWF and NASA respectively.  Monthly mean zonal mean U and T data are required.

Age of air data (Andrews et al., 2001; Engel et al., 2009) is as provided in age_of_air.py of the stratospheric area of auto_assess.

For UKMO analysis data, contact the Met Office.


Sample Plots and metrics
------------------------

===============================================     ================
Metric name                                         UKESM1-0-LL (historical)
                                                    value
===============================================     ================
Polar night jet: northern hem (January)             40.326
Polar night jet: southern hem (July)                84.867
Easterly jet: southern hem (January)                24.854
Easterly jet: northern hem (July)                   29.870
QBO period at 30 hPa                                41.500
QBO amplitude at 30 hPa (westward)                  27.383
QBO amplitude at 30 hPa (eastward)                  17.316
50 hPa temperature: 60N-90N (DJF)                   26.753
50 hPa temperature: 60N-90N (MAM)                   40.946
50 hPa temperature: 90S-60S (JJA)                   11.103
50 hPa temperature: 90S-60S (SON)                   23.299
100 hPa equatorial temp (annual mean)               15.292
100 hPa equatorial temp (annual cycle strength)      1.668
100 hPa 10Sto10N temp (annual mean)                 15.435
100 hPa 10Sto10N temp (annual cycle strength)        1.623
70 hPa 10Sto10N wv (annual mean)                     5.743
===============================================     ================


.. figure:: /recipes/figures/autoassess_stratosphere/metrics.png
   :scale: 50 %
   :alt: metrics.png

   Standard metrics plot


.. figure:: /recipes/figures/autoassess_stratosphere/t100_vs_q70.png
   :scale: 50 %
   :alt: t100_vs_q70.png

   Biases in tropical tropopause temperature (100hPa, 10S-10N) and lower stratospheric humidity (70hPa, 10S-10N)


.. figure:: /recipes/figures/autoassess_stratosphere/qbo_30hpa.png
   :scale: 50 %
   :alt: qbo_30hpa.png

   QBO at 30hPa comparison between UKESM1-0-LL (piControl and historical).
