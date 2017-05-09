Performance metrics for essential climate parameters
====================================================

Overview
--------

The goal is to create a standard namelist for the calculation of performance metrics to quantify the ability of the models to reproduce the
climatological mean annual cycle for selected "Essential Climate Variables" (ECVs) plus some additional corresponding diagnostics and plots to better
understand and interpret the results. The namelist can be used to calculate performance metrics at different vertical levels (e.g., 5, 30, 200, 850
hPa as in Gleckler et al., 2008) and in four regions (global, tropics 20°N-20°S, northern extratropics 20°-90°N, southern extratropics 20°-90°S). As
an additional reference, we consider the Righi et al. (2015) paper.

Available Namelists and Diagnostics
-----------------------------------

Namelists are stored in nml/

* namelist_perfmetrics_CMIP5.xml

Diagnostics are stored in diag_scripts/

* perfmetrics_grading.ncl: calculates grades according to a given metric, with different options for normalization. It requires fields
precalculated by perfmetrics_main.ncl.
* perfmetrics_grading_collect.ncl: collects results from metrics previously calculated by perfmetrics_grading.ncl and passes them to the plotting functions.
* perfmetrics_main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields from input monthly 2-d or 3-d ("T2M", "T3Ms") data. The calculated fields can be also plotted as difference w.r.t. a given reference model. They are also used as input to calculate grading metrics (see perfmetrics_grading.ncl).
* perfmetrics_taylor.ncl: calculates grades according to a given metric, with different options for normalization. It requires fields precalculated by perfmetrics_main.ncl.
* perfmetrics_taylor_collect.ncl: collects results from metrics previously calculated by perfmetrics_taylor.ncl and passes them to the plotting functions.

User settings
-------------

TBD


Variables
---------

TBD


Observations and Reformat Scripts
---------------------------------

TBD



References
----------

TBD


Example plots
-------------

TBD

.. figure:: ../../source/namelists/figures/TBDNAMELIST/TBDFIG.png
   :scale: 50 %
   :alt: xxxx
   
   CAPTION CAN GO HERE














