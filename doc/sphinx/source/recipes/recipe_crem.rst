Cloud Regime Error Metric (CREM)
================================

Overview
--------

The radiative feedback from clouds remains the largest source of uncertainty
in determining the climate sensitivity. Traditionally, cloud has been
evaluated in terms of its impact on the mean top of atmosphere fluxes.
However it is quite possible to achieve good performance on these criteria
through compensating errors, with boundary layer clouds being too reflective
but having insufficient horizontal coverage being a common example (e.g.,
Nam et al., 2012). Williams and Webb (2009) (WW09) propose a Cloud Regime
Error Metric (CREM) which critically tests the ability of a model to
simulate both the relative frequency of occurrence and the radiative
properties correctly for a set of cloud regimes determined by the daily
mean cloud top pressure, cloud albedo and fractional coverage at each
grid-box. WW09 describe in detail how to calculate their metrics and we
have included the CREMpd metric from their paper in ESMValTool, with clear
references in the lodged code to tables in their paper. This has been
applied to those CMIP5 models who have submitted the required diagnostics
for their AMIP simulation (see Figure 8 below). As documented by WW09, a
perfect score with respect to ISCCP would be zero. WW09 also compared
MODIS/ERBE to ISCCP in order to provide an estimate of observational
uncertainty. This was found to be 0.96 and this is marked on Figure 8,
hence a model with a CREM similar to this value could be considered to have
an error comparable with observational uncertainty, although it should be
noted that this does not necessarily mean that the model lies within the
observations for each regime. A limitation of the metric is that it requires
a model to be good enough to simulate each regime. If a model is that poor
that the simulated frequency of occurrence of a particular regime is zero,
then a NaN will be returned from the code and a bar not plotted on the
figure for that model.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_williams09climdyn_CREM.yml

Diagnostics are stored in diag_scripts/crem/

* ww09_esmvaltool.py



User settings
-------------

None.


Variables
---------

* albisccp (atmos, daily mean, longitude latitude time)
* cltisccp (atmos, daily mean, longitude latitude time)
* pctisccp (atmos, daily mean, longitude latitude time)
* rlut (atmos, daily mean, longitude latitude time)
* rlutcs (atmos, daily mean, longitude latitude time)
* rsut (atmos, daily mean, longitude latitude time)
* rsutcs (atmos, daily mean, longitude latitude time)
* sic/siconc (seaice, daily mean, longitude latitude time)
* snc (atmos, daily mean, longitude latitude time)

If snc is not available then snw can be used instead. For AMIP simulations,
sic/siconc is often not submitted as it a boundary condition and effectively
the same for every model. In this case the same daily sic data set can be
used for each model.

**Note: in case of using sic/siconc data from a different model (AMIP), it has to
be checked by the user that the calendar definitions of all data sets are
compatible, in particular whether leap days are included or not.**



Observations and reformat scripts
---------------------------------

All observational data have been pre-processed and included within the
routine. These are ISCCP, ISCCP-FD, MODIS, ERBE. No additional observational
data are required at runtime.



References
----------

* Nam, C., Bony, S., Dufresne, J.-L., and Chepfer, H.: The 'too few, too bright'
  tropical low-cloud problem in CMIP5 models, Geophys. Res. Lett., 39, L21801,
  doi: 10.1029/2012GL053421, 2012.
* Williams, K.D. and Webb, M.J.: A quantitative performance assessment of
  cloud regimes in climate models. Clim. Dyn. 33, 141-157, doi:
  10.1007/s00382-008-0443-1, 2009.


Example plots
-------------

.. figure:: /recipes/figures/crem/crem_error_metric.png
   :width: 10cm
   :alt: xxxxx

   Cloud Regime Error Metrics (CREMpd) from William and Webb (2009) applied
   to those CMIP5 AMIP simulations with the required data in the archive. A
   perfect score with respect to ISCCP is zero; the dashed red line is an
   indication of observational uncertainty.
