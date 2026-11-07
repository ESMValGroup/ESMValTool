.. _recipe_russell18jgr_cmip7:

Southern Ocean metrics with CMIP7
=================================

Overview
--------

``recipe_russell18jgr_cmip7.yml`` applies the diagnostics described in
:doc:`recipe_russell18jgr` to CMIP7 historical simulations for 1986--2005.
It is a separate companion to the original CMIP5 recipe and
:doc:`recipe_russell18jgr_cmip6`.
The three recipes retain the same diagnostic methods.

The selection reflects published fields available in the Gadi ``im55`` and
``ib12`` archives at the time of testing.
These model generations are not interchangeable with the paper's CMIP5 models.
The recipe includes 13 diagnostic groups; Figures 9a--c are omitted.

.. list-table:: Model selection
   :header-rows: 1
   :widths: 25 20 55

   * - CMIP7 model
     - Member
     - Diagnostic coverage
   * - ACCESS-ESM1-6
     - ``r1i1p1f1``
     - All included groups except Figure 8 (pH).
   * - CanESM6-0-MR
     - ``r1i1p1f1``
     - Figures 1, 2, 3b, 3b-2, 4 (velocity only), 5 and 5g.
   * - UKESM1-3-LL
     - ``r1i1p1f2``
     - All included groups except Figures 5g and 7i.

Historical ``ph`` was unavailable for ACCESS-ESM1-6 and CanESM6-0-MR.
CanESM6-0-MR also lacked ``fgco2`` and ``volcello`` in the inspected archive.
UKESM1-3-LL lacked ``areacello``, so its area-integrated diagnostics
are omitted.
The CanESM6-0-MR fixed cell areas are published under ``1pctCO2`` on the same
``g147`` ocean grid used by its historical sea ice.
Only that fixed geometric field is taken from ``1pctCO2``.

Only ACCESS-ESM1-6 has the complete archived historical carbon-flux and
cell-area inputs for Figures 9b/c in this selection.
The recipe therefore omits all Figure 9 regressions; it cannot reproduce the
paper's multi-model relationships or significance from this coverage.

Input conventions
-----------------

CMIP7 variables use the ``atmos``, ``ocean``, ``seaIce`` and ``ocnBgchem``
tables, explicit branding suffixes, and ``region: glb``.
Monthly fields specify ``frequency: mon`` and fixed areas ``frequency: fx``.
The recipe selects each field's published grid explicitly.
Atmospheric wind stress is land masked by the preprocessor; ocean fields
retain their published masks.

.. list-table:: Published grids
   :header-rows: 1
   :widths: 25 20 20 20 15

   * - Model
     - Atmospheric stress
     - Ocean tracers / sea ice
     - Zonal velocity / stress
     - Meridional velocity
   * - ACCESS-ESM1-6
     - ``g109``
     - ``g118``
     - ``g116``
     - ``g116``
   * - CanESM6-0-MR
     - ``g150``
     - ``g147``
     - ``g148``
     - ``g149``
   * - UKESM1-3-LL
     - ``g108``
     - ``g126``
     - ``g125``
     - ``g124``

No additional regridding is introduced.
The diagnostics use the Matplotlib colour cycle when a CMIP7 model style
file is unavailable.

Scientific interpretation
-------------------------

The caveats in :doc:`recipe_russell18jgr_cmip6` also apply here.
Figure 1 uses atmospheric ``tauu`` while Figure 2 uses ocean surface ``tauuo``.
Figure 5 uses March for the summer extent, following the legacy recipe;
the paper uses February.

Both ACCESS-ESM1-6 and UKESM1-3-LL supply monthly ``volcello``.
Figures 4 and 6 average these volumes and other input fields before computing
transport, retaining the original climatological method.
In general, the product of time means is not the mean of the product.
These calculations omit monthly covariance and changes in density-layer
membership, and use tracer-cell volumes to approximate velocity-face areas.
The explicit grid identifiers document the staggering; they do not establish
that pairing velocity and tracer cells by array index is conservative.

The outputs are legacy climatological transport estimates, not validated
model-native monthly volume or heat transports.
A small effect of volume variability on net transport for one model does not
validate individual density layers or heat transport, or establish the same
result for another model.
A successful recipe run checks execution and output consistency; it does not
remove these scientific limitations.

Gadi execution check
--------------------

The complete recipe passed on Gadi using ``analysis3-26.08`` and ESMValCore
2.15.0 (PBS job ``183400498``, 7 November 2026).
All 35 tasks completed, producing 17 figures and 30 diagnostic NetCDF files
across 13 groups, with model coverage and provenance checked.
All 32 preprocessed files retaining time contain the full 240 monthly samples.
Sea-ice and cell-area coordinates match exactly for ACCESS and CanESM.
The scientific limitations above still apply.

Example output
--------------

.. figure:: /recipes/figures/russell18jgr/cmip7/Fig1_page1.png
   :width: 60%

   Annual mean atmospheric zonal wind stress for UKESM1-3-LL and
   ACCESS-ESM1-6, historical 1986--2005, generated on Gadi.
   The other diagnostic types are illustrated in :doc:`recipe_russell18jgr`.

References
----------

* Russell, J. L., et al. (2018), Metrics for the Evaluation of the Southern
  Ocean in Coupled Climate Models and Earth System Models,
  Journal of Geophysical Research: Oceans, 123, 3120--3143,
  https://doi.org/10.1002/2017JC013461.
