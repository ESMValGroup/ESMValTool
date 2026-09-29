.. _recipe_cmip_hydrography:

CMIP ocean hydrography and density compensation
===============================================

The `hydrographic benchmark recipe <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/recipe_cmip_hydrographic_benchmark.yml>`_
compares CMIP6 ``thetao`` and ``so`` with WOA18 at 10, 100, 500, 1000, and
2000 m. It reports area-weighted historical-minus-WOA mean differences, root
mean square differences, and paired wet-cell coverage by depth and latitude
band. It also estimates piControl drift from annual fixed-depth fields and
from archived whole-ocean ``thetaoga``. The
`density compensation recipe <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/recipe_cmip_density_compensation.yml>`_
uses the same hydrographic comparison to attribute the density difference to
temperature and salinity with TEOS-10. Both recipes save numeric NetCDF outputs
and figures. With two or more models, they also save a comparison using cells
that are valid in WOA and every model.

Input and interpretation
------------------------

* The historical comparison uses 1981–2010 annual model means. The intended
  reference is the WOA18 ``decav81B0`` 1981–2010 climate normal, supplied as
  one climatological field with representative year 2000 in OBS6. Confirm that
  the local WOA ``thetao`` and ``so`` files were both derived from this normal:
  the OBS6 filenames and global attributes do not identify the original WOA
  source files or prove their averaging period.
* WOA ``thetao`` is derived from in-situ temperature by the OBS6 CMORizer. The
  diagnostic converts it to potential temperature with WOA practical salinity
  before comparing it with CMIP6 ``thetao``.
* Choose each model's own contiguous piControl segment in the hydrographic
  recipe. The example years are for ACCESS-CM2 and should be checked against
  the local archive. Add one historical and one control dataset entry per
  model in the shared list used by ``thetao_model`` and ``so_model``, and add
  its ``thetaoga`` control entry to ``whole_ocean_drift``. For the density
  recipe, add historical entries only.
* Differences against WOA are climatological reference differences. Fixed-depth
  area means are neither heat content nor volume-weighted metrics. WOA wet-cell
  coverage is not observational uncertainty. The density attribution is exact
  for the equation of state evaluated on the two climatological fields; it
  does not recover the time mean of nonlinear density from monthly fields.
* The control slope comes from the chosen 40-year segment. It is sensitive to
  internal variability and is not, by itself, a long-term drift estimate.
  Compare longer or multiple control segments before using it for model
  assessment. Fixed-depth differences can also reflect displaced water masses;
  sections or isopycnal diagnostics help with that interpretation.

The diagnostic scripts require ``gsw`` from the ESMValTool ``all`` optional
dependency set. Confirm the WOA source filenames and run both complete recipes
on the target archive before interpreting scientific results.

Gadi validation and example figures
-----------------------------------

Both recipes completed with ACCESS-CM2 historical 1981–2010, piControl
0950–0989, and the local WOA18 OBS6 temperature and salinity files on Gadi.
The run used ``analysis3-26.09`` (ESMValTool and ESMValCore 2.15.0, ``gsw``
3.6.23) and 28 Dask distributed workers for preprocessing. The figures below
come from the recipe diagnostics, not the companion notebooks. Provenance and
numeric metrics are saved with the recipe output. The original WOA source
filenames could not be recovered from the OBS6 file metadata, so the
1981–2010 reference-period identity remains to be confirmed before using
these results as a published model benchmark.

.. list-table:: Example global, area-weighted ACCESS-CM2 minus local WOA18 values
   :header-rows: 1

   * - Depth (m)
     - Potential temperature (°C)
     - Practical salinity (1e-3)
     - In-situ density (kg m⁻³)
     - Paired WOA-valid area
   * - 10
     - +0.016
     - −0.287
     - −0.208
     - 99.3%
   * - 100
     - +0.187
     - −0.274
     - −0.259
     - 97.6%
   * - 500
     - +2.336
     - −0.006
     - −0.418
     - 97.4%
   * - 1000
     - +2.854
     - +0.065
     - −0.393
     - 96.8%
   * - 2000
     - +1.370
     - +0.038
     - −0.176
     - 94.1%

The strong 500–1000 m warm difference contributes to a lower model density.
At 10 m, the global mean temperature difference is small despite a spatial
root mean square difference of 1.19 °C. The northern extratropical salinity
root mean square difference is much larger than the global value and merits
regional investigation. The archived ``thetaoga`` control slope over
0950–0989 is +0.060 °C per century; this short-window estimate is not a
long-term drift constraint.

.. figure:: /recipes/figures/cmip_hydrography/hydrography_ACCESS-CM2_thetao_reference.png
   :width: 90%

   ACCESS-CM2 historical-minus-WOA18 potential temperature by region and depth.
   The right panel is spatial, not temporal, root mean square difference.

.. figure:: /recipes/figures/cmip_hydrography/hydrography_ACCESS-CM2_so_reference.png
   :width: 90%

   Practical salinity differences and spatial root mean square differences.

.. figure:: /recipes/figures/cmip_hydrography/density_compensation_ACCESS-CM2.png
   :width: 90%

   Density difference and the area-and-effect-weighted local cancellation
   fraction of temperature and salinity contributions.

.. figure:: /recipes/figures/cmip_hydrography/hydrography_ACCESS-CM2_thetao_drift.png
   :width: 90%

   piControl temperature anomalies relative to 0950 at five sampled depths.
   Each row is a separate depth level, not a continuous vertical section.

.. figure:: /recipes/figures/cmip_hydrography/hydrography_ACCESS-CM2_so_drift.png
   :width: 90%

   piControl practical salinity anomalies at the same sampled depths.

.. figure:: /recipes/figures/cmip_hydrography/hydrography_ACCESS-CM2_global_mean_drift.png
   :width: 80%

   Archived whole-ocean ``thetaoga`` anomaly and fitted 40-year control slope.
