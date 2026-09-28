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

* The historical comparison uses 1981–2010 annual model means and the WOA18
  ``decav81B0`` 1981–2010 climate normal. WOA supplies one climatological
  field, represented by year 2000 in OBS6. Confirm that the local WOA
  ``thetao`` and ``so`` files were both derived from this normal.
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

The diagnostic scripts require ``gsw`` from the ESMValTool ``all`` optional
dependency set. Confirm the WOA source filenames and run both complete recipes
on the target archive before interpreting scientific results.
