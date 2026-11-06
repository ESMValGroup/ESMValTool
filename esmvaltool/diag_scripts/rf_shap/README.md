# rf_shap diagnostic

Python diagnostic for ExtraTrees / SHAP evaluation of collocated ESM and
observation fields (AI4PEX WP4, ESMValTool issue #4609).

## Runtime layout

The scientific kernel is **vendored** under ``_vendor/rf_shap`` so this
diagnostic runs inside an ESMValTool checkout without a separate pip install of
the development repository.

Optional extras (install in the ESMValTool environment if you enable them in
the recipe):

* ``shap`` — TreeSHAP / LinearExplainer backends
* ``fasttreeshap`` — faster TreeSHAP when available

Core stack already expected by ESMValTool: ``numpy``, ``pandas``, ``xarray``,
``scikit-learn``, ``matplotlib``, ``netCDF4``.

## Recipes

* ``recipe_rf_shap_demo.yml`` — synthetic cubes (no climate files)
* ``recipe_rf_shap_gpp.yml`` — CMIP6 GPP vs FluxCom/ERA5 drivers
