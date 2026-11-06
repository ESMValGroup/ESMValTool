.. _recipes_rf_shap:

Random-forest / SHAP diagnosis of collocated ESM fields
=======================================================

Overview
--------

``rf_shap`` learns a mapping :math:`Y = f(X)` between collocated Earth-system
model (or observational) predictor cubes and a target cube, then explains the
mapping with tree-ensemble importance, SHAP, and related plots. Predictions and
SHAP attributions are restored to ``(time, latitude, longitude)`` NetCDF so they
fit the usual ESMValTool provenance and plotting workflow.

This diagnostic is developed in AI4PEX WP4 and tracked upstream in
`ESMValTool issue #4609 <https://github.com/ESMValGroup/ESMValTool/issues/4609>`_.

Available recipes and diagnostics
---------------------------------

Recipes are stored in ``esmvaltool/recipes/``

* ``recipe_rf_shap_demo.yml`` — synthetic cubes (install / CI smoke test)
* ``recipe_rf_shap_gpp.yml`` — CMIP6 GPP versus FluxCom/ERA5 drivers (land example)
* Additional overlay recipes (cloud, ocean, generic tas) may follow in later PRs

Diagnostics are stored in ``esmvaltool/diag_scripts/rf_shap/``

* ``diagnostic.py``: ExtraTrees (or related) fit, metrics, SHAP, optional ALE
* ``compare.py``: compare two pipeline explanations (e.g. CMIP vs reference)

User settings in recipe
-----------------------

#. Script ``rf_shap/diagnostic.py``

 *Required settings for script*

 * ``rf_shap``: nested configuration for transforms, split, sampling, model, and
   explain options (see demo recipe for a minimal working block)

 *Optional settings for script*

 * ``use_synthetic``: if true, ignore input files and build in-memory cubes
 * ``target`` / ``predictors``: short names selecting which preprocessed variables
   are :math:`Y` and :math:`X`
 * ``caption``: figure / NetCDF provenance caption

 *Required settings for variables*

 * At least one target field and one or more predictor fields on a shared grid
   after the preprocessor (demo recipe needs none)

 *Optional settings for preprocessor*

 * ``regrid`` to a common grid (e.g. ``1x1``) before the diagnostic

Variables
---------

* Demo: synthetic predictors / target (no CMOR variables)
* GPP recipe: ``gpp`` (land, monthly), drivers ``tas``, ``hurs``, ``rsds``,
  ``mrsos``, ``pr``

Observations and reformat scripts
---------------------------------

* FluxCom GPP + ERA5 drivers can be staged as OBS6 with the helper scripts in
  the development repository (``prepare_obs6_gpp.py``). Official obs4MIPs /
  OBS6 entries should be preferred once available in the shared OBS tree.

References
----------

* AI4PEX project (https://www.ai4pex.org/)
* Lundberg, S. M., and Lee, S.-I.: A unified approach to interpreting model
  predictions, NeurIPS, 2017.
* Eyring, V., et al.: ESMValTool v2.0, Geosci. Model Dev., 13, 1179–1197, 2020.

Example plots
-------------

.. _fig_rf_shap_1:
.. figure:: /recipes/figures/rf_shap/shap_summary.png
   :align: center
   :width: 70%

   SHAP summary for the synthetic demo recipe (install smoke test).

.. _fig_rf_shap_2:
.. figure:: /recipes/figures/rf_shap/rf_importance.png
   :align: center
   :width: 60%

   Impurity importance from the ExtraTrees model (synthetic demo).

.. _fig_rf_shap_3:
.. figure:: /recipes/figures/rf_shap/pred_vs_obs.png
   :align: center
   :width: 60%

   Predicted versus target scatter for the held-out synthetic samples.
