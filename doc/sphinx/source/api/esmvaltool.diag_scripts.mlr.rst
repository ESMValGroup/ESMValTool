.. _api.esmvaltool.diag_scripts.mlr:

Machine Learning Regression (MLR) diagnostics
=============================================

This module provides various tools to create and evaluate MLR models for
arbitrary input variables.


Examples
--------

* :ref:`recipes_schlund20jgr`: Use Gradient Boosted Regression Tree (GBRT)
  algorithm to constrain projected Gross Primary Production (GPP) in RCP 8.5
  scenario using observations of process-based predictors.


Diagnostic scripts
------------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.mlr/evaluate_residuals
   esmvaltool.diag_scripts.mlr/main
   esmvaltool.diag_scripts.mlr/mmm
   esmvaltool.diag_scripts.mlr/plot
   esmvaltool.diag_scripts.mlr/postprocess
   esmvaltool.diag_scripts.mlr/preprocess
   esmvaltool.diag_scripts.mlr/rescale_with_emergent_constraint


Auxiliary scripts
-----------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.mlr/init
   esmvaltool.diag_scripts.mlr/custom_sklearn
   esmvaltool.diag_scripts.mlr/models
   esmvaltool.diag_scripts.mlr/models.gbr_base
   esmvaltool.diag_scripts.mlr/models.linear_base


.. _availableMLRModels:

Available MLR models
--------------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.mlr/models.gbr_sklearn
   esmvaltool.diag_scripts.mlr/models.gbr_xgboost
   esmvaltool.diag_scripts.mlr/models.gpr_sklearn
   esmvaltool.diag_scripts.mlr/models.huber
   esmvaltool.diag_scripts.mlr/models.krr
   esmvaltool.diag_scripts.mlr/models.lasso
   esmvaltool.diag_scripts.mlr/models.lasso_cv
   esmvaltool.diag_scripts.mlr/models.lasso_lars_cv
   esmvaltool.diag_scripts.mlr/models.linear
   esmvaltool.diag_scripts.mlr/models.rfr
   esmvaltool.diag_scripts.mlr/models.ridge
   esmvaltool.diag_scripts.mlr/models.ridge_cv
   esmvaltool.diag_scripts.mlr/models.svr
