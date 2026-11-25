.. _recipe_portrait:

Portrait plot
=============


Overview
--------
Portrait plots are a flexible way to visualize performance metrics for multiple
datasets and up to four references. In this recipe ``recipe_portrait_CMIP.yml``
the normalized Root Mean Squared Deviation (RMSD) of global mean seasonal
climatologies is calculated for a selection of CMIP models.
In the example recipe, for each variable up to two observation based datasets
are used as reference.
See :ref:`variables` for complete list of references.
The recipe uses preprocessor functions (distance metrics, global mean,
climate statistics) to calculate a scalar metric for each combination of
dataset, variable and reference, which is plotted by the ``portrait_plot.py``
diagnostic script.


User settings in recipe
-----------------------

By default cells are plotted for combinations of ``short_name``,
``dataset``, ``project`` and ``split``,
where ``split`` is an optional extra_facet for variables.
However, this can be customized using the ``x_by``,
``y_by``, ``group_by`` and ``split_by`` script settings.
For a complete and detailed list of settings, see the
:doc:`diagnostic documentation </api/esmvaltool.diag_scripts.portrait_plot>`.
While this allows very flexible use for any kind of data, there are some
limitations as well: The grouping (subplots) and normalization is always
applied along the x-axis.
With default settings this means normalizing all metrics for each variable
and grouping all datasets by project.

To plot distance metrics like RMSE, pearson R, bias etc. the
:func:`distance_metric <esmvalcore.preprocessor.distance_metric>` preprocessor
or custom diagnostics can be used.



.. _variables:

Variables and Datasets
------------------------

.. note::

   The recipe generally works for any variable that is preprocessed correctly.
   To use different preprocessors or reference datasets it could be useful
   to create different variable groups and link them with the same extra_facet
   like ``variable_name``. See recipe for examples. Listed below are the variables
   used to produce the example figure.


The following list shows which observational dataset is used as reference for
each variable in this recipe. All variables are atmospheric monthly means.
For 3D variables the selected pressure level is specified.

* clt (Ref1: ESACCI-CLOUD, Ref2: PATMOS-x)
* pr (Ref1: GPCP-V2.2)
* rlut, rsut (Ref1: CERES-EBAF)
* tas (Ref1: ERA-Interim, Ref2: NCEP-NCAR-R1)
* ua (200 hPa, Ref1: ERA-Interim, Ref2: NCEP-NCAR-R1)
* zg (500 hPa, Ref1: ERA-Interim, Ref2: NCEP-NCAR-R1)


References
----------

* Gleckler, P. J., K. E. Taylor, and C. Doutriaux, Performance metrics for climate models, J.
  Geophys. Res., 113, D06104, doi: 10.1029/2007JD008972 (2008).

* Righi, M., Eyring, V., Klinger, C., Frank, F., Gottschaldt, K.-D., Jöckel, P.,
  and Cionni, I.: Quantitative evaluation of ozone and selected climate parameters in a set of EMAC simulations,
  Geosci. Model Dev., 8, 733, doi: 10.5194/gmd-8-733-2015 (2015).


Example plots
-------------

.. _fig_portrait_plot:

.. figure:: /recipes/figures/portrait/portrait_plot.png
   :width: 90%
   :align: center


   Relative space-time root-mean-square deviation (RMSD) calculated from the climatological
   seasonal cycle of CMIP5 and CMIP6 simulations. A relative performance is displayed, with blue shading
   indicating better and red shading indicating worse performance than the median of all model results.
   A diagonal split of a grid square shows the relative error with respect to the reference data set.
