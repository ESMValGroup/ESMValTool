.. _recipes_bock24acp:

Cloud properties and their projected changes in CMIP models with low to high climate sensitivity
================================================================================================

Overview
--------

The recipes recipe_bock24acp_*.yml reproduce figures from the publication Bock and Lauer, 2024.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/clouds

    * recipe_bock24acp_fig1_2_ecs_cloud_feedback.yml
    * recipe_bock24acp_fig3-4_maps.yml
    * recipe_bock24acp_fig6_zonal.yml
    * recipe_bock24acp_fig7_boxplots.yml

Diagnostics are stored in diag_scripts/

    Fig.1:

    * mlr/plot.py

    Fig.2:

    * clouds/clouds_ecs_scatterplot.py

    Fig. 3 and 4:

    * clouds/clouds_ecs_groups_maps.py: Geographical maps of the multi-year annual means for group means of historical CMIP simulations from all three ECS groups.

    Fig. 6:

    * clouds/clouds_ecs_groups_zonal.py: Zonally averaged group means.

    Fig. 7:

    * clouds/clouds_ecs_groups_boxplots.py: Boxplots of relative changes for all groups.


User settings in recipe
-----------------------

#. Script mlr/plot.py

   * :ref:`mlr/plot.py<api.esmvaltool.diag_scripts.mlr.plot>`

#. Script clouds/clouds_ecs_scatterplot.py

   *Required settings (scripts)*

   * ``file1:`` contains the feedback parameters for all cmip models from the 
     ``climate_metrics/feedback_parameters.py`` diagnostic
   * ``file2:`` contains the ECS values for all cmip models from the 
     ``climate_metrics/ecs.py`` diagnostic
   * ``output_file_name``: Set output file name.
   * ``xlabel:`` label on x-axis
   * ``ylabel:`` label on y-axis
   * ``x_range:`` Range for the x-axis of the plot
   * ``y_range:`` Range for the y-axis of the plot

   *Optional settings (scripts)*

   * ``seaborn_settings``, *dict*, optional: Options for
     :func:`seaborn.set_theme` (affects all plots).
   * ``print_corr``: Print and save Pearson correlation coefficient between all datasets at the
     end.  Requires identical shapes for all datasets. (default: False)

#. Script clouds_ecs_groups_maps.py 

   *Required settings (scripts)*

   reference: if true, a reference dataset is given within 'variable_group' equal 'OBS'

   *Optional settings (scripts)*

   plot_each_model: one figure for each single model


#. Script clouds/clouds_ecs_groups_zonal.py

   *Required settings (scripts)*

   group_by: list of 'variable_group's to have the order
   plot_type: 'zonal' and 'height' plots are available 

   *Optional settings (scripts)*

   filename_attach: attachment to the output files


#. Script clouds/clouds_ecs_groups_boxplots.py

   *Required settings (scripts)*

   exclude_datasets: list of datasets which are not used for the statistics,
                     default is ['MultiModelMean', 'MultiModelP5', 'MultiModelP95']
   group_by: list of 'variable_group's to have the order
   plot_type: 'zonal' and 'height' plots are available 

   *Optional settings (scripts)*

   filename_attach: attachment to the output files
   title: set title of figure
   y_range: set range of the y-axes


##. Script clouds_ipcc.ncl
#
#   See :ref:`here<clouds_ipcc.ncl>`.


Variables
---------

* clt (atmos, monthly, longitude latitude time)
* iwp (atmos, monthly, longitude latitude time)
* lwp (atmos, monthly, longitude latitude time)
* rlut (atmos, monthly, longitude latitude time)
* rsut (atmos, monthly, longitude latitude time)
* rlutcs (atmos, monthly, longitude latitude time)
* rsutcs (atmos, monthly, longitude latitude time)
* tas (atmos, monthly, longitude latitude time)


Observations and reformat scripts
---------------------------------

* CERES-EBAF - CERES TOA radiation fluxes (used for calculation of
  cloud forcing)

  *Reformat script:* cmorizers/data/formatters/datasets/ceres_ebaf.py


References
----------

* Bock, L. and Lauer, A.: Cloud properties and their projected changes in CMIP
  models with low to high climate sensitivity, Atmos. Chem. Phys., 24, 1587–1605,
  https://doi.org/10.5194/acp-24-1587-2024, 2024.


Example plots
-------------

.. _fig_bock24acp_1:

.. figure::  /recipes/figures/bock24acp/ecs_netcre.png
   :align:   center

   Scatterplot of the global mean net cloud feedback (x axis) and ECS (y axis) of the CMIP
   models, with a regression line, including the confidence interval of the regression of
   95 %. Dashed horizontal lines indicate separations of the three ECS groups (Fig. 1a).

.. _fig_bock24acp_2:
.. figure::  /recipes/figures/bock24acp/map_prediction_output___high_ECS_netcre.png
   :align:   center

   Geographical maps of net cloud feedback for high-ECS groups (Fig. 2a).

.. _fig_bock24acp_3:
.. figure::  /recipes/figures/bock24acp/map_netcre.png
   :align:   center

   Geographical map of the multi-year annual mean net cloud radiative effect from 
   (a) CERES–EBAF Ed4.2 (OBS) and (b–d) group means of historical CMIP simulations
   from all three ECS groups (Fig. 4).

.. _fig_bock24acp_4:
.. figure::  /recipes/figures/bock24acp/map_netcre.png
   :align:   center

   Each labelled subfigure contains two panels, namely the upper panels and lower
   panels. Upper panels show the zonally averaged group means of (a) total cloud
   fraction, (b) liquid water path, (c) ice water path, and (d) net, (e) shortwave,
   and (f) longwave cloud radiative effect from historical simulations (solid lines)
   and RCP8.5/SSP5-8.5 scenarios (dashed lines) for the three different ECS groups.
   The reference dataset CERES–EBAF Ed4.2 is shown as solid black lines in panels
   (d)–(f). Lower panels show the corresponding relative differences of all zonally
   averaged group means between the RCP8.5/SSP5-8.5 scenarios and the corresponding
   historical simulations. Shading indicates the 5 % and 95 % quantiles of the single
   model results (Fig. 6).

.. _fig_bock24acp_5:
.. figure::  /recipes/figures/bock24acp/map_netcre.png
   :align:   center

   Relative change (calculated as the difference between the scenario value and the
   historical value divided by the historical value) of total cloud fraction (clt),
   ice water path (iwp), liquid water path (lwp), and net cloud radiative effect
   (netcre) per degree of warming averaged over selected regions over the ocean.
   (a) Arctic (70–90∘ N), (b) Southern Ocean (30–65∘ S), (c) tropical ocean
   (30∘ N–30∘ S), (d) Pacific ITCZ (0–12∘ N, 135∘ E–85∘ W), and (e) the three
   stratocumulus regions of the southeast Pacific (10–30∘ S, 75–95∘ W), southeast
   Atlantic (10–30∘ S, 10∘ W–10∘ E), and northeast Pacific (15–35∘ N, 120–140∘ W)
   (see also Fig. 5). In the box plot, each box indicates the range from the first
   quartile to the third quartile, the vertical line shows the median, and the
   whiskers the minimum and maximum values, excluding the outliers. Outliers are
   defined as being outside 1.5 times the interquartile range (Fig. 7).

