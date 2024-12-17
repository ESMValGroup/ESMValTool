.. _recipes_bock24acp:

Cloud properties and their projected changes in CMIP models with low to high climate sensitivity
================================================================================================

Overview
--------

The recipes recipe_bock24acp_*.yml reproduce figures (Fig. 3, 4, 6 and 7) from the publication `Bock and Lauer, 2024`_ investigating cloud properties and their projected changes in CMIP models with low to high climate sensitivity.

.. _`Bock and Lauer, 2024`: https://doi.org/10.5194/acp-24-1587-2024

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/clouds

    * recipe_bock24acp_fig3-4_maps.yml
    * recipe_bock24acp_fig6_zonal.yml
    * recipe_bock24acp_fig7_boxplots.yml

Diagnostics are stored in diag_scripts/

    Fig. 3 and 4:

    * clouds/clouds_ecs_groups_maps.py: Geographical maps of the multi-year annual means for group means of historical CMIP simulations from all three ECS groups.

    Fig. 6:

    * clouds/clouds_ecs_groups_zonal.py: Zonally averaged group means.

    Fig. 7:

    * clouds/clouds_ecs_groups_boxplots.py: Boxplots of relative changes for all groups.


User settings in recipe
-----------------------

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

   exclude_datasets: list of datasets which are not used for the statistics, default is ['MultiModelMean', 'MultiModelP5', 'MultiModelP95']
   group_by: list of 'variable_group's to have the order
   plot_type: 'zonal' and 'height' plots are available 

   *Optional settings (scripts)*

   filename_attach: attachment to the output files
   title: set title of figure
   y_range: set range of the y-axes


Variables
---------

* clt (atmos, monthly, longitude latitude time)
* clivi (atmos, monthly, longitude latitude time)
* clwvi (atmos, monthly, longitude latitude time)
* rlut (atmos, monthly, longitude latitude time)
* rsut (atmos, monthly, longitude latitude time)
* rlutcs (atmos, monthly, longitude latitude time)
* rsutcs (atmos, monthly, longitude latitude time)
* tas (atmos, monthly, longitude latitude time)


Observations and reformat scripts
---------------------------------

* CERES-EBAF (Ed4.2) - TOA radiation fluxes (used for calculation of
  the cloud radiative effects)

  *Reformat script:* cmorizers/data/formatters/datasets/ceres_ebaf.py


References
----------

* Bock, L. and Lauer, A.: Cloud properties and their projected changes in CMIP
  models with low to high climate sensitivity, Atmos. Chem. Phys., 24, 1587–1605,
  https://doi.org/10.5194/acp-24-1587-2024, 2024.


Example plots
-------------

.. _fig_bock24acp_1:
.. figure::  /recipes/figures/bock24acp/map_netcre.png
   :align:   center

   Geographical map of the multi-year annual mean net cloud radiative effect from 
   (a) CERES–EBAF Ed4.2 (OBS) and (b–d) group means of historical CMIP simulations
   from all three ECS groups (Fig. 4).

.. _fig_bock24acp_2:
.. figure::  /recipes/figures/bock24acp/zonal_diff_clt_ssp585.png
   :align:   center

   The upper panel show the zonally averaged group means of total cloud
   fraction from historical simulations (solid lines)
   and RCP8.5/SSP5-8.5 scenarios (dashed lines) for the three different ECS groups.
   Lower panels show the corresponding relative differences of all zonally
   averaged group means between the RCP8.5/SSP5-8.5 scenarios and the corresponding
   historical simulations. Shading indicates the 5 % and 95 % quantiles of the single
   model results (Fig. 6a).

.. _fig_bock24acp_3:
.. figure::  /recipes/figures/bock24acp/boxplot_ssp585_south_oc.png
   :align:   center

   Relative change (calculated as the difference between the scenario value and the
   historical value divided by the historical value) of total cloud fraction (clt),
   ice water path (iwp), liquid water path (lwp), and net cloud radiative effect
   (netcre) per degree of warming averaged over the Southern Ocean (30–65°S). In the
   box plot, each box indicates the range from the first
   quartile to the third quartile, the vertical line shows the median, and the
   whiskers the minimum and maximum values, excluding the outliers. Outliers are
   defined as being outside 1.5 times the interquartile range (Fig. 7b).

