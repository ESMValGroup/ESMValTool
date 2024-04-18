.. _recipes_bock24acp:

Cloud properties and their projected changes in CMIP models with low to high climate sensitivity
================================================================================================

Overview
--------

The recipes recipe_bock24acp_*.yml reproduce figures from the publication Bock and Lauer, 2024.


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

   none

   *Optional settings (scripts)*

   * title_key:
   * plot_each_model:

   *Required settings (variables)*

   none

   *Optional settings (variables)*

   none

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

.. _fig_bock24acp_4:
.. figure::  /recipes/figures/bock20jgr/tas_Global_CMIP6_historical_anom_1850-2014.png
   :align:   center

   Geographical map of the multi-year annual mean net cloud radiative effect from 
   (a) CERES–EBAF Ed4.2 (OBS) and (b–d) group means of historical CMIP simulations
   from all three ECS groups (Fig. 4).

.. _fig_bock24acp_6:
.. figure::  /recipes/figures/bock20jgr/tas_Global_CMIP6_historical_anom_1850-2014.png
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

.. _fig_bock24acp_7:
.. figure::  /recipes/figures/bock20jgr/tas_Global_CMIP6_historical_anom_1850-2014.png
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

