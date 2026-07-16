.. _recipes_impact:

Quick insights for climate impact researchers
=============================================

Overview
--------

Many impact researchers do not have the time and finances to use a large
ensemble of climate model runs for their impact analysis. To get an idea of the
range of impacts of climate change it also suffices to use a small number of
climate model runs. In case a system is only sensitive to annual temperature,
one can select a run with a high change and one with a low change of annual
temperature, preferably both with a low bias.

This recipe calculates the bias with respect to observations, and the change
with respect to a reference period, for a wide range of (CMIP) models. These
metrics are tabulated and also visualized in a diagram.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_impact.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * impact/bias_and_change.py: tabulate and visualize bias and change.


User settings in recipe
-----------------------

#. Script ``impact.py``

   *Required settings for variables*

   * tag: ``'model'`` or ``'observations'``, so the diagnostic script knows which datasets to use for the bias calculation. This must be specified for each dataset.

   *Optional settings for preprocessor*

   * Region and time settings (both for the future and reference period) can be changed at will.


Variables
---------

* tas (atmos, mon, longitude latitude time)
* pr (atmos, mon, longitude latitude time)
* any other variables of interest


Observations and reformat scripts
---------------------------------

* ERA5 data can be used via the native6 project.

References
----------

* None

Example plots
-------------

.. _fig_impact_1:
.. figure::  /recipes/figures/impact/bias_vs_change.png
   :align:   center

   Bias and change for each variable

.. raw:: html
    <embed>

        <style>#T_caeac_ .index_name{text-align:right}#T_caeac_
        .row_heading{text-align:right}#T_caeac_ td{padding:3px
        25px}#T_caeac_row0_col0,#T_caeac_row16_col1,#T_caeac_row16_col3{background-color:#bd1726;color:#f1f1f1}#T_caeac_row0_col1{background-color:#fffdbc;color:#000}#T_caeac_row0_col2,#T_caeac_row9_col3{background-color:#f7844e;color:#000}#T_caeac_row0_col3,#T_caeac_row10_col1{background-color:#e14430;color:#f1f1f1}#T_caeac_row1_col0{background-color:#e95538;color:#000}#T_caeac_row1_col1{background-color:#fdc171;color:#000}#T_caeac_row1_col2{background-color:#f99153;color:#000}#T_caeac_row15_col0,#T_caeac_row1_col3{background-color:#e24731;color:#f1f1f1}#T_caeac_row2_col0{background-color:#dc3b2c;color:#f1f1f1}#T_caeac_row10_col3,#T_caeac_row2_col1{background-color:#ea5739;color:#000}#T_caeac_row2_col2{background-color:#fed07e;color:#000}#T_caeac_row2_col3{background-color:#f88c51;color:#000}#T_caeac_row10_col0,#T_caeac_row3_col0{background-color:#b10b26;color:#f1f1f1}#T_caeac_row3_col1,#T_caeac_row8_col2,#T_caeac_row8_col3,#T_caeac_row9_col0{background-color:#feffbe;color:#000}#T_caeac_row3_col2{background-color:#fdb163;color:#000}#T_caeac_row11_col1,#T_caeac_row3_col3{background-color:#d93429;color:#f1f1f1}#T_caeac_row4_col0{background-color:#ab0626;color:#f1f1f1}#T_caeac_row4_col1,#T_caeac_row7_col2{background-color:#f67c4a;color:#000}#T_caeac_row4_col2{background-color:#fca55d;color:#000}#T_caeac_row15_col1,#T_caeac_row4_col3{background-color:#fa9857;color:#000}#T_caeac_row13_col3,#T_caeac_row5_col0{background-color:#ed5f3c;color:#000}#T_caeac_row5_col1{background-color:#dd3d2d;color:#f1f1f1}#T_caeac_row5_col2{background-color:#fdb365;color:#000}#T_caeac_row12_col1,#T_caeac_row5_col3,#T_caeac_row6_col3{background-color:#f67a49;color:#000}#T_caeac_row11_col2,#T_caeac_row6_col0{background-color:#f47044;color:#000}#T_caeac_row6_col1{background-color:#fba35c;color:#000}#T_caeac_row14_col0,#T_caeac_row17_col1,#T_caeac_row17_col3,#T_caeac_row6_col2,#T_caeac_row8_col0{background-color:#a50026;color:#f1f1f1}#T_caeac_row13_col0,#T_caeac_row7_col0,#T_caeac_row8_col1{background-color:#ad0826;color:#f1f1f1}#T_caeac_row16_col0,#T_caeac_row7_col1{background-color:#b50f26;color:#f1f1f1}#T_caeac_row7_col3{background-color:#d62f27;color:#f1f1f1}#T_caeac_row9_col1{background-color:#f57748;color:#000}#T_caeac_row9_col2{background-color:#ec5c3b;color:#000}#T_caeac_row10_col2{background-color:#e75337;color:#000}#T_caeac_row11_col0{background-color:#e54e35;color:#000}#T_caeac_row11_col3,#T_caeac_row13_col2{background-color:#d22b27;color:#f1f1f1}#T_caeac_row12_col0{background-color:#af0926;color:#f1f1f1}#T_caeac_row12_col2{background-color:#d42d27;color:#f1f1f1}#T_caeac_row12_col3{background-color:#e65036;color:#000}#T_caeac_row13_col1{background-color:#f16640;color:#000}#T_caeac_row14_col1,#T_caeac_row16_col2{background-color:#c62027;color:#f1f1f1}#T_caeac_row14_col2,#T_caeac_row14_col3{background-color:#f7814c;color:#000}#T_caeac_row15_col2{background-color:#fff3ac;color:#000}#T_caeac_row15_col3{background-color:#f36b42;color:#000}#T_caeac_row17_col0{background-color:#a70226;color:#f1f1f1}#T_caeac_row17_col2{background-color:#ca2427;color:#f1f1f1}</style><table
        id=T_caeac_><thead><tr><th class="level0 index_name">metric<th
        class="level0 col_heading col0"colspan=2>Bias (RMSD of all
        gridpoints)<th class="level0 col_heading col2"colspan=2>Mean change
        (Future - Reference)<tr><th class="level1 index_name">variable<th
        class="col0 col_heading level1">Temperature (K)<th class="col1
        col_heading level1">Precipitation (kg/m2/s)<th class="col2 col_heading
        level1">Temperature (K)<th class="col3 col_heading level1">Precipitation
        (kg/m2/s)<tr><th class="level0 index_name">dataset<th class=blank><th
        class=blank><th class=blank><th class=blank><tbody><tr><th class="level0
        row_heading row0"id=T_caeac_level0_row0>CMIP5_ACCESS1-0_r1i1p1<td
        class="data col0 row0"id=T_caeac_row0_col0>3.19e+00<td class="data col1
        row0"id=T_caeac_row0_col1>1.96e-05<td class="data col2
        row0"id=T_caeac_row0_col2>2.36e+00<td class="data col3
        row0"id=T_caeac_row0_col3>8.00e-09<tr><th class="level0 row_heading
        row1"id=T_caeac_level0_row1>CMIP5_BNU-ESM_r1i1p1<td class="data col0
        row1"id=T_caeac_row1_col0>4.08e+00<td class="data col1
        row1"id=T_caeac_row1_col1>1.87e-05<td class="data col2
        row1"id=T_caeac_row1_col2>2.44e+00<td class="data col3
        row1"id=T_caeac_row1_col3>2.96e-08<tr><th class="level0 row_heading
        row2"id=T_caeac_level0_row2>CMIP6_ACCESS-CM2_r1i1p1f1<td class="data
        col0 row2"id=T_caeac_row2_col0>3.75e+00<td class="data col1
        row2"id=T_caeac_row2_col1>1.77e-05<td class="data col2
        row2"id=T_caeac_row2_col2>2.87e+00<td class="data col3
        row2"id=T_caeac_row2_col3>6.63e-07<tr><th class="level0 row_heading
        row3"id=T_caeac_level0_row3>CMIP6_ACCESS-ESM1-5_r1i1p1f1<td class="data
        col0 row3"id=T_caeac_row3_col0>3.01e+00<td class="data col1
        row3"id=T_caeac_row3_col1>1.96e-05<td class="data col2
        row3"id=T_caeac_row3_col2>2.63e+00<td class="data col3
        row3"id=T_caeac_row3_col3>-1.39e-07<tr><th class="level0 row_heading
        row4"id=T_caeac_level0_row4>CMIP6_AWI-CM-1-1-MR_r1i1p1f1<td class="data
        col0 row4"id=T_caeac_row4_col0>2.91e+00<td class="data col1
        row4"id=T_caeac_row4_col1>1.80e-05<td class="data col2
        row4"id=T_caeac_row4_col2>2.56e+00<td class="data col3
        row4"id=T_caeac_row4_col3>7.67e-07<tr><th class="level0 row_heading
        row5"id=T_caeac_level0_row5>CMIP6_BCC-CSM2-MR_r1i1p1f1<td class="data
        col0 row5"id=T_caeac_row5_col0>4.22e+00<td class="data col1
        row5"id=T_caeac_row5_col1>1.74e-05<td class="data col2
        row5"id=T_caeac_row5_col2>2.64e+00<td class="data col3
        row5"id=T_caeac_row5_col3>5.02e-07<tr><th class="level0 row_heading
        row6"id=T_caeac_level0_row6>CMIP6_CAMS-CSM1-0_r1i1p1f1<td class="data
        col0 row6"id=T_caeac_row6_col0>4.43e+00<td class="data col1
        row6"id=T_caeac_row6_col1>1.84e-05<td class="data col2
        row6"id=T_caeac_row6_col2>1.48e+00<td class="data col3
        row6"id=T_caeac_row6_col3>4.89e-07<tr><th class="level0 row_heading
        row7"id=T_caeac_level0_row7>CMIP6_CESM2-WACCM_r1i1p1f1<td class="data
        col0 row7"id=T_caeac_row7_col0>2.95e+00<td class="data col1
        row7"id=T_caeac_row7_col1>1.69e-05<td class="data col2
        row7"id=T_caeac_row7_col2>2.33e+00<td class="data col3
        row7"id=T_caeac_row7_col3>-1.91e-07<tr><th class="level0 row_heading
        row8"id=T_caeac_level0_row8>CMIP6_CanESM5_r1i1p1f1<td class="data col0
        row8"id=T_caeac_row8_col0>2.81e+00<td class="data col1
        row8"id=T_caeac_row8_col1>1.69e-05<td class="data col2
        row8"id=T_caeac_row8_col2>3.36e+00<td class="data col3
        row8"id=T_caeac_row8_col3>2.10e-06<tr><th class="level0 row_heading
        row9"id=T_caeac_level0_row9>CMIP6_FGOALS-g3_r1i1p1f1<td class="data col0
        row9"id=T_caeac_row9_col0>6.74e+00<td class="data col1
        row9"id=T_caeac_row9_col1>1.80e-05<td class="data col2
        row9"id=T_caeac_row9_col2>2.13e+00<td class="data col3
        row9"id=T_caeac_row9_col3>5.95e-07<tr><th class="level0 row_heading
        row10"id=T_caeac_level0_row10>CMIP6_FIO-ESM-2-0_r1i1p1f1<td class="data
        col0 row10"id=T_caeac_row10_col0>3.02e+00<td class="data col1
        row10"id=T_caeac_row10_col1>1.75e-05<td class="data col2
        row10"id=T_caeac_row10_col2>2.07e+00<td class="data col3
        row10"id=T_caeac_row10_col3>1.89e-07<tr><th class="level0 row_heading
        row11"id=T_caeac_level0_row11>CMIP6_MIROC6_r1i1p1f1<td class="data col0
        row11"id=T_caeac_row11_col0>4.00e+00<td class="data col1
        row11"id=T_caeac_row11_col1>1.74e-05<td class="data col2
        row11"id=T_caeac_row11_col2>2.25e+00<td class="data col3
        row11"id=T_caeac_row11_col3>-2.45e-07<tr><th class="level0 row_heading
        row12"id=T_caeac_level0_row12>CMIP6_MPI-ESM1-2-HR_r1i1p1f1<td
        class="data col0 row12"id=T_caeac_row12_col0>2.98e+00<td class="data
        col1 row12"id=T_caeac_row12_col1>1.80e-05<td class="data col2
        row12"id=T_caeac_row12_col2>1.84e+00<td class="data col3
        row12"id=T_caeac_row12_col3>1.18e-07<tr><th class="level0 row_heading
        row13"id=T_caeac_level0_row13>CMIP6_MPI-ESM1-2-LR_r1i1p1f1<td
        class="data col0 row13"id=T_caeac_row13_col0>2.95e+00<td class="data
        col1 row13"id=T_caeac_row13_col1>1.78e-05<td class="data col2
        row13"id=T_caeac_row13_col2>1.82e+00<td class="data col3
        row13"id=T_caeac_row13_col3>2.52e-07<tr><th class="level0 row_heading
        row14"id=T_caeac_level0_row14>CMIP6_MRI-ESM2-0_r1i1p1f1<td class="data
        col0 row14"id=T_caeac_row14_col0>2.81e+00<td class="data col1
        row14"id=T_caeac_row14_col1>1.71e-05<td class="data col2
        row14"id=T_caeac_row14_col2>2.36e+00<td class="data col3
        row14"id=T_caeac_row14_col3>5.75e-07<tr><th class="level0 row_heading
        row15"id=T_caeac_level0_row15>CMIP6_NESM3_r1i1p1f1<td class="data col0
        row15"id=T_caeac_row15_col0>3.90e+00<td class="data col1
        row15"id=T_caeac_row15_col1>1.83e-05<td class="data col2
        row15"id=T_caeac_row15_col2>3.22e+00<td class="data col3
        row15"id=T_caeac_row15_col3>3.60e-07<tr><th class="level0 row_heading
        row16"id=T_caeac_level0_row16>CMIP6_NorESM2-LM_r1i1p1f1<td class="data
        col0 row16"id=T_caeac_row16_col0>3.08e+00<td class="data col1
        row16"id=T_caeac_row16_col1>1.70e-05<td class="data col2
        row16"id=T_caeac_row16_col2>1.74e+00<td class="data col3
        row16"id=T_caeac_row16_col3>-4.97e-07<tr><th class="level0 row_heading
        row17"id=T_caeac_level0_row17>CMIP6_NorESM2-MM_r1i1p1f1<td class="data
        col0 row17"id=T_caeac_row17_col0>2.86e+00<td class="data col1
        row17"id=T_caeac_row17_col1>1.67e-05<td class="data col2
        row17"id=T_caeac_row17_col2>1.76e+00<td class="data col3
        row17"id=T_caeac_row17_col3>-7.65e-07</table>

    </embed>
