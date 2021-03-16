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

    * impact/quick_insights.py: tabulate and visualize bias and change.


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

   "Bias and change for each variable"

.. raw:: html
    <embed>

        <style type="text/css" > #T_22a8e_row0_col0{ background-color: #e24731; color:
        #f1f1f1;
        }#T_22a8e_row0_col1,#T_22a8e_row8_col3,#T_22a8e_row9_col0,#T_22a8e_row15_col2{
        background-color: #eef8a8; color: #000000; }#T_22a8e_row0_col2{
        background-color: #fcaa5f; color: #000000;
        }#T_22a8e_row0_col3,#T_22a8e_row12_col2{ background-color: #f36b42; color:
        #000000; }#T_22a8e_row1_col0{ background-color: #f7844e; color: #000000;
        }#T_22a8e_row1_col1{ background-color: #fee28f; color: #000000;
        }#T_22a8e_row1_col2{ background-color: #fdb567; color: #000000;
        }#T_22a8e_row1_col3,#T_22a8e_row17_col2{ background-color: #f16640; color:
        #000000; }#T_22a8e_row2_col0,#T_22a8e_row16_col2{ background-color: #f47044;
        color: #000000; }#T_22a8e_row2_col1,#T_22a8e_row9_col2{ background-color:
        #fa9656; color: #000000; }#T_22a8e_row2_col2{ background-color: #feeb9d; color:
        #000000; }#T_22a8e_row2_col3,#T_22a8e_row4_col1,#T_22a8e_row12_col1{
        background-color: #fdb365; color: #000000;
        }#T_22a8e_row3_col0,#T_22a8e_row10_col0{ background-color: #db382b; color:
        #f1f1f1; }#T_22a8e_row3_col1{ background-color: #eff8aa; color: #000000;
        }#T_22a8e_row3_col2,#T_22a8e_row5_col2{ background-color: #fed884; color:
        #000000; }#T_22a8e_row3_col3{ background-color: #ef633f; color: #000000;
        }#T_22a8e_row4_col0,#T_22a8e_row17_col0{ background-color: #d83128; color:
        #f1f1f1; }#T_22a8e_row4_col2{ background-color: #fed481; color: #000000;
        }#T_22a8e_row4_col3{ background-color: #fdbb6c; color: #000000;
        }#T_22a8e_row5_col0{ background-color: #fa9b58; color: #000000;
        }#T_22a8e_row5_col1,#T_22a8e_row10_col1,#T_22a8e_row15_col0{ background-color:
        #f67f4b; color: #000000;
        }#T_22a8e_row5_col3,#T_22a8e_row6_col0,#T_22a8e_row15_col3{ background-color:
        #fca55d; color: #000000; }#T_22a8e_row6_col1{ background-color: #fece7c; color:
        #000000;
        }#T_22a8e_row6_col2,#T_22a8e_row8_col0,#T_22a8e_row14_col0,#T_22a8e_row17_col1,#T_22a8e_row17_col3{
        background-color: #d22b27; color: #f1f1f1;
        }#T_22a8e_row6_col3,#T_22a8e_row13_col1{ background-color: #fba05b; color:
        #000000; }#T_22a8e_row7_col0,#T_22a8e_row12_col0,#T_22a8e_row13_col0{
        background-color: #d93429; color: #f1f1f1; }#T_22a8e_row7_col1{
        background-color: #de402e; color: #f1f1f1; }#T_22a8e_row7_col2{
        background-color: #fdb768; color: #000000;
        }#T_22a8e_row7_col3,#T_22a8e_row13_col2{ background-color: #ed5f3c; color:
        #000000; }#T_22a8e_row8_col1,#T_22a8e_row16_col0,#T_22a8e_row16_col3{
        background-color: #dd3d2d; color: #f1f1f1; }#T_22a8e_row8_col2{
        background-color: #f7fcb4; color: #000000; }#T_22a8e_row9_col1{
        background-color: #fdaf62; color: #000000;
        }#T_22a8e_row9_col3,#T_22a8e_row14_col3{ background-color: #fdb163; color:
        #000000; }#T_22a8e_row10_col2,#T_22a8e_row12_col3{ background-color: #f88950;
        color: #000000; }#T_22a8e_row10_col3{ background-color: #f8864f; color:
        #000000; }#T_22a8e_row11_col0{ background-color: #f7814c; color: #000000;
        }#T_22a8e_row11_col1{ background-color: #f26841; color: #000000;
        }#T_22a8e_row11_col2{ background-color: #fdad60; color: #000000;
        }#T_22a8e_row11_col3{ background-color: #ee613e; color: #000000;
        }#T_22a8e_row13_col3{ background-color: #fa9857; color: #000000;
        }#T_22a8e_row14_col1{ background-color: #e95538; color: #000000;
        }#T_22a8e_row14_col2,#T_22a8e_row15_col1{ background-color: #fdbd6d; color:
        #000000; }#T_22a8e_row16_col1{ background-color: #e54e35; color: #000000;
        }</style><table id="T_22a8e_" ><thead> <tr> <th class="index_name level0"
        >metric</th> <th class="col_heading level0 col0" colspan="2">bias</th> <th
        class="col_heading level0 col2" colspan="2">change</th> </tr> <tr> <th
        class="index_name level1" >variable</th> <th class="col_heading level1 col0"
        >tas</th> <th class="col_heading level1 col1" >pr</th> <th class="col_heading
        level1 col2" >tas</th> <th class="col_heading level1 col3" >pr</th> </tr> <tr>
        <th class="index_name level0" >dataset</th> <th class="blank" ></th> <th
        class="blank" ></th> <th class="blank" ></th> <th class="blank" ></th>
        </tr></thead><tbody> <tr> <th id="T_22a8e_level0_row0" class="row_heading
        level0 row0" >CMIP5_ACCESS1-0</th> <td id="T_22a8e_row0_col0" class="data row0
        col0" >3.21e+00</td> <td id="T_22a8e_row0_col1" class="data row0 col1"
        >1.95e-05</td> <td id="T_22a8e_row0_col2" class="data row0 col2" >2.46e+00</td>
        <td id="T_22a8e_row0_col3" class="data row0 col3" >9.21e-09</td> </tr> <tr> <th
        id="T_22a8e_level0_row1" class="row_heading level0 row1" >CMIP5_BNU-ESM</th>
        <td id="T_22a8e_row1_col0" class="data row1 col0" >4.04e+00</td> <td
        id="T_22a8e_row1_col1" class="data row1 col1" >1.86e-05</td> <td
        id="T_22a8e_row1_col2" class="data row1 col2" >2.55e+00</td> <td
        id="T_22a8e_row1_col3" class="data row1 col3" >-4.89e-08</td> </tr> <tr> <th
        id="T_22a8e_level0_row2" class="row_heading level0 row2" >CMIP6_ACCESS-CM2</th>
        <td id="T_22a8e_row2_col0" class="data row2 col0" >3.76e+00</td> <td
        id="T_22a8e_row2_col1" class="data row2 col1" >1.78e-05</td> <td
        id="T_22a8e_row2_col2" class="data row2 col2" >3.12e+00</td> <td
        id="T_22a8e_row2_col3" class="data row2 col3" >7.23e-07</td> </tr> <tr> <th
        id="T_22a8e_level0_row3" class="row_heading level0 row3"
        >CMIP6_ACCESS-ESM1-5</th> <td id="T_22a8e_row3_col0" class="data row3 col0"
        >3.00e+00</td> <td id="T_22a8e_row3_col1" class="data row3 col1" >1.94e-05</td>
        <td id="T_22a8e_row3_col2" class="data row3 col2" >2.88e+00</td> <td
        id="T_22a8e_row3_col3" class="data row3 col3" >-5.78e-08</td> </tr> <tr> <th
        id="T_22a8e_level0_row4" class="row_heading level0 row4"
        >CMIP6_AWI-CM-1-1-MR</th> <td id="T_22a8e_row4_col0" class="data row4 col0"
        >2.92e+00</td> <td id="T_22a8e_row4_col1" class="data row4 col1" >1.80e-05</td>
        <td id="T_22a8e_row4_col2" class="data row4 col2" >2.85e+00</td> <td
        id="T_22a8e_row4_col3" class="data row4 col3" >8.19e-07</td> </tr> <tr> <th
        id="T_22a8e_level0_row5" class="row_heading level0 row5"
        >CMIP6_BCC-CSM2-MR</th> <td id="T_22a8e_row5_col0" class="data row5 col0"
        >4.33e+00</td> <td id="T_22a8e_row5_col1" class="data row5 col1" >1.76e-05</td>
        <td id="T_22a8e_row5_col2" class="data row5 col2" >2.88e+00</td> <td
        id="T_22a8e_row5_col3" class="data row5 col3" >5.69e-07</td> </tr> <tr> <th
        id="T_22a8e_level0_row6" class="row_heading level0 row6"
        >CMIP6_CAMS-CSM1-0</th> <td id="T_22a8e_row6_col0" class="data row6 col0"
        >4.46e+00</td> <td id="T_22a8e_row6_col1" class="data row6 col1" >1.84e-05</td>
        <td id="T_22a8e_row6_col2" class="data row6 col2" >1.50e+00</td> <td
        id="T_22a8e_row6_col3" class="data row6 col3" >5.30e-07</td> </tr> <tr> <th
        id="T_22a8e_level0_row7" class="row_heading level0 row7"
        >CMIP6_CESM2-WACCM</th> <td id="T_22a8e_row7_col0" class="data row7 col0"
        >2.97e+00</td> <td id="T_22a8e_row7_col1" class="data row7 col1" >1.70e-05</td>
        <td id="T_22a8e_row7_col2" class="data row7 col2" >2.57e+00</td> <td
        id="T_22a8e_row7_col3" class="data row7 col3" >-1.28e-07</td> </tr> <tr> <th
        id="T_22a8e_level0_row8" class="row_heading level0 row8" >CMIP6_CanESM5</th>
        <td id="T_22a8e_row8_col0" class="data row8 col0" >2.82e+00</td> <td
        id="T_22a8e_row8_col1" class="data row8 col1" >1.69e-05</td> <td
        id="T_22a8e_row8_col2" class="data row8 col2" >3.54e+00</td> <td
        id="T_22a8e_row8_col3" class="data row8 col3" >2.20e-06</td> </tr> <tr> <th
        id="T_22a8e_level0_row9" class="row_heading level0 row9" >CMIP6_FGOALS-g3</th>
        <td id="T_22a8e_row9_col0" class="data row9 col0" >6.64e+00</td> <td
        id="T_22a8e_row9_col1" class="data row9 col1" >1.80e-05</td> <td
        id="T_22a8e_row9_col2" class="data row9 col2" >2.31e+00</td> <td
        id="T_22a8e_row9_col3" class="data row9 col3" >6.84e-07</td> </tr> <tr> <th
        id="T_22a8e_level0_row10" class="row_heading level0 row10"
        >CMIP6_FIO-ESM-2-0</th> <td id="T_22a8e_row10_col0" class="data row10 col0"
        >3.02e+00</td> <td id="T_22a8e_row10_col1" class="data row10 col1"
        >1.76e-05</td> <td id="T_22a8e_row10_col2" class="data row10 col2"
        >2.22e+00</td> <td id="T_22a8e_row10_col3" class="data row10 col3"
        >2.84e-07</td> </tr> <tr> <th id="T_22a8e_level0_row11" class="row_heading
        level0 row11" >CMIP6_MIROC6</th> <td id="T_22a8e_row11_col0" class="data row11
        col0" >4.02e+00</td> <td id="T_22a8e_row11_col1" class="data row11 col1"
        >1.74e-05</td> <td id="T_22a8e_row11_col2" class="data row11 col2"
        >2.47e+00</td> <td id="T_22a8e_row11_col3" class="data row11 col3"
        >-8.52e-08</td> </tr> <tr> <th id="T_22a8e_level0_row12" class="row_heading
        level0 row12" >CMIP6_MPI-ESM1-2-HR</th> <td id="T_22a8e_row12_col0" class="data
        row12 col0" >2.96e+00</td> <td id="T_22a8e_row12_col1" class="data row12 col1"
        >1.81e-05</td> <td id="T_22a8e_row12_col2" class="data row12 col2"
        >2.01e+00</td> <td id="T_22a8e_row12_col3" class="data row12 col3"
        >3.13e-07</td> </tr> <tr> <th id="T_22a8e_level0_row13" class="row_heading
        level0 row13" >CMIP6_MPI-ESM1-2-LR</th> <td id="T_22a8e_row13_col0" class="data
        row13 col0" >2.94e+00</td> <td id="T_22a8e_row13_col1" class="data row13 col1"
        >1.79e-05</td> <td id="T_22a8e_row13_col2" class="data row13 col2"
        >1.91e+00</td> <td id="T_22a8e_row13_col3" class="data row13 col3"
        >4.48e-07</td> </tr> <tr> <th id="T_22a8e_level0_row14" class="row_heading
        level0 row14" >CMIP6_MRI-ESM2-0</th> <td id="T_22a8e_row14_col0" class="data
        row14 col0" >2.82e+00</td> <td id="T_22a8e_row14_col1" class="data row14 col1"
        >1.72e-05</td> <td id="T_22a8e_row14_col2" class="data row14 col2"
        >2.62e+00</td> <td id="T_22a8e_row14_col3" class="data row14 col3"
        >6.93e-07</td> </tr> <tr> <th id="T_22a8e_level0_row15" class="row_heading
        level0 row15" >CMIP6_NESM3</th> <td id="T_22a8e_row15_col0" class="data row15
        col0" >3.98e+00</td> <td id="T_22a8e_row15_col1" class="data row15 col1"
        >1.82e-05</td> <td id="T_22a8e_row15_col2" class="data row15 col2"
        >3.64e+00</td> <td id="T_22a8e_row15_col3" class="data row15 col3"
        >5.64e-07</td> </tr> <tr> <th id="T_22a8e_level0_row16" class="row_heading
        level0 row16" >CMIP6_NorESM2-LM</th> <td id="T_22a8e_row16_col0" class="data
        row16 col0" >3.10e+00</td> <td id="T_22a8e_row16_col1" class="data row16 col1"
        >1.71e-05</td> <td id="T_22a8e_row16_col2" class="data row16 col2"
        >2.03e+00</td> <td id="T_22a8e_row16_col3" class="data row16 col3"
        >-4.65e-07</td> </tr> <tr> <th id="T_22a8e_level0_row17" class="row_heading
        level0 row17" >CMIP6_NorESM2-MM</th> <td id="T_22a8e_row17_col0" class="data
        row17 col0" >2.92e+00</td> <td id="T_22a8e_row17_col1" class="data row17 col1"
        >1.68e-05</td> <td id="T_22a8e_row17_col2" class="data row17 col2"
        >1.97e+00</td> <td id="T_22a8e_row17_col3" class="data row17 col3"
        >-6.64e-07</td> </tr> </tbody></table>

    </embed>
