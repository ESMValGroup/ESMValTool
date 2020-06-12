.. _recipe_kcs:

KNMI Climate Scenarios
======================

Description
-----------

This recipe implements the method described in `Lenderink et al., 2014 <https://doi.org/10.1088/1748-9326/9/11/115008>`_, to prepare the 2014 KNMI Climate Scenarios (KCS) for the Netherlands. A set of 8 global climate projections from EC-Earth were downscaled with the RACMO regional climate model. Since the EC-Earth ensemble is not readily representative for the spread in the full CMIP ensemble, this method recombines 5-year segments from the EC-Earth ensemble to obtain a large suite of "resamples". Subsequently, 8 new resamples are selected that cover the spread in CMIP much better than the original set.

The original method created 16 resampled datasets:

* 2 main scenarios: Moderate (M) and Warm (W) (Lenderink 2014 uses "G" instead of "M").
* 2 'sub'scenarios: Relatively high (H) or low (L) changes in seasonal temperature and precipitation
* 2 time horizons: Mid-century (MOC; 2050) and end-of-century (EOC; 2085)
* 2 periods: Control (1981-2010) and future (variable)

The configuration settings for these scenarios can be found in table 1 of Lenderink 2014's supplement.

Implementation
--------------

The implementation is such that application to other datasets, regions, etc. is relatively straightforward. However, the description below focuses on the reference use case of Lenderink et al., 2014, where the target model is EC-Earth.

In the first diagnostic, the spread of the full CMIP ensemble is used to obtain 4 values of a *global* :math:`{\Delta}T_{CMIP}`, corresponding to the 10th and 90th percentiles for the M and W scenarios, respectively, for both MOC and EOC. Subsequently, for each of these 4 *steering parameters*, 30-year periods are selected from the EC-Earth ensemble, where :math:`{\Delta}T_{ECEarth}{\approx}{\Delta}T_{CMIP}`.

In the second diagnostic, for for both the control and future periods, the 8 EC-Earth ensemble members are split into 6 segments of 5 years each. Out of all :math:`8^6` possible re-combinations of these 5-year segments, eventually 8 new 'resamples' are selected based on *local* changes in seasonal temperature and precipitation. This is done in the following steps:

1. Select 1000 samples for the control period, and 2 x 1000 samples for the future period (one for each subscenario). Step 1 poses a constraint on winter precipitation. For the control period, winter precipitation must still closely represent the average of the original ensemble. For the two future periods, the change in winter precipitation with respect to the control period must approximately equal :math:`4{\Delta}T` % (subscenario L) or  :math:`8{\Delta}T` % (subscenario H).
2. Further constrain the selection by picking samples that represent either high or low changes in summer precipitation and summer and winter temperature, by limiting the remaining samples to certain percentile ranges: relatively wet/cold in the control and dry/warm in the future, or vice versa. The percentile ranges are listed in table 1 of Lenderink 2014's supplement. This should result is approximately 50 remaining samples for each scenario, for both control and future.
3. Use a monte-carlo method to make a final selection of 8 resamples with minimal reuse of the same ensemble member/segment.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

- recipe_kcs.yml

Diagnostics are stored in diag_scripts/kcs/

- global_matching.py
- local_resampling.py


User settings
-------------

Datasets: Datasets have been split in two parts, such that the CMIP ensemble and the target model datasets can be treated separately. The recipe can work with a target model that is not part of CMIP. In that case the target model data should be stored alongside the CMIP data in the exact same format. One use case for this recipe is to compare between CMIP5 and CMIP6.

Preprocessors: We've tried to use built-in preprocessor functions as much as possible. The first diagnostic requires global mean temperature anomalies for each dataset, both CMIP and the target model, and some multimodel statistics. The second diagnostic requires focuses on a point in the Netherlands. The :code:`extract_point` preprocessor can be changed to :code:`extract_shape` or :code:`extract_region`, in conjunction with an area mean. And of course, the coordinates can be changed to analyze a different region.

Diagnostics:

* global_matching

  * :code:`scenario_years`: a list of time horizons. Default: :code:`[2050, 2085]`
  * :code:`scenario_percentiles`: a list of percentiles for the steering table. Default: :code:`[p10, p90]`
  * :code:`target_model`: the name of the target model. Default: :code:`EC-Earth`

* local_resampling

  * :code:`target_model`: the name of the target model. Default: :code:`EC-Earth`
  * :code:`scenarios`: a scenario name and list of options. The default setting is a single scenario:

    .. code-block:: yaml

        scenarios:
          ML:
            description: "Moderate / low changes in seasonal temperature & precipitation"
            global_dT: 1.0
            scenario_year: 2050
            resampling_period: [2021, 2050]
            dpr_winter: 4
            pr_summer_control: [25, 55]
            pr_summer_future: [45, 75]
            tas_winter_control: [50, 80]
            tas_winter_future: [20, 50]
            tas_summer_control: [0, 100]
            tas_summer_future: [0, 50]

    These values are taken from table 1 in the Lenderink 2014's supplement. For new applications, :code:`global_dT`, :code:`resampling_period` and :code:`dpr_winter` are informed by the output of the first diagnostic. The percentile bounds are to be tuned until a satisfactory selection is achieved. Multiple scenarios can be processed at once by appending more configurations below the default one.

Example output
--------------

The diagnostic :code:`global matching` produces a scenarios table like the one below

.. code-block:: python

         year percentile   cmip_dt period_bounds           target_dt  pattern_scaling_factor
      0  2050       Mean  1.494123  [2055, 2084]   1.496281513828891                0.998557
      1  2050     Median  1.434082  [2051, 2080]  1.4258415170339225                1.005779
      2  2085       Mean  1.916866  [2085, 2114]   1.877538437763913                1.020946
      3  2085     Median  1.383674  [2049, 2078]  1.3914389868132169                0.994419

which is printed to the log file and also saved as a csv-file :code:`scenarios.csv`.
Additionally, a figure is created showing the CMIP spread in global temperature change,
AND highlighting the selected steering parameters and resampling periods:

.. _fig_kcs_global_matching:
.. figure::  /recipes/figures/kcs/global_matching.png
   :align:   center


.. TODO Add second diagnostic output



References
----------

* `Lenderink et al. 2014, Environ. Res. Lett., 9, 115008 <https://doi.org/10.1088/1748-9326/9/11/115008>`_.
