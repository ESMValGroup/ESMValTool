.. _recipe_kcs:

KNMI Climate Scenarios
======================

Overview
--------

This recipe implements the method described in Lenderink et al., 2014, to prepare the 2014 KNMI Climate Scenarios (KCS) for the Netherlands.
For these scenarios, a set of 8 global climate projections from EC-Earth were downscaled with the RACMO regional climate model. Since the EC-Earth ensemble is not readily representative for the spread in the full CMIP ensemble, this method recombines 5-year segments from the EC-Earth ensemble to obtain a large suite of recombined climates. From these new climates, 8 combinations are selected that together are much more representative for the CMIP spread than the original set.

The implementation is such, that application to other datasets, regions, etc. is relatively straightforward. However, the description below focuses on the reference use case of Lenderink et al., 2014.

The original method creates 16 resampled datasets:
- 2 main scenarios: Moderate (M) and Warm (W) (in Lenderink's paper is "G" is used instead of "M", for "gematigd" in Dutch).
- 2 'sub'scenarios: Relatively high (H) or low (L) changes in seasonal temperature and precipitation
- 2 time horizons: Mid-century (MOC; 2050) and end-of-century (EOC; 2085)
- 2 periods: Control (1981-2010) and future (variable)

The configuration settings for each of these scenarios can be found in table 1 of Lenderink 2014's supplement.

In the first diagnostic, the spread of the full CMIP ensemble is used to get a 'CMIP delta T' corresponding to the 10th and 90th percentile for the moderate and warm scenarios. Subsequently, 30-year periods are selected where the delta_T in the target model (EC-Earth) roughly matches the 'CMIP delta T' for each scenario (W/M, MOC/EOC).

In the second diagnostic, for for both the control and future periods, the 8 ensemble members of the target model are segmented 6 segments of 5 years each. Out of all 8^6 possible re-combinations of these 5-year segments, eventually 8 new 're-samples' are selected based on changes in seasonal temperature and precipitation. This is done in the following steps:

1. Select 1000 samples for the control period, and 2 x 1000 samples for the future period. The two future periods will be used to  obtaining (sub)scenarios with either relatively low or high seasonal changes in temperature and precipitation. For now, only pose a constraint on winter precipitation. For the control period, winter precipitation must still closely represent the average winter precipitation of the 8 original members. For the two future periods, the change in winter precipitation with respect to the control period must approximately equal 4 x Delta T % or  8 x Delta T %, respectively.
2. Further constrain the selection by picking samples that represent either high or low change in summer precipitation and summer and winter temperature. This is done by limiting the remaining samples for both the control and the future periods to a certain percentile range: relatively wet/cold in the control and dry/warm in the future, or vice versa. The percentile ranges are listed in table 1 of Lenderink 2014's supplement. This should result is approximately 50 remaining samples for each scenario, for both control and future.
3. Out of the remaining samples, a final selection of 8 recombinations is made that has a minimal reuse of the same ensemble member/segment. This is done using a monte-carlo method.


References
----------

Lenderink et al., 2014: https://doi.org/10.1088/1748-9326/9/11/115008


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

- recipe_kcs.yml

Diagnostics are stored in diag_scripts/kcs/

- scenarios_table.py
- resampling.py

User settings
-------------

Datasets: Datasets have been split in two parts. One part is for the entire CMIP5 ensemble. Datasets for both historical rcp's have been combined manually, so that they are concatenated before all the preprocessing steps are performed. Sometimes, the same historical dataset is used multiple times. This should not be a problem as long as the multi-model statistics for the historical period are interpreted with caution. The other part is for the target model. The recipe can work with a target model that is not part of CMIP. In that case the target model data should be stored alongside the CMIP data in the exact same format.


Preprocessors:


Scripts:


Example output
--------------

