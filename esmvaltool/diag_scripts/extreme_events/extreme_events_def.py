"""
Created on Tue May 19 17:23:06 2020

@author: bmueller
"""
import yaml

index_definition = yaml.safe_load("""
# %% annual_number_of_frost_days                     
annual_number_of_frost_days:
    name: frost days
    period: annual
    required:
        - tasmin
    threshold:
        value: 273.15
        unit: K
        logic: less
    cf_name: annual_number_of_frost_days
# %% annual_number_of_summer_days
annual_number_of_summer_days:
    name: summer days
    period: annual
    required:
        - tasmax
    threshold:
        value: 298.15
        unit: K
        logic: greater
    cf_name: number_of_days_with_air_temperature_above_25_degree_Celsius
 # %% annual_number_of_icing_days
annual_number_of_icing_days:
    name: icing days
    period: annual
    required:
        - tasmax
    threshold:
        value: 273.15
        unit: K
        logic: less
    cf_name: number_of_days_where_air_temperature_remains_below_freezing_point
# %% annual_number_of_tropical_nights
annual_number_of_tropical_nights:
    name: tropical nights
    period: annual
    required:
        - tasmin
    threshold:
        value: 293.15
        unit: K
        logic: greater
    cf_name: number_of_days_where_air_temperature_remains_above_20_degre_Celsius
# %% annual_number_of_days_where_cumulative_precipitation_is_above_10_mm
annual_number_of_days_where_cumulative_precipitation_is_above_10_mm:
    name: R10mm
    period: annual
    required:
        - pr
    threshold:
        value: 10
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_10_mm
# %% annual_number_of_days_where_cumulative_precipitation_is_above_20_mm
annual_number_of_days_where_cumulative_precipitation_is_above_20_mm:
    name: R20mm
    period: annual
    required:
        - pr
    threshold:
        value: 20
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_20_mm
# %% annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm
annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm:
    name: R{}mm
    period: annual
    required:
        - pr
    threshold:
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_{}_mm
# %% monthly_maximum_value_of_daily_maximum_temperature
monthly_maximum_value_of_daily_maximum_temperature:
    name: TXx
    period: monthly
    required:
        - tasmax
    logic: max
    cf_name: monthly_maximum_value_of_daily_maximum_temperature
  # %% monthly_maximum_value_of_daily_minimum_temperature
monthly_maximum_value_of_daily_minimum_temperature:
    name: TNx
    period: monthly
    required:
        - tasmin
    logic: max
    cf_name: monthly_maximum_value_of_daily_minimum_temperature
# %% monthly_minimum_value_of_daily_maximum_temperature
monthly_minimum_value_of_daily_maximum_temperature:
    name: TXn
    period: monthly
    required:
        - tasmax
    logic: min
    cf_name: monthly_minimum_value_of_daily_maximum_temperature
# %% monthly_minimum_value_of_daily_minimum_temperature
monthly_minimum_value_of_daily_minimum_temperature:
    name: TNn
    period: monthly
    required:
        - tasmin
    logic: min
    cf_name: monthly_minimum_value_of_daily_minimum_temperature
# %% monthly_maximum_1day_precipitation
monthly_maximum_1day_precipitation:
    name: Rx1day
    period: monthly
    required:
        - pr
    spell:
        value: 1
        unit: day
    logic: max
    cf_name: monthly_maximum_1day_precipitation
# %% monthly_maximum_5day_precipitation
monthly_maximum_5day_precipitation:
    name: Rx5day
    period: monthly
    required:
        - pr
    spell:
        value: 5
        unit: day
    logic: max
    cf_name: monthly_maximum_5day_precipitation
# %% annual_total_precipitation_in_wet_days
annual_total_precipitation_in_wet_days:
    name: PRCPTOT
    period: annual
    required:
        - pr
    logic: sum
    cf_name: annual_total_precipitation_in_wet_days
# %% annual_total_precipitation_in_wet_days_where_daily_precipitation_above_95%
annual_total_precipitation_in_wet_days_where_daily_precipitation_above_95%:
    name: R95pTOT
    period: annual
    required:
        - pr
    threshold:
        value: 95
        unit: percent
        logic: greater
    logic: sum
    cf_name: annual_total_precipitation_in_wet_days_where_daily_precipitation_above_95%
# %% annual_total_precipitation_in_wet_days_where_daily_precipitation_above_99%
annual_total_precipitation_in_wet_days_where_daily_precipitation_above_99%:
    name: R99pTOT
    period: annual
    required:
        - pr
    threshold:
        value: 99
        unit: percent
        logic: greater
    logic: sum
    cf_name: annual_total_precipitation_in_wet_days_where_daily_precipitation_above_99%
# %% monthly_number_of_days_where_daily_minimum_temperature_below_10%
monthly_number_of_days_where_daily_minimum_temperature_below_10%:
    name: TN10p
    period: monthly
    required:
        - tasmin
    threshold:
        value: 10
        unit: percent
        logic: less
    cf_name: monthly_number_of_days_where_daily_minimum_temperature_below_10%
# %% monthly_number_of_days_where_daily_minimum_temperature_above_90%
monthly_number_of_days_where_daily_minimum_temperature_above_90%:
    name: TN90p
    period: monthly
    required:
        - tasmin
    threshold:
        value: 90
        unit: percent
        logic: greater
    cf_name: monthly_number_of_days_where_daily_minimum_temperature_above_90%
# %% monthly_number_of_days_where_daily_maximum_temperature_below_10%
monthly_number_of_days_where_daily_maximum_temperature_below_10%:
    name: TX10p
    period: monthly
    required:
        - tasmax
    threshold:
        value: 10
        unit: percent
        logic: less
    cf_name: monthly_number_of_days_where_daily_maximum_temperature_below_10%
# %% monthly_number_of_days_where_daily_maximum_temperature_above_90%
monthly_number_of_days_where_daily_maximum_temperature_above_90%:
    name: TX90p
    period: monthly
    required:
        - tasmax
    threshold:
        value: 90
        unit: percent
        logic: greater
    cf_name: monthly_number_of_days_where_daily_maximum_temperature_above_90%
# %% annual_warm_spell_duration_index
annual_warm_spell_duration_index:
    name: WSDI
    period: annual
    required:
        - tasmax
    threshold:
        value: 90
        unit: percent
        logic: greater
    spell:
        value: 6
        unit: days
        logic: greater_equal
    cf_name: annual_warm_spell_duration_index 
# %% annual_cold_spell_duration_index
annual_cold_spell_duration_index:
    name: CSDI
    period: annual
    required:
        - tasmin
    threshold:
        value: 10
        unit: percent
        logic: less
    spell:
        value: 6
        unit: days
        logic: greater_equal
    cf_name: annual_cold_spell_duration_index
# %% annual_maximum_length_of_dry_spell
annual_maximum_length_of_dry_spell:
    name: CDD
    period: annual
    required:
        - pr
    threshold:
        value: 1
        unit: mm day-1
        logic: less
    cf_name: annual_maximum_length_of_dry_spell
# %% annual_maximum_length_of_wet_spell
annual_maximum_length_of_wet_spell:
    name: CWD
    period: annual
    required:
        - pr
    threshold:
        value: 1
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_maximum_length_of_wet_spell
# %% monthly_simple_precipitation_intensity_index
monthly_simple_precipitation_intensity_index:
    name: SDII
    period: monthly
    required:
        - pr
    threshold:
        value: 1
        unit: mm day-1
        logic: greater_equal
    cf_name: monthly_simple_precipitation_intensity_index
# %% daily_temperature_range
daily_temperature_range:
    name: DTR
    period: daily
    required:
        - tasmin
        - tasmax
    logic: diff
    cf_name: daily_temperature_range
# %% annual_growing_season_length
annual_growing_season_length:
    name: GSL
    required:
        - tas
    start:
        threshold:
            value: 278.15
            unit: K
            logic: greater
        spell:
            value: 6
            unit: days
            logic: equal
        time:
            delay: 0
            len: 6
            unit: month
    end:
        threshold:
            value: 278.15
            unit: K
            logic: less
        spell:
            value: 6
            unit: days
            logic: equal
        time:
            delay: 6
            len: 6
            unit: month
    spatial_subsets:
        NH:
            latitude: [0, 90]
            longitude: [0, 360]
        SH:
            latitude: [0, -90]
            longitude: [0, 360]
    cf_name: annual_growing_season_length
""")
print("INDEX_DEFINITION:")
print(yaml.dump(index_definition))
