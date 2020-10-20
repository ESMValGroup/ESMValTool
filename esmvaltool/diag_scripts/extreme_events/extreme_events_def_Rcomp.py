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
        value: 0.0
        unit: celsius
        logic: less
    cf_name: fdETCCDI
# %% annual_number_of_summer_days
annual_number_of_summer_days:
    name: summer days
    period: annual
    required:
        - tasmax
    threshold:
        value: 25.0
        unit: celsius
        logic: greater
    cf_name: suETCCDI
 # %% annual_number_of_icing_days
annual_number_of_icing_days:
    name: icing days
    period: annual
    required:
        - tasmax
    threshold:
        value: 0.0
        unit: celsius
        logic: less
    cf_name: idETCCDI
# %% annual_number_of_tropical_nights
annual_number_of_tropical_nights:
    name: tropical nights
    period: annual
    required:
        - tasmin
    threshold:
        value: 20.0 
        unit: celsius
        logic: greater
    cf_name: trETCCDI
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
    cf_name: r10mmETCCDI
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
    cf_name: r20mmETCCDI
# %% annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm
annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm:
    name: R{}mm
    period: annual
    required:
        - pr
    threshold:
        unit: mm day-1
        logic: greater_equal
    cf_name: r{}mmETCCDI
# %% monthly_maximum_value_of_daily_maximum_temperature
monthly_maximum_value_of_daily_maximum_temperature:
    name: TXx
    period: monthly
    required:
        - tasmax
    logic: max
    cf_name: txxETCCDI
# %% annual_maximum_value_of_daily_maximum_temperature
annual_maximum_value_of_daily_maximum_temperature:
    name: TXx
    period: annual
    required:
        - tasmax
    logic: max
    cf_name: txxETCCDI
# %% monthly_maximum_value_of_daily_minimum_temperature
monthly_maximum_value_of_daily_minimum_temperature:
    name: TNx
    period: monthly
    required:
        - tasmin
    logic: max
    cf_name: tnxETCCDI
# %% annual_maximum_value_of_daily_minimum_temperature
annual_maximum_value_of_daily_minimum_temperature:
    name: TNx
    period: annual
    required:
        - tasmin
    logic: max
    cf_name: tnxETCCDI
# %% monthly_minimum_value_of_daily_maximum_temperature
monthly_minimum_value_of_daily_maximum_temperature:
    name: TXn
    period: monthly
    required:
        - tasmax
    logic: min
    cf_name: txnETCCDI
# %% annual_minimum_value_of_daily_maximum_temperature
annual_minimum_value_of_daily_maximum_temperature:
    name: TXn
    period: annual
    required:
        - tasmax
    logic: min
    cf_name: txnETCCDI
# %% monthly_minimum_value_of_daily_minimum_temperature
monthly_minimum_value_of_daily_minimum_temperature:
    name: TNn
    period: monthly
    required:
        - tasmin
    logic: min
    cf_name: tnnETCCDI
# %% annual_minimum_value_of_daily_minimum_temperature
annual_minimum_value_of_daily_minimum_temperature:
    name: TNn
    period: annual
    required:
        - tasmin
    logic: min
    cf_name: tnnETCCDI
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
    cf_name: rx1dayETCCDI
# %% annual_maximum_1day_precipitation
annual_maximum_1day_precipitation:
    name: Rx1day
    period: annual
    required:
        - pr
    spell:
        value: 1
        unit: day
    logic: max
    cf_name: rx1dayETCCDI
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
    cf_name: rx5dayETCCDI
# %% annual_maximum_5day_precipitation
annual_maximum_5day_precipitation:
    name: Rx5day
    period: annual
    required:
        - pr
    spell:
        value: 5
        unit: day
    logic: max
    cf_name: rx5dayETCCDI
# %% annual_total_precipitation_in_wet_days
annual_total_precipitation_in_wet_days:
    name: PRCPTOT
    period: annual
    required:
        - pr
    logic: sum
    cf_name: prcptotETCCDI
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
    cf_name: r95pETCCDI
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
    cf_name: r99pETCCDI
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
    cf_name: tn10pETCCDI
# %% annual_number_of_days_where_daily_minimum_temperature_below_10%
annual_number_of_days_where_daily_minimum_temperature_below_10%:
    name: TN10p
    period: annual
    required:
        - tasmin
    threshold:
        value: 10
        unit: percent
        logic: less
    cf_name: tn10pETCCDI
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
    cf_name: tn90pETCCDI
# %% annual_number_of_days_where_daily_minimum_temperature_above_90%
annual_number_of_days_where_daily_minimum_temperature_above_90%:
    name: TN90p
    period: annual
    required:
        - tasmin
    threshold:
        value: 90
        unit: percent
        logic: greater
    cf_name: tn90pETCCDI
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
    cf_name: tx10pETCCDI
# %% annual_number_of_days_where_daily_maximum_temperature_below_10%
annual_number_of_days_where_daily_maximum_temperature_below_10%:
    name: TX10p
    period: annual
    required:
        - tasmax
    threshold:
        value: 10
        unit: percent
        logic: less
    cf_name: tx10pETCCDI
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
    cf_name: tx90pETCCDI
# %% annual_number_of_days_where_daily_maximum_temperature_above_90%
annual_number_of_days_where_daily_maximum_temperature_above_90%:
    name: TX90p
    period: annual
    required:
        - tasmax
    threshold:
        value: 90
        unit: percent
        logic: greater
    cf_name: tx90pETCCDI
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
    cf_name: wsdiETCCDI 
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
    cf_name: cddiETCCDI
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
    cf_name: cddETCCDI
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
    cf_name: cwdETCCDI
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
    cf_name: sdiiETCCDI
# %% annual_simple_precipitation_intensity_index
annual_simple_precipitation_intensity_index:
    name: SDII
    period: annual
    required:
        - pr
    threshold:
        value: 1
        unit: mm day-1
        logic: greater_equal
    cf_name: sdiiETCCDI
# %% daily_temperature_range
daily_temperature_range:
    name: DTR
    period: daily
    required:
        - tasmin
        - tasmax
    logic: diff
    cf_name: dtrETCCDI
# %% annual_growing_season_length
annual_growing_season_length:
    name: GSL
    required:
        - tas
    start:
        threshold:
            value: 5
            unit: celsius
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
            value: 5
            unit: celsius
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
    cf_name: gslETCCDI
""")
#print("INDEX_DEFINITION:")
#print(yaml.dump(index_definition))
