diag_script_info = True

; used by clouds_taylor.ncl
diag_script_info@styleset = "CMIP5"         ; "CMIP5"

; --- multi-model mean parameters -----------------------------
; calculate multi-model mean?
diag_script_info@calcmm = True

; --- time averaging ---
diag_script_info@timemean = "annualclim"   ; options: 1) annualclim (1 plot: annual mean)
                                           ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest", "1x1")
diag_script_info@grid = "1x1"

; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5

; mask T < 272 K as sea ice for variable "ts"
diag_script_info@mask_ts_sea_ice = True

; -------------------------------------------------------------

diag_script_info@embracelegend = True

; data sets to be skipped when calculating the multi-model mean
; (= observations)
diag_script_info@models_to_skip = (/"UWisc", "CERES-EBAF", "MODIS", "MODIS-CFMIP", "GPCP-SG", \
             "ERA-Interim", "NCEP", "HadISST", "NIWA", "PATMOS", "CLARA-A2", "ESACCI-CLOUD",  \
             "ESACCI-OZONE", "ESACCI-GHG", "ESACCI-SST", "ESACCI-SOILMOISTURE", "ESACCI-AEROSOL"/)

