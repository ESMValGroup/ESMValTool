diag_script_info = True

; used by clouds_ipcc.ncl
diag_script_info@projection = "Robinson"   ; map projection
diag_script_info@styleset = "CMIP5"         ; "CMIP5"

; target grid for regridding: use grid of reference model
diag_script_info@grid = "1x1"

; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5

; --- time averaging ---
diag_script_info@timemean = "seasonalclim"   ; options: 1) annualclim (1 plots: annual mean)
                                             ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; reference model

if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"ESACCI-CLOUD"/)
end if

; data sets to be skipped when calculating the multi-model mean
; (= observations)

diag_script_info@models_to_skip = (/"UWisc", "CERES-EBAF", "MODIS", "MODIS-CFMIP", "GPCP-SG", \
             "ERA-Interim", "PATMOS", "CLARA-A2", "ESACCI-CLOUD"/)
