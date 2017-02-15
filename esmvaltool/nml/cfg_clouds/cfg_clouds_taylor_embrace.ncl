diag_script_info = True

; used by clouds_taylor.ncl
;diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@ncdf = "default" ;enable to output to netCDF; either use "default" or give a full file name

; --- multi-model mean parameters -----------------------------
; calculate multi-model mean?
diag_script_info@calcmm = True

; --- time averaging ---
diag_script_info@timemean = "annualclim"   ; options: 1) annualclim (1 plots: annual mean)
                                           ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest")
diag_script_info@grid = "ref"

; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5
; -------------------------------------------------------------

diag_script_info@embracelegend = True

; reference model

if (variables(0).eq."lwp") then
    diag_script_info@ref_model = (/"UWisc"/)
;    diag_script_info@ref_model = (/"MODIS-L3"/)
else if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut", "SW_CRE", "LW_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
else if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"MODIS"/)
else if (any(variables(0).eq.(/"clivi"/))) then
;    diag_script_info@ref_model = (/"MODIS-CFMIP"/)
    diag_script_info@ref_model = (/"MODIS"/)
else if (any(variables(0).eq.(/"pr","pr-mmday"/))) then
    diag_script_info@ref_model = (/"GPCP-SG"/)
else if (any(variables(0).eq.(/"prw"/))) then
    diag_script_info@ref_model = (/"SSMI-MERIS"/)
end if
end if
end if
end if
end if
end if
