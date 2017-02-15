diag_script_info = True

; used by clouds_taylor.ncl
diag_script_info@styleset = "CMIP5"         ; "CMIP5"
diag_script_info@ncdf = "default" ;enable to output to netCDF; either use "default" or give a full file name

; --- multi-model mean parameters -----------------------------
; calculate multi-model mean?
diag_script_info@calcmm = True

; --- time averaging ---
diag_script_info@timemean = "annualclim"   ; options: 1) annualclim (1 plots: annual mean)
                                           ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; target grid for regridding
; (valid options are: "1x1")
diag_script_info@grid = "1x1"

; region
diag_script_info@region = "Arctic"

; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5
; -------------------------------------------------------------

diag_script_info@embracelegend = True

; reference model

if (variables(0).eq."lwp") then
    diag_script_info@ref_model = (/"UWisc"/)
;    diag_script_info@ref_model = (/"MODIS-L3"/)
end if
if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut", "SW_CRE", "LW_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
end if
if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"ESACCI-CLOUD"/)
end if
if (any(variables(0).eq.(/"clivi"/))) then
    diag_script_info@ref_model = (/"MODIS-L3"/)
end if
if (any(variables(0).eq.(/"pr","pr-mmday"/))) then
    diag_script_info@ref_model = (/"GPCP-SG"/)
end if
if (any(variables(0).eq.(/"prw"/))) then
    diag_script_info@ref_model = (/"SSMI-MERIS"/)
end if
if (any(variables(0).eq.(/"od550aer", "od870aer"/))) then
    diag_script_info@ref_model = (/"ESACCI-AEROSOL"/)
end if
if (any(variables(0).eq.(/"dos"/))) then
    diag_script_info@ref_model = (/"ESACCI-SOILMOISTURE"/)
end if
if (any(variables(0).eq.(/"tas"/))) then
    diag_script_info@ref_model = (/"ERA-Interim"/)
end if
if (any(variables(0).eq.(/"toz"/))) then
    diag_script_info@ref_model = (/"ESACCI-OZONE"/)
end if
if (any(variables(0).eq.(/"sic"/))) then
   diag_script_info@ref_model = (/"ESACCI-SSMI"/)
end if


