diag_script_info = True

;; Used by: clouds.ncl
diag_script_info@projection = "CylindricalEquidistant"   ; map projection, e.g., Mollweide, Mercator
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@colormap = "WhiteBlueGreenYellowRed"  ; e.g., WhiteBlueGreenYellowRed, rainbow
diag_script_info@ncdf = "default" ;enable to output to netCDF; either use "default" or give a full file name

; --- multi-model mean parameters -----------------------------
; calculate multi-model mean?
; valid options: True, False
diag_script_info@calcmm = False
diag_script_info@showdiff = True

diag_script_info@embracesetup = False      ; True  = 2 plots per line
                                           ; False = 4 plots per line

; --- time averaging ---
diag_script_info@timemean = "annualclim"   ; options: 1) annualclim (1 plots: annual mean)
                                           ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest")
diag_script_info@grid = "ref"

; reference model (to be skipped when calculating the multi-model mean)

diag_script_info@ref_model = "n.a."

if (variables(0).eq."lwp") then
    diag_script_info@ref_model = (/"UWisc"/)
end if
if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut","SW_CRE","LW_CRE","NET_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
end if
if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"MODIS"/)
end if
if (any(variables(0).eq.(/"clivi"/))) then
;    diag_script_info@ref_model = (/"MODIS-CFMIP"/)
    diag_script_info@ref_model = (/"MODIS"/)
end if
if (any(variables(0).eq.(/"pr","pr-mmday"/))) then
    diag_script_info@ref_model = (/"GPCP-SG"/)
end if
if (any(variables(0).eq.(/"prc","prc-mmh"/))) then
    diag_script_info@ref_model = (/"ERA-Interim"/)
end if
if (any(variables(0).eq.(/"tas", "tas-degC"/))) then
    diag_script_info@ref_model = (/"ERA-Interim"/)
end if
if (any(variables(0).eq.(/"prw"/))) then
    diag_script_info@ref_model = (/"SSMI-MERIS"/)
end if
; -------------------------------------------------------------

