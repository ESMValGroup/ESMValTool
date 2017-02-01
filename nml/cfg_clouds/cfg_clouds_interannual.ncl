diag_script_info = True

;; Used by: clouds.ncl
diag_script_info@projection = "CylindricalEquidistant"   ; map projection, e.g., Mollweide, Mercator
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@colormap = "WhiteBlueGreenYellowRed"  ; e.g., WhiteBlueGreenYellowRed, rainbow
diag_script_info@ncdf = "default" ;enable to output to netCDF; either use "default" or give a full file name

; --- multi-model mean parameters -----------------------------
; calculate multi-model mean?
; valid options: True, False
diag_script_info@calcmm = True

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest")
diag_script_info@grid = "ref"

; reference model (to be skipped when calculating the multi-model mean)

if (variables(0).eq."lwp") then
    diag_script_info@ref_model = (/"UWisc"/)
else if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut","SW_CRE","LW_CRE","NET_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
else if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"MODIS"/)
else if (any(variables(0).eq.(/"clivi"/))) then
    diag_script_info@ref_model = (/"MODIS-CFMIP"/)
else if (any(variables(0).eq.(/"pr","pr-mmday"/))) then
    diag_script_info@ref_model = (/"GPCP-SG"/)
end if
end if
end if
end if
end if
; -------------------------------------------------------------

