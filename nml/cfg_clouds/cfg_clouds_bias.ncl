diag_script_info = True

; used by clouds_bias.ncl
;diag_script_info@projection = "CylindricalEquidistant"   ; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   ; map projection
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@colormap = "WhiteBlueGreenYellowRed"  ; e.g., WhiteBlueGreenYellowRed, rainbow
diag_script_info@ncdf = "default" ;enable to output to netCDF; either use "default" or give a full file name

; target grid for regridding: use grid of reference model
diag_script_info@grid = "ref"

; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5

; --- time averaging ---
diag_script_info@timemean = "annualclim"   ; options: 1) annualclim (1 plots: annual mean)
                                           ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; reference model

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
else if (any(variables(0).eq.(/"tas", "tas-degC"/))) then
    diag_script_info@ref_model = (/"ERA-Interim"/)
end if
end if
end if
end if
end if
end if

