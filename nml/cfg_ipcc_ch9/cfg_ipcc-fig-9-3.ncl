diag_script_info = True

diag_script_info@projection = "Robinson"   ; map projection

; target grid for regridding: use grid of reference model
diag_script_info@grid = "ref"

; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5

; reference model

if (any(variables(0).eq.(/"tas", "tas-degC"/))) then
    diag_script_info@ref_model = (/"ERA-Interim"/)
end if

