diag_script_info = True

;; Selected grid ("0.5deg", "1deg", "2deg")
diag_script_info@grid = "2deg"

;; Define region for the spatial average
diag_script_info@region = "Southern Hemisphere"

;; Reference model
if (any(variables(0).eq.(/"tas-degC", "pr-mmday"/))) then
    diag_script_info@ref_model = "CRU"
end if
if (variables(0).eq."GPP") then
    diag_script_info@ref_model = "MTE"
end if
if (variables(0).eq."LAI") then
    diag_script_info@ref_model = "LAI3g"
end if
if (variables(0).eq."NBP") then
    diag_script_info@ref_model = "JMA-TRANSCOM"
end if
