;
;  Requires: none
;
variable_info = True
variable_info@derived = False

; Countour levels for zonalmean_profile.ncl
variable_info@zonal_ref_cmap = "temp_19lev"
variable_info@zonal_ref_Levels = ispan(-15, 15, 2)/10.
variable_info@zonal_diff_Levels = ispan(-10, 10, 1)/10.
variable_info@long_name = "Northward Wind (m s-1)"
