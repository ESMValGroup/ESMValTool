;
;  Requires: none
;
variable_info = True
variable_info@derived = False

; Countour levels for zonalmean_profile.ncl
variable_info@zonal_ref_cmap = "temp_19lev"
variable_info@zonal_ref_Levels = ispan(-70, 70, 10)
variable_info@zonal_diff_Levels = ispan(-20, 20, 5)
variable_info@long_name = "Eastward Wind (m s-1)"
