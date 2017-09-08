;
;  Requires: none
;
variable_info = True
variable_info@derived = False
variable_info@long_name = "surface temperature"
variable_info@skip_long_name_in_subheader = True

; Contour levels
variable_info@map_ref_Levels = ispan(240, 300, 3)
variable_info@map_diff_Levels = ispan(-50, 50, 10) * 0.1
