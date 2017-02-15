;
;  Requires: none
;
variable_info = True
variable_info@derived = False

; Countour levels for zonalmean_profile.ncl
variable_info@zonal_ref_Levels = ispan(0, 15, 1)
variable_info@zonal_diff_Levels = ispan(-5, 5, 1)

; Minimum pressure for zonalmean_profile.ncl
variable_info@zonal_ymin = 300.

;; Plot units
variable_info@plot_units = "g/kg"
