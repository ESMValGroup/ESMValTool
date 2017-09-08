;
;  Requires: none
;
variable_info = True
variable_info@derived = False
variable_info@long_name = "LWdnSFCclr"
variable_info@units = "W/m2"

;ymin,ymax info for DiurnalCycle_box.ncl diag_script
variable_info@minmax = (/280.,460./)

; Contour levels
variable_info@map_ref_Levels = ispan(100, 320, 20)
variable_info@map_diff_Levels = ispan(-30, 30, 5)

; Add global average to the plot
variable_info@show_global_avg = True
