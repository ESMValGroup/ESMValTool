;
;  Requires: none
;
variable_info = True
variable_info@derived = False
variable_info@long_name = "TOA Outgoing Shortwave Radiation Clear Sky"

;ymin,ymax info for DiurnalCycle_box.ncl diag_script
variable_info@minmax = (/0,500./)

; Contour levels
variable_info@map_ref_Levels = ispan(100, 320, 20)
variable_info@map_diff_Levels = ispan(-30, 30, 5)

; Add global average to the plot
variable_info@show_global_avg = True
