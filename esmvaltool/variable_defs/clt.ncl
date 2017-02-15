;
; Requires: none
;
variable_info = True
variable_info@derived = False

; Contour levels
variable_info@map_ref_Levels = ispan(0, 100, 10)
variable_info@map_diff_Levels = ispan(-30, 30, 5)

; Add global average to the plot
variable_info@show_global_avg = True

variable_info@long_name = "total cloud cover"
variable_info@units = "%"
variable_info@title = "total CA"

