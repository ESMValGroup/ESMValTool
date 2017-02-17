;
; Requires: none
;
variable_info = True
variable_info@derived = False
variable_info@long_name = "condensed water path"
variable_info@units = "kg m-2"

; Contour levels (in g m-2)
; note that variable will be plotted in kg m-2, so contour level will be divided by 1000
variable_info@map_ref_Levels = ispan(10, 200, 10)
variable_info@map_diff_Levels = ispan(-50, 50, 10)

; Add global average to the plot
variable_info@show_global_avg = True

; misc
variable_info@title = "LWP+IWP"
