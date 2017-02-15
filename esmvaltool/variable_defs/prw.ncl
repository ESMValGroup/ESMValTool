;
; Requires: none
;
variable_info = True
variable_info@derived = False
variable_info@long_name = "water vapor path"
variable_info@short_name = "WVP"
variable_info@units = "kg m-2"

; Contour levels
; note that variable will be plotted in kg m-2
variable_info@map_ref_Levels = ispan(2, 50, 25)
variable_info@map_diff_Levels = fspan(-5, 5, 21)
variable_info@colormap = "BlAqGrYeOrReVi200"
variable_info@acp_cnlvl = fspan(2., 50., 13)

; Add global average to the plot
variable_info@show_global_avg = True

; misc
variable_info@title = "WVP"
