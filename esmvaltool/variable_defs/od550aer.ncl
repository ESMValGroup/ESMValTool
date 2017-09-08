;
; Requires: none
;
variable_info          = True
variable_info@derived  = False

; Contour levels for absolute value maps
variable_info@contours_abs = fspan(0,0.5,21)

; Color table for absolute value maps
;; variable_info@colortable_abs = "WhiteBlueGreenYellowRed"

; Contour levels for difference maps
variable_info@contours_diff = fspan(-0.25,0.25,21)

; Color table for difference maps
;; variable_info@colortable_diff = "hotcolr_19lev"

; For perfmetrics
variable_info@map_ref_cmap = "WhiteBlueGreenYellowRed"
variable_info@map_ref_Levels = ispan(0,100,5)/100.
variable_info@map_diff_Levels = ispan(-50,50,5)/100.
