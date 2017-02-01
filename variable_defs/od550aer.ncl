;
; Requires: none
;
variable_info          = True
variable_info@derived  = False

; Contour levels for absolute value maps
variable_info@contours_abs = fspan(0,0.5,11)

; Color table for absolute value maps
variable_info@colortable_abs = "WhBlGrYeRe"

; Contour levels for difference maps
variable_info@contours_diff = ispan(-30,30,5)/100.

; Color table for difference maps
variable_info@colortable_diff = "hotcolr_19lev"
