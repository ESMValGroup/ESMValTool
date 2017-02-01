diag_script_info = True

;; Used by: MyDiag.ncl
diag_script_info@projection = "Mollweide"   ; map projection, e.g., Mollweide, Mercator
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@colormap = "WhiteBlueGreenYellowRed"  ; e.g., WhiteBlueGreenYellowRed, rainbow

;; Or read rgb-file for colormap
; delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                  ;; the colormap here
; colormap_file = "diag_scripts/lib/ncl/rgb/red-blue.rgb"
; diag_script_info@colormap = read_colormap_file(colormap_file)

;; Set the contour levels explicitly
; diag_script_info@explicit_cn_levels = (/205, 206, 207, 208, 209, 210,\
;                                       211, 212, 213, 214, 215, 216,\
;                                       217, 218, 219/)
