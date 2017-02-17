plot_type_info = True

;; Used by: MyDiag.ncl
plot_type_info@projection = "Mollweide"   ; map projection, e.g., Mollweide, Mercator
plot_type_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
plot_type_info@colormap = "WhiteBlueGreenYellowRed"  ; e.g., WhiteBlueGreenYellowRed, rainbow

;; Or read rgb-file for colormap
; delete(plot_type_info@colormap) ;; if defined above it is necessary to delete
                                  ;; the colormap here
; colormap_file = "./diag_scripts/lib/ncl/rgb/red-blue.rgb"
; plot_type_info@colormap = read_colormap_file(colormap_file)

;; Set the contour levels explicitly
; plot_type_info@explicit_cn_levels = (/205, 206, 207, 208, 209, 210,\
;                                       211, 212, 213, 214, 215, 216,\
;                                       217, 218, 219/)
