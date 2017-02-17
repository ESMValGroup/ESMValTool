plot_type_info = True

;; Used by: all SeaIce routines
plot_type_info@region = "Antarctic"       ; "Arctic" or "Antarctic"; entire hemisphere will be evaluated
                               ; Make sure to specify data for the correct hemisphere in the main nml.
plot_type_info@month = "3"                ; A = annual mean, 3 = March, 9 = September
plot_type_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"

;; Used by: SeaIce_tsline.ncl, SeaIce_ancyc.ncl
;; see style_FUNCTIONS/style_$project$.ncl for which data shall be included in the statistics
plot_type_info@multi_model_mean = "y"     ; "y" = plot multi-model mean & stddev; "n" = don't

;; Used by: SeaIce_polcon, SeaIce_polcon_diff
plot_type_info@contour_extent = "y"       ; y = draw a red contour line for sic extent in polar stereographic plots
plot_type_info@max_vert = 1               ; max allowed number of rows on a panel page (vertical)
plot_type_info@max_hori = 1               ; max allowed number of columns on a panel page (horizontal)
plot_type_info@max_lat = -50.             ; Antarctic plotted from -90° up to this latitude
plot_type_info@min_lat = 60.              ; Arctic plotted from 90° up to this latitude
;plot_type_info@PanelTop = 0.99            ; tune down to get full title of uppermost panel row (default 0.99)

;; Used by: SeaIce_polcon_diff.ncl
plot_type_info@ref_model = "MPI-ESM-LR"   ; reference model, as specified in annotations -> if the specified string is
                                     ;  not found, the routine will print a list of valid strings before crashing
;plot_type_info@grid_min    = -0.3   ;; min contour value, default: -1.
;plot_type_info@grid_max    =  0.3   ;; max contour value, default: 1.
;plot_type_info@grid_step   =  0.05  ;; step between contours, default: 0.2
;plot_type_info@grid_center =  0.0  ;; value to center the color bar, default: 0.0
