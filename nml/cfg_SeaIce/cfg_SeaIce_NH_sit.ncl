diag_script_info = True

;; Used by: all SeaIce routines
diag_script_info@region = "Arctic"          ; "Arctic" or "Antarctic"; entire hemisphere will be evaluated
	                                    ; Make sure to specify data for the correct hemisphere in the main nml.
diag_script_info@month = "9"                ; A = annual mean, 3 = March, 9 = September
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"

;; Used by: SeaIce_polcon, SeaIce_polcon_diff, SeaIce_ancyc
diag_script_info@range_option = 1           ; 0 = use each model's whole time range as specified in namelist
                                            ; 1 = use only intersection of all time ranges 

;; Used by: SeaIce_tsline, SeaIce_ancyc
;; see style_FUNCTIONS/style_$project$.ncl for which data shall be included in the statistics
diag_script_info@multi_model_mean = "n"     ; "y" = plot multi-model mean & stddev; "n" = don't
diag_script_info@EMs_in_lg = False          ; create legend label for each individual ensemble member
diag_script_info@fill_pole_hole = True      ; Fill pole hole in data
diag_script_info@legend_outside = False

;; Used by: SeaIce_polcon, SeaIce_polcon_diff
diag_script_info@contour_extent = "n"       ; y = draw a red contour line for sic extent in polar stereographic plots
diag_script_info@max_vert = 1               ; max allowed number of rows on a panel page (vertical)
diag_script_info@max_hori = 1               ; max allowed number of columns on a panel page (horizontal)
diag_script_info@max_lat = -50.             ; Antarctic plotted from -90° up to this latitude
diag_script_info@min_lat = 60.              ; Arctic plotted from 90° up to this latitude
diag_script_info@showunits = True           ; Show units in title
;diag_script_info@PanelTop = 0.99            ; tune down to get full title of uppermost panel row (default 0.99)

;; Used by: SeaIce_polcon_diff.ncl
diag_script_info@ref_model = "PIOMAS"   ; reference model, as specified in annotations -> if the specified string is
                                        ; not found, the routine will print a list of valid strings before crashing
diag_script_info@dst_grid = "./diag_scripts/aux/CDO/cdo_dst_grid_g010" 
                                        ; path to destination grid file, required by cdo_remapdis
;diag_script_info@grid_min    = -0.3   	; min contour value, default: -1.
;diag_script_info@grid_max    =  0.3   	; max contour value, default: 1.
;diag_script_info@grid_step   =  0.05  	; step between contours, default: 0.2
;diag_script_info@grid_center =  0.0  	; value to center the color bar, default: 0.0
