diag_script_info = True

;; Used by: all SeaIce routines
diag_script_info@region       = "Antarctic" ; "Arctic" or "Antarctic"; entire hemisphere will be evaluated
                                            ; Make sure to specify data for the correct hemisphere in the main nml.
diag_script_info@month        = "A"         ; A = annual mean, 3 = March, 9 = September
diag_script_info@styleset     = "CMIP5"     ; "CMIP5", "DEFAULT"

;; Used by: SeaIce_polcon, SeaIce_polcon_diff
diag_script_info@contour_extent = "n"       ; y = draw a red contour line for sic extent in polar stereographic plots
diag_script_info@max_vert     =    2        ; max allowed number of rows on a panel page (vertical)
diag_script_info@max_hori     =    1        ; max allowed number of columns on a panel page (horizontal)
diag_script_info@max_lat      =  -40.       ; Antarctic plotted from -90 degrees up to this latitude
diag_script_info@min_lat      =   60.       ; Arctic plotted from 90 degrees up to this latitude
diag_script_info@PanelTop     =    0.99     ; tune down to get full title of uppermost panel row (default 0.99)
diag_script_info@grid_min     =   -0.16     ; min contour value, default: -1.
diag_script_info@grid_max     =    0.16     ; max contour value, default: 1.
diag_script_info@grid_step    =    0.02     ; step between contours, default: 0.2
diag_script_info@grid_center  =    0.0      ; value to center the color bar, default: 0.0
diag_script_info@colormap     = "nrl_sirkes"; name of colour table to use
diag_script_info@PanelLabelBar= True        ; use single label bar per page
diag_script_info@showunits    = True        ; display units in figure title
diag_script_info@range_option = 1           ; 0 = use each model's time range specified in namelist
                                            ; 1 = use only intersection of all time ranges 

;; Used by: SeaIce_polcon_diff.ncl
