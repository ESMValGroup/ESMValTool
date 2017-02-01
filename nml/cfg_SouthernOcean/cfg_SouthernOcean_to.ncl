diag_script_info              = True

;; Used by: SouthernOcean_areamean_vertconplot.ncl
diag_script_info@region       = "Southern Ocean"                        ; Region name
diag_script_info@West         =   0.0                                   ; Western limit for area
diag_script_info@East         = 360.0                                   ; Eastern limit for area
diag_script_info@South        = -70.0                                   ; Southern limit for area
diag_script_info@North        = -60.0                                   ; Northern limit for area
diag_script_info@styleset     = "CMIP5"                                 ; Plot style: "CMIP5", "DEFAULT"
diag_script_info@grid_min     = -1.7                                    ; min for plot scale
diag_script_info@grid_max     =  1.7                                    ; max for plot scale
diag_script_info@grid_step    =  0.1                                    ; step size for plot scale
diag_script_info@colormap     = "temp_diff_18lev"                       ; e.g., WhiteBlueGreenYellowRed, rainbow
diag_script_info@dbg_plots    = False                                   ; create additional plots for debugging purposes

;; Or read rgb-file for colormap
; delete(diag_script_info@colormap)                                       ; if defined above it is necessary to delete the colormap here
; colormap_file               = "diag_scripts/lib/ncl/red-blue.rgb"
; diag_script_info@colormap   = read_colormap_file(colormap_file)

