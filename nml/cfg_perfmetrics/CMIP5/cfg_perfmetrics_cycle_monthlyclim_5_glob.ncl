diag_script_info = True

;; Used by: perfmetrics_main.ncl

;; Plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon))
diag_script_info@plot_type = "cycle"

;; Time average ("opt" argument of time_operations.ncl)
diag_script_info@time_avg = "monthlyclim"

;; Required fraction of valid values
diag_script_info@valid_fraction = 0.95

;; Selected level (hPA, "all" for no selection)
diag_script_info@level = 5.

;; Selected region ("Global", "Tropics", "NH extratropics", "SH extratropics")
diag_script_info@region = "Global"

;; Regridding option
diag_script_info@grid = "ref"

;; Draw plots
diag_script_info@draw_plots = True

;; Draw difference plots
diag_script_info@plot_diff = True

;; Plot standard deviation ("all", "none", "ref_model" or given model name)
diag_script_info@plot_stddev = "ref_model"

;; Plot legend in a separate file (for cycle)
diag_script_info@legend_outside = True

;; Plot style (for cycle)
diag_script_info@styleset = "CMIP5"

;; Calculate t-test in difference plot (for zonal and lat-lon)
diag_script_info@t_test = False

;; Confidence level for the t-test (for zonal and lat-lon)
diag_script_info@conf_level = 0.95

;; Reference model, as specified in annotations
;; If the attribute is a 2D array, the first model specified the
;; reference, the second model the alternative
if (variables(0).eq."hus") then
    diag_script_info@ref_model = (/"AIRS", "ERA-Interim"/)
else if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut"/))) then
    diag_script_info@ref_model = (/"SRB","CERES-EBAF"/)
else if (any(variables(0).eq.(/"ta", "ua", "va", "zg"/))) then
    diag_script_info@ref_model = (/"ERA-Interim", "NCEP"/)
end if
end if
end if
