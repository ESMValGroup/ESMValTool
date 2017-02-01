diag_script_info = True

;; Used by: perfmetrics_main.ncl

;; Plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon))
diag_script_info@plot_type = "cycle_latlon"

;; Time average ("opt" argument of time_operations.ncl)
diag_script_info@time_avg = "monthlyclim"

;; Required fraction of valid values
diag_script_info@valid_fraction = 0.95

;; Selected level (hPA, "all" for no selection)
diag_script_info@level = 500.

;; Selected region (Global, Tropics, NH extratropics, SH extratropics)
diag_script_info@region = "Global"

;; Selected grid ("ref", "coarest", "finest")
diag_script_info@grid = "ref"

;; Draw plots
diag_script_info@draw_plots = False

;; Draw difference plots
diag_script_info@plot_diff = True

;; Plot standard deviation ("all", "none", "ref_model" or given model name)
diag_script_info@plot_stddev = "all"

;; Plot legend in a separate file, for "cycle"
diag_script_info@legend_outside = False

;; Plot style (EMAC, CMIP5, DEFAULT), for "cycle"
diag_script_info@styleset = "CMIP5"

;; Calculate t-test in difference plot, for "zonal" and "lat-lon"
diag_script_info@t_test = False

;; Confidence level for the t-test
diag_script_info@conf_level = 0.95

;; Reference model, as specified in annotations
;; If the attribute is a 2D array, the first model specified the
;; reference, the second model the alternative
if (variables(0).eq."hus") then
    diag_script_info@ref_model = (/"AIRS", "ERA-Interim"/)
else if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut", "SW_CRE", "LW_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
else if (any(variables(0).eq.(/"ta", "ua", "va", "zg", "tas"/))) then
    diag_script_info@ref_model = (/"ERA-Interim", "NCEP"/)
else if (variables(0).eq."pr") then
    diag_script_info@ref_model = (/"GPCP-SG"/)
else if (variables(0).eq."toz") then
    diag_script_info@ref_model = (/"ESACCI-OZONE", "NIWA"/)
end if
end if
end if
end if
end if

;; Used by: perfmetrics_grading.ncl

;; MultiModelMean
diag_script_info@MultiModelMean = True

;; MultiModelMedian
diag_script_info@MultiModelMedian = True

;; Metric ("RMSD", "BIAS")
diag_script_info@metric = "RMSD"  

;; Normalization ("mean", "median", "centered_median")
diag_script_info@normalization = "centered_median"
