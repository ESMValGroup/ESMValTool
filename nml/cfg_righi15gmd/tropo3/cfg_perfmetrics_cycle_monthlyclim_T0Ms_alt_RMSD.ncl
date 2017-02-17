diag_script_info = True

;; Used by: perfmetrics_main.ncl

;; Plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon))
diag_script_info@plot_type = "cycle"

;; Time average ("opt" argument of time_operations.ncl)
diag_script_info@time_avg = "monthlyclim"

;; Required fraction of valid values
diag_script_info@valid_fraction = 0.50

;; Selected level (hPA, "all" for no selection)
diag_script_info@level = "all"

;; Annotation
diag_script_info@location = "Alert (82N 62W)"
diag_script_info@altitude = "210 m"                
;diag_script_info@aux_info = ""

;; Draw plots
diag_script_info@draw_plots = True

;; Draw difference plots
diag_script_info@plot_diff = True

;; Plot standard deviation ("all", "none", "ref_model" or given model name)
diag_script_info@plot_stddev = "all"

;; Plot legend in a separate file (for cycle)
diag_script_info@legend_outside = True

;; Plot style (for cycle)
diag_script_info@styleset = "righi15gmd"

;; Calculate t-test in difference plot (for zonal and lat-lon)
diag_script_info@t_test = False

;; Confidence level for the t-test (for zonal and lat-lon)
diag_script_info@conf_level = 0.95

;; Used by: perfmetrics_grading.ncl

;; Metric ("RMSD", "BIAS")
diag_script_info@metric = "RMSD"

;; Normalization ("mean", "median", "centered_median")
diag_script_info@normalization = "mean"
