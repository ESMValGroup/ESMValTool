diag_script_info = True

;; Used by: perfmetrics_main.ncl

;; Plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon))
diag_script_info@plot_type = "cycle_latlon"

;; Time average ("opt" argument of time_operations.ncl)
diag_script_info@time_avg = "monthlyclim"

;; Required fraction of valid values
diag_script_info@valid_fraction = 0.95

;; Selected level (hPA, "all" for no selection)
diag_script_info@level = 200.

;; Selected region ("Global", "Tropics", "NH extratropics", "SH extratropics")
diag_script_info@region = "NH extratropics"

;; Regridding option
diag_script_info@grid = "coarsest"

;; Draw plots
diag_script_info@draw_plots = False

;; Draw difference plots
diag_script_info@plot_diff = False

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

;; Used by perfmetrics_grading.ncl

;; Metric ("RMSD", "BIAS")
diag_script_info@metric = "BIAS"  

;; Normalization ("mean", "median", "centered_median")
diag_script_info@normalization = "mean"
