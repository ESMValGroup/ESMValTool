diag_script_info = True

;; Used by: perfmetrics_main.ncl

;; Plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon))
diag_script_info@plot_type = "cycle_latlon"

;; Time average ("opt" argument of time_operations.ncl)
diag_script_info@time_avg = "monthlyclim"

;; Required fraction of valid values
diag_script_info@valid_fraction = 0.1

;; Selected level (hPA, "all" for no selection)
diag_script_info@level = "all"

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

;; Used by: perfmetrics_grading.ncl

;; MultiModelMean
diag_script_info@MultiModelMean = True

;; MultiModelMedian
diag_script_info@MultiModelMedian = True

;; Metric ("RMSD", "BIAS")
diag_script_info@metric = "RMSD"  

;; Normalization ("mean", "median", "centered_median")
diag_script_info@normalization = "centered_median"
