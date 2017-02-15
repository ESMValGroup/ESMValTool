; Settings used by daily WAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJA"

diag_script_info@destgrid = "2.5x2.5" 

; Style
diag_script_info@styleset = "CMIP5"

; Settings specific for the contour plot limits
diag_script_info@latrange_basic = (/-50, 50/)
diag_script_info@lonrange_basic = (/-180, 180/)

; Attributes with dashes requires the special syntax: $"...string..."$
diag_script_info@cn_levels_stddev_basic = ispan(0,30,3) 
diag_script_info@cn_levels_stddev_diff_basic = ispan(-15,15,3)

diag_script_info@plot_norm = False
diag_script_info@cn_levels_stddev_norm_basic  = (/0.00, 0.13, 0.26, 0.39, 0.52, 0.65, 0.78, 0.91, 1.04, 1.17, 1.30, 1.43, 1.56, 1.69, 1.82, 1.95/)
diag_script_info@cn_levels_stddev_norm_diff_basic = (/-2.00, -1.87, -1.74, -1.61, -1.48, -1.35, -1.22, -1.09, -0.96, -0.83, -0.70, -0.57, -0.44, -0.31, -0.18, -0.05, 0.08, 0.21, 0.34, 0.47, 0.60, 0.73, 0.86, 0.99, 1.12, 1.25, 1.38, 1.51, 1.64, 1.77, 1.90/)


