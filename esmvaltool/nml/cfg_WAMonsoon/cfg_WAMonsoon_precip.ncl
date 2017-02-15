; Settings used by most WAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJAS"

; Settings specific for the precip basic plot_type
diag_script_info@latrange_basic = (/-10, 25/)
diag_script_info@lonrange_basic = (/-30, 50/)

diag_script_info@diff_colormap = "diag_scripts/aux/WAMonsoon/cmap_difference.rgb"

diag_script_info@cn_levels_mean_basic = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19/)
diag_script_info@cn_levels_mean_diff_basic = ispan(-11, 13, 2)

diag_script_info@cn_levels_stddev_basic = (/0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00/)
diag_script_info@cn_levels_stddev_diff_basic = (/-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0/)

; Settings specific for precip annual cycle
diag_script_info@latrange_seasonal = (/10, 20/)
diag_script_info@lonrange_seasonal = (/-10,10/)

diag_script_info@my_region = "WA"
diag_script_info@styleset = "CMIP5"

diag_script_info@multi_model_mean = "y"
diag_script_info@supporting_gridlines = "y"
diag_script_info@plot_stddev = False
