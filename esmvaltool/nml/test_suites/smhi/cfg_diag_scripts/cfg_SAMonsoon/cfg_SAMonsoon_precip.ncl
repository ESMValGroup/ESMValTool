; Settings used by most SAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJAS"

; Settings specific for the precip basic plot_type
diag_script_info@latrange_basic = (/-20, 40/)
diag_script_info@lonrange_basic = (/30, 150/)

diag_script_info@cn_levels_mean_basic = (/2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5/)
diag_script_info@cn_levels_mean_diff_basic = (/-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15/)
diag_script_info@cn_levels_stddev_basic = (/0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0/)
diag_script_info@cn_levels_stddev_diff_basic = (/-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7/)
diag_script_info@cn_levels_stddev_norm_basic  = (/0.00, 0.13, 0.26, 0.39, 0.52, 0.65, 0.78, 0.91, 1.04, 1.17, 1.30, 1.43, 1.56, 1.69, 1.82, 1.95/)
diag_script_info@cn_levels_stddev_norm_diff_basic  = (/-2.00, -1.87, -1.74, -1.61, -1.48, -1.35, -1.22, -1.09, -0.96, -0.83, -0.70, -0.57, -0.44, -0.31, -0.18, -0.05, 0.08, 0.21, 0.34, 0.47, 0.60, 0.73, 0.86, 0.99, 1.12, 1.25, 1.38, 1.51, 1.64, 1.77, 1.90/)

; Settings specific for the precip seasonal plot_type
diag_script_info@latrange_seasonal = (/5, 30/)
diag_script_info@lonrange_seasonal = (/65, 95/)
diag_script_info@colors_seasonal = (/"red",              "green",          "blue",\
                                     "yellow",           "purple",         "white",\
                                     "pink",             "violet",         "gray",\
                                     "orange",           "black",          "burlywood",\
                                     "DarkOrange4",      "DarkSeaGreen3",  "DodgerBlue2",\
                                     "gainsboro",        "goldenrod3",     "honeydew2",\
                                     "IndianRed2",       "khaki3",         "LemonChiffon3",\
                                     "LightCyan3",       "LightPink1",     "LightSkyBlue1",\
                                     "LightSteelBlue4",  "magenta4",       "MediumOrchid2",\
                                     "OliveDrab4",       "OrangeRed4",     "PaleGreen3",\
                                     "PaleVioletRed3",   "pink2",          "purple1",\
                                     "RosyBrown1",       "salmon"/)
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@multi_model_mean = "y"
diag_script_info@supporting_gridlines = "y"

; Settings shared by the precip global domain/intensity plot_type
diag_script_info@summer_season = "MJJAS"
diag_script_info@winter_season = "NDJFM"

; Settings specific for the precip global domain plot_type
diag_script_info@latrange_global = (/-50, 50/)
diag_script_info@cn_levels_global_domain_mean = (/0.01, 1.01/)
diag_script_info@cn_levels_global_domain_diff = (/-0.5, 0.5, 0.51/)

diag_script_info@colors_global_domain_mean = (/"white", "blue"/)
diag_script_info@colors_global_domain_diff = (/"red", "white", "blue"/)

diag_script_info@high_intensity_cutoff_mmday = 2.0


; Settings specific for the precip intensity plots
diag_script_info@cn_levels_intensity_mean = (/-16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16/)
diag_script_info@cn_levels_intensity_diff = (/-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8/)
diag_script_info@abs_cmap = "diag_scripts/aux/SAMonsoon/precip_domain_absolute.rgb"
diag_script_info@diff_cmap = "diag_scripts/aux/SAMonsoon/precip_domain_difference.rgb"

diag_script_info@my_region = "SA"
