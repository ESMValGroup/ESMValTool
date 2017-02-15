; Settings used by daily SAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJAS"

; Style
diag_script_info@styleset = "CMIP5"

; Settings specific for the contour plot limits
diag_script_info@latrange_basic = (/-20, 40/)
diag_script_info@lonrange_basic = (/30, 150/)

; Attributes with dashes requires the special syntax: $"...string..."$
diag_script_info@cn_levels_mean_basic = (/2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5/)
diag_script_info@cn_levels_mean_diff_basic = (/-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15/)
diag_script_info@cn_levels_stddev_basic = (/0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0/)
diag_script_info@cn_levels_stddev_diff_basic = (/-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7/)
diag_script_info@cn_levels_stddev_norm_basic  = (/0.00, 0.13, 0.26, 0.39, 0.52, 0.65, 0.78, 0.91, 1.04, 1.17, 1.30, 1.43, 1.56, 1.69, 1.82, 1.95/)
diag_script_info@cn_levels_stddev_norm_diff_basic = (/-2.00, -1.87, -1.74, -1.61, -1.48, -1.35, -1.22, -1.09, -0.96, -0.83, -0.70, -0.57, -0.44, -0.31, -0.18, -0.05, 0.08, 0.21, 0.34, 0.47, 0.60, 0.73, 0.86, 0.99, 1.12, 1.25, 1.38, 1.51, 1.64, 1.77, 1.90/)
diag_script_info@cn_levels_propagation = (/-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7/)

; Settings specific for the daily precip plots
diag_script_info@isv_latrange = (/-20, 40/)
diag_script_info@isv_lonrange = (/0, 180/)
diag_script_info@bob_latrange = (/10, 20/)
diag_script_info@bob_lonrange = (/80, 100/)
diag_script_info@eio_latrange = (/-10, 10/)
diag_script_info@eio_lonrange = (/80, 100/)
diag_script_info@filter_min = 30
diag_script_info@filter_max = 50
diag_script_info@filter_weights = 401

; Settings specific for the daily propagation plots
diag_script_info@prop_isv_latrange = (/-10, 30/)
diag_script_info@prop_isv_lonrange = (/60, 160/)
diag_script_info@prop_lag_latrange = (/-10, 10/)
diag_script_info@prop_lag_lonrange = (/80, 100/)
diag_script_info@prop_ref_latrange = (/-10, 5/)
diag_script_info@prop_ref_lonrange = (/75, 100/)
diag_script_info@cn_levels_rp_timelat_basic = (/-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,\
                                                -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,\
                                                 0.7, 0.8, 0.9, 1.0/)
diag_script_info@cn_levels_rp_timelon_basic = diag_script_info@cn_levels_rp_timelat_basic

diag_script_info@prop_filter_min = 30
diag_script_info@prop_filter_max = 80
diag_script_info@prop_filter_weights = 41
diag_script_info@prop_filter_mxlag = 25
diag_script_info@ihp = 2                     ; bpf=>band pass filter
diag_script_info@nWgt = 161
diag_script_info@sigma = 1.0                 ; Lanczos sigma

diag_script_info@max_days = 80.
diag_script_info@min_days = 30.
diag_script_info@fca = 1./diag_script_info@max_days
diag_script_info@fcb = 1./diag_script_info@min_days
