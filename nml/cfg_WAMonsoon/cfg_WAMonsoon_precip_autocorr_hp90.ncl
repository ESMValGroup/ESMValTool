; Settings used by daily WAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJAS"

diag_script_info@destgrid = "1x1"

; Style
diag_script_info@styleset = "CMIP5"

; Settings specific for the contour plot limits
diag_script_info@latrange_basic = (/-10, 25/)
diag_script_info@lonrange_basic = (/-30, 50/)

; Attributes with dashes requires the special syntax: $"...string..."$
diag_script_info@cn_levels_corr_basic = (/-1., -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6, 0.8, 1., 1.2/)
diag_script_info@cn_levels_corr_diff_basic = (/-1., -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6, 0.8, 1/)

; Settings specific for the daily precip plots
diag_script_info@filter_type = "hp" ; hp or bp
diag_script_info@filter_hp = 90
;diag_script_info@filter_min = 10
;diag_script_info@filter_max = 90
diag_script_info@filter_weights = 301

diag_script_info@sahel_latrange = (/5, 15/)
diag_script_info@sahel_lonrange = (/-10, 10/)

diag_script_info@my_region = "WA"
