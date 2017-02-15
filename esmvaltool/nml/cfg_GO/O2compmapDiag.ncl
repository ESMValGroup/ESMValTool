diag_script_info = True
;;diag_script_info@styleset     = "GO"    
diag_script_info@GO_season = "JFM"
diag_script_info@GO_detrend = False
diag_script_info@latrange_basic = (/-90, 90/)
diag_script_info@lonrange_basic = (/0., 360./)
diag_script_info@depthrange_basic =(/98., 106./)

diag_script_info@regrid_for_these_projects = (/"CMIP5_gridfile", "GO_gridfile", "CMIP5_fx"/)
diag_script_info@models_using_i_j_indices = (/"GFDL-ESM2M"/)

diag_script_info@cn_levels_mean_basic = fspan(0, .5, 10)
diag_script_info@cn_levels_mean_diff_basic = fspan(-.16, .16, 32)


diag_script_info@cn_levels_stddev_basic = (/0.0, .002, .004, .006, 0.008, .01, .012, .014, .016, .018, .02 /)
diag_script_info@cn_levels_stddev_diff_basic = (/-.01, -.008,-.006,-.004,-.002, 0, .002, .004, .006, .008, .01/)


;;diag_script_info@plottype_lbres = False
diag_script_info@regrid_for_these_projects = (/"CMIP5_gridfile", "GO_gridfile", "CMIP5_fx"/)
