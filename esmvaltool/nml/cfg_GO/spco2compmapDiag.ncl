diag_script_info = True
;;diag_script_info@styleset     = "GO"    
diag_script_info@GO_season = "JFMAMJJASOND"
diag_script_info@GO_detrend = True
diag_script_info@latrange_basic = (/-90, 90/)
diag_script_info@lonrange_basic = (/0.,360./)
;;diag_script_info@depth_basic = (/100./)
;;diag_script_info@loncen_basic =(/180./)

diag_script_info@regrid_for_these_projects = (/"CMIP5_gridfile", "GO_gridfile", "CMIP5_fx"/)
diag_script_info@models_using_i_j_indices = (/""/)

;;diag_script_info@cn_levels_mean_basic = (/0, .2,.4,.6,.8,1.,1.2,1.4,1.6,1.8,2./) 
;;diag_script_info@cn_levels_mean_diff_basic = (/ -20,  -15,   -10, -5, 0., 5, 10, 15, 20./)
;;diag_script_info@cn_levels_stddev_basic = (/0.0, 0.05, .1, .15, .2, .25, .3, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 10., 15., 20., 25/)
;;diag_script_info@cn_levels_stddev_diff_basic = (/-23, -20, -15, -10, -6, -5, -4, -3, -2, -1, 0, 1, 2/)

;;diag_script_info@cn_levels_mean_basic = (/0,50,100,150,200,250,300,350,400,450/) 
diag_script_info@cn_levels_mean_basic = fspan(20,70,20) 
diag_script_info@cn_levels_mean_diff_basic = (/ -1,  -.8,  -.6, -.4, -.2, 0., .2, .4, .6, .8, 1./)


diag_script_info@cn_levels_stddev_basic = fspan(0,2,20)
;diag_script_info@cn_levels_stddev_basic = fspan(0,2,20)
diag_script_info@cn_levels_stddev_diff_basic = (/ -1,  -.8,  -.6, -.4, -.2, 0., .2, .4, .6, .8, 1./)

;;diag_script_info@cn_levels_stddev_basic = (/0.0, 0.2, .4, .6, .8, 1.2, 1.4, 1.6, 1.8, 2.0/)
;;diag_script_info@cn_levels_stddev_diff_basic = (/-2, -1.8, -1.6, -1.4, -1.2, -1, -.8, -.6, -.4, -.2, 0, .2, .4, .6,.8, 1,1.2,1.4,1.6,1.8, 2/)
;;diag_script_info@plottype_lbres = False
