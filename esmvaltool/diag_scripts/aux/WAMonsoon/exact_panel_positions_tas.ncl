;;########################################################################
;;######                                                            ######
;;######      EXPLICIT PANEL POSITIONS FOR THREE BY ONE PLOTS       ######
;;######                                                            ######
;;########################################################################
three_by_one_first_vpXF = 0.29
three_by_one_first_vpYF = 0.94
three_by_one_first_vpWidthF = 0.33
three_by_one_first_vpHeightF = 0.24

three_by_one_ver_space = 0.07

;;########################################################################
;;######                                                            ######
;;######   DUMMY VALUES FOR PANEL POSITIONS FOR TWO BY ONE PLOTS    ######
;;######                                                            ######
;;########################################################################
two_by_one_first_vpXF      = default_fillvalue("float")
two_by_one_first_vpYF      = default_fillvalue("float")
two_by_one_first_vpWidthF  = default_fillvalue("float")
two_by_one_first_vpHeightF = default_fillvalue("float")
two_by_one_ver_space       = default_fillvalue("float")

two_by_one_first_vpXF@_FillValue      = default_fillvalue("float")
two_by_one_first_vpYF@_FillValue      = default_fillvalue("float")
two_by_one_first_vpWidthF@_FillValue  = default_fillvalue("float")
two_by_one_first_vpHeightF@_FillValue = default_fillvalue("float")
two_by_one_ver_space@_FillValue       = default_fillvalue("float")

;;########################################################################
;;######                                                            ######
;;######        EXPLICIT PANEL POSITIONS FOR N BY ONE PLOTS         ######
;;######                                                            ######
;;########################################################################
n_by_cols_1_first_vpXF = 0.27
n_by_cols_1_first_vpYF = 0.9
n_by_cols_1_first_vpWidthF = 0.40
n_by_cols_1_first_vpHeightF = 0.19

n_by_cols_1_ver_space = 0.07
n_by_cols_1_hor_space = 0.00

;;########################################################################
;;######                                                            ######
;;######        EXPLICIT PANEL POSITIONS FOR N BY TWO PLOTS         ######
;;######                                                            ######
;;########################################################################
n_by_cols_2_first_vpXF = 0.20
n_by_cols_2_first_vpYF = 0.91
n_by_cols_2_first_vpWidthF = 0.30
n_by_cols_2_first_vpHeightF = 0.26

n_by_cols_2_ver_space = 0.03
n_by_cols_2_hor_space = 0.02

n_by_cols_2_lb_offset = -0.02

;;########################################################################
;;######                                                            ######
;;######       EXPLICIT PANEL POSITIONS FOR N BY THREE PLOTS        ######
;;######                                                            ######
;;########################################################################
n_by_cols_3_first_vpXF = 0.039
n_by_cols_3_first_vpYF = 0.91
n_by_cols_3_first_vpWidthF = 0.28
n_by_cols_3_first_vpHeightF = 0.26

n_by_cols_3_ver_space = 0.020
n_by_cols_3_hor_space = 0.027

n_by_cols_3_lb_offset = 0.00

; How to distribute figures on pages
max_cols = 3  ; Max number of columns on multipanel page
max_figures_pp = 9  ; pp = per_page
min_figures_pp = 4   ; pp = per_page

; Manually adjust label bar
three_by_one_lb_amres_adjust = -0.16
