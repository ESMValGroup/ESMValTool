;;#############################################################################
;; SPECIFIC ROUTINES FOR EXPLICIT PLACEMENT OF
;; PROPAGATION PRECIP MULTIPANEL PLOTS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;
;;#############################################################################

;;########################################################################
;;######                                                            ######
;;######       EXPLICIT PANEL POSITIONS FOR TWO BY ONE PLOTS        ######
;;######                                                            ######
;;########################################################################
two_by_one_first_vpXF = 0.270517
two_by_one_first_vpYF = 0.923324
two_by_one_first_vpWidthF = 0.48
two_by_one_first_vpHeightF = 0.35
two_by_one_ver_space = 0.08098

;;########################################################################
;;######                                                            ######
;;######      EXPLICIT PANEL POSITIONS FOR THREE BY ONE PLOTS       ######
;;######                                                            ######
;;########################################################################
three_by_one_first_vpXF = 0.270517
three_by_one_first_vpYF = 0.943324
three_by_one_first_vpWidthF = 0.48
three_by_one_first_vpHeightF = 0.24

three_by_one_ver_space = 0.07098

;;########################################################################
;;######                                                            ######
;;######       EXPLICIT PANEL POSITIONS FOR N BY ONE PLOTS          ######
;;######                                                            ######
;;########################################################################
n_by_cols_1_first_vpXF = 0.27
n_by_cols_1_first_vpYF = 0.9
n_by_cols_1_first_vpWidthF = 0.44
n_by_cols_1_first_vpHeightF = 0.21

n_by_cols_1_ver_space = 0.08
n_by_cols_1_hor_space = 0.00

n_by_cols_1_tiMainFontHeightF = 0.013
n_by_cols_1_tiMainOffsetYF = -0.015

;;########################################################################
;;######                                                            ######
;;######       EXPLICIT PANEL POSITIONS FOR N BY TWO PLOTS          ######
;;######                                                            ######
;;########################################################################
n_by_cols_2_first_vpXF = 0.08
n_by_cols_2_first_vpYF = 0.910164
n_by_cols_2_first_vpWidthF = 0.42
n_by_cols_2_first_vpHeightF = 0.2

n_by_cols_2_ver_space = 0.062756
n_by_cols_2_hor_space = 0.0500824

n_by_cols_2_tiMainFontHeightF = 0.012
n_by_cols_2_tiMainOffsetYF = -0.013
n_by_cols_2_lb_offset = 0.015

;;########################################################################
;;######                                                            ######
;;######      EXPLICIT PANEL POSITIONS FOR N BY THREE PLOTS         ######
;;######                                                            ######
;;########################################################################
n_by_cols_3_first_vpXF = 0.050
n_by_cols_3_first_vpYF = 0.866538
n_by_cols_3_first_vpWidthF = 0.297685
n_by_cols_3_first_vpHeightF = 0.148842

n_by_cols_3_ver_space = 0.051825
n_by_cols_3_hor_space = 0.0166486

n_by_cols_3_tiMainFontHeightF = 0.010
n_by_cols_3_tmFontHeightF = 0.008
n_by_cols_3_tiMainOffsetYF = -0.008
n_by_cols_3_lb_offset = 0.000

; How to distribute figures on pages
max_cols = 3  ; Max number of columns on multipanel page
max_figures_pp = 12  ; pp = per_page
min_figures_pp = 4   ; pp = per_page

; Manually adjust label bar
three_by_one_lb_amres_adjust = -0.16
two_by_one_lb_amres_adjust = 0.1
n_by_cols_lb_amres_adjust = 0.1
