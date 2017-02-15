;;#############################################################################
;; SPECIFIC ROUTINES FOR EXPLICIT PLACEMENT OF
;; INTRASEASONAL PRECIP MULTIPANEL PLOTS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;
;;#############################################################################

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
;;######      EXPLICIT PANEL POSITIONS FOR THREE BY ONE PLOTS       ######
;;######                                                            ######
;;########################################################################
three_by_one_first_vpXF = 0.170517
three_by_one_first_vpYF = 0.923324
three_by_one_first_vpWidthF = 0.68
three_by_one_first_vpHeightF = 0.23

three_by_one_ver_space = 0.075

;;########################################################################
;;######                                                            ######
;;######       EXPLICIT PANEL POSITIONS FOR N BY TWO PLOTS          ######
;;######                                                            ######
;;########################################################################
n_by_cols_2_first_vpXF = 0.0376816
n_by_cols_2_first_vpYF = 0.910164
n_by_cols_2_first_vpWidthF = 0.45
n_by_cols_2_first_vpHeightF = 0.15

n_by_cols_2_ver_space = 0.062756
n_by_cols_2_hor_space = 0.0500824

n_by_cols_2_lb_offset = 0.0

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
n_by_cols_3_hor_space = 0.0356486
n_by_cols_3_hor_space = 0.0166486

; A-laue_ax+
n_by_cols_3_lb_offset = 0.0
; A-laue_ax-

;;########################################################################
;;######                                                            ######
;;######       EXPLICIT PANEL POSITIONS FOR N BY ONE PLOTS          ######
;;######                                                            ######
;;########################################################################
n_by_cols_1_first_vpXF = 0.17
n_by_cols_1_first_vpYF = 0.9
n_by_cols_1_first_vpWidthF = 0.66
n_by_cols_1_first_vpHeightF = 0.25

n_by_cols_1_ver_space = 0.03
n_by_cols_1_hor_space = 0.00

; How to distribute figures on pages
max_cols = 2  ; Max number of columns on multipanel page
max_figures_pp = 8  ; pp = per_page
min_figures_pp = 3   ; pp = per_page

; Manually adjust label bar
three_by_one_lb_amres_adjust = -0.16
