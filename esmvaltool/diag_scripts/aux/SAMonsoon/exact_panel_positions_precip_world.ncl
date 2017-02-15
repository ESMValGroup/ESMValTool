;;#############################################################################
;; SPECIFIC ROUTINES FOR EXPLICIT PLACEMENT OF
;; GLOBAL PRECIP MULTIPANEL PLOTS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;
;;#############################################################################

;;########################################################################
;;######                                                            ######
;;######        EXPLICIT PANEL POSITIONS FOR N BY ONE PLOTS         ######
;;######                                                            ######
;;########################################################################
n_by_cols_1_first_vpXF = 0.176483
n_by_cols_1_first_vpYF = 0.9
n_by_cols_1_first_vpWidthF = 0.667565
n_by_cols_1_first_vpHeightF = 0.185435

n_by_cols_1_ver_space = 0.064565
n_by_cols_1_hor_space = 0

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
;;######       EXPLICIT PANEL POSITIONS FOR N BY TWO PLOTS          ######
;;######                                                            ######
;;########################################################################
n_by_cols_2_first_vpXF = 0.0237663
n_by_cols_2_first_vpYF = 0.9
n_by_cols_2_first_vpWidthF = 0.466824
n_by_cols_2_first_vpHeightF = 0.129673

n_by_cols_2_ver_space = 0.05515
n_by_cols_2_hor_space = 0.03318

n_by_cols_2_lb_offset = 0.0

;;########################################################################
;;######                                                            ######
;;######   VARIABLES BELOW ARE NOT USED (BUT NEEDS TO BE DEFINED)   ######
;;######                                                            ######
;;########################################################################
three_by_one_first_vpXF = 0.0
three_by_one_first_vpYF = 0.0
three_by_one_first_vpWidthF = 0.0
three_by_one_first_vpHeightF = 0.0

three_by_one_ver_space = 0.0

; A-laue_ax+
n_by_cols_3_lb_offset = 0.0
; A-laue_ax-

;;########################################################################
;;######                                                            ######
;;######   VARIABLES BELOW ARE NOT USED (BUT NEEDS TO BE DEFINED)   ######
;;######                                                            ######
;;########################################################################
n_by_cols_3_first_vpXF = 0.0
n_by_cols_3_first_vpYF = 0.0
n_by_cols_3_first_vpWidthF = 0.0
n_by_cols_3_first_vpHeightF = 0.0

n_by_cols_3_ver_space = 0.0
n_by_cols_3_hor_space = 0.0

max_cols = 2  ; Max number of columns on multipanel page
max_figures_pp = 8  ; pp = per_page
min_figures_pp = 4  ; pp = per_page
