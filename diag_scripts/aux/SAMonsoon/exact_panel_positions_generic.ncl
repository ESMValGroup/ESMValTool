;;#############################################################################
;; GENERAL ROUTINES FOR EXPLICIT PLACEMENT OF MULTIPANEL PLOTS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;
;;#############################################################################
begin
    loadscript(exact_panel_position_file)
end

;;########################################################################
;;######                                                            ######
;;######        RELATIVE PANEL POSITIONS FOR TWO BY ONE PLOTS       ######
;;######                                                            ######
;;########################################################################
;; Second plot
two_by_one_second_vpXF = two_by_one_first_vpXF
two_by_one_second_vpYF = two_by_one_first_vpYF - two_by_one_first_vpHeightF - two_by_one_ver_space
two_by_one_second_vpWidthF = two_by_one_first_vpWidthF
two_by_one_second_vpHeightF = two_by_one_first_vpHeightF

;;########################################################################
;;######                                                            ######
;;######       RELATIVE PANEL POSITIONS FOR THREE BY ONE PLOTS      ######
;;######                                                            ######
;;########################################################################
;; Second plot
three_by_one_second_vpXF = three_by_one_first_vpXF
three_by_one_second_vpYF = three_by_one_first_vpYF - three_by_one_first_vpHeightF - three_by_one_ver_space
three_by_one_second_vpWidthF = three_by_one_first_vpWidthF
three_by_one_second_vpHeightF = three_by_one_first_vpHeightF

;; Third plot
three_by_one_third_vpXF = three_by_one_second_vpXF
three_by_one_third_vpYF = three_by_one_second_vpYF - three_by_one_second_vpHeightF - three_by_one_ver_space
three_by_one_third_vpWidthF = three_by_one_second_vpWidthF
three_by_one_third_vpHeightF = three_by_one_second_vpHeightF

;;########################################################################
;;######                                                            ######
;;######        RELATIVE PANEL POSITIONS FOR N BY ONE PLOTS         ######
;;######                                                            ######
;;########################################################################
n_by_cols_1_second_vpXF      = n_by_cols_1_first_vpXF
n_by_cols_1_second_vpYF      = n_by_cols_1_first_vpYF - n_by_cols_1_first_vpHeightF - n_by_cols_1_ver_space
n_by_cols_1_second_vpWidthF  = n_by_cols_1_first_vpWidthF
n_by_cols_1_second_vpHeightF = n_by_cols_1_first_vpHeightF

n_by_cols_1_third_vpXF       = n_by_cols_1_second_vpXF
n_by_cols_1_third_vpYF       = n_by_cols_1_second_vpYF - n_by_cols_1_second_vpHeightF - n_by_cols_1_ver_space
n_by_cols_1_third_vpWidthF   = n_by_cols_1_second_vpWidthF
n_by_cols_1_third_vpHeightF  = n_by_cols_1_second_vpHeightF

n_by_cols_1_fourth_vpXF      = n_by_cols_1_third_vpXF
n_by_cols_1_fourth_vpYF      = n_by_cols_1_third_vpYF - n_by_cols_1_third_vpHeightF - n_by_cols_1_ver_space
n_by_cols_1_fourth_vpWidthF  = n_by_cols_1_third_vpWidthF
n_by_cols_1_fourth_vpHeightF = n_by_cols_1_third_vpHeightF

;;########################################################################
;;######                                                            ######
;;######        RELATIVE PANEL POSITIONS FOR N BY TWO PLOTS         ######
;;######                                                            ######
;;########################################################################
n_by_cols_2_second_vpXF      = n_by_cols_2_first_vpXF + n_by_cols_2_first_vpWidthF + n_by_cols_2_hor_space
n_by_cols_2_second_vpYF      = n_by_cols_2_first_vpYF
n_by_cols_2_second_vpWidthF  = n_by_cols_2_first_vpWidthF
n_by_cols_2_second_vpHeightF = n_by_cols_2_first_vpHeightF

n_by_cols_2_third_vpXF       = n_by_cols_2_first_vpXF
n_by_cols_2_third_vpYF       = n_by_cols_2_first_vpYF - n_by_cols_2_first_vpHeightF - n_by_cols_2_ver_space
n_by_cols_2_third_vpWidthF   = n_by_cols_2_first_vpWidthF
n_by_cols_2_third_vpHeightF  = n_by_cols_2_first_vpHeightF

n_by_cols_2_fourth_vpXF      = n_by_cols_2_third_vpXF + n_by_cols_2_third_vpWidthF + n_by_cols_2_hor_space
n_by_cols_2_fourth_vpYF      = n_by_cols_2_third_vpYF
n_by_cols_2_fourth_vpWidthF  = n_by_cols_2_first_vpWidthF
n_by_cols_2_fourth_vpHeightF = n_by_cols_2_first_vpHeightF

n_by_cols_2_fifth_vpXF       = n_by_cols_2_third_vpXF
n_by_cols_2_fifth_vpYF       = n_by_cols_2_third_vpYF - n_by_cols_2_third_vpHeightF - n_by_cols_2_ver_space
n_by_cols_2_fifth_vpWidthF   = n_by_cols_2_first_vpWidthF
n_by_cols_2_fifth_vpHeightF  = n_by_cols_2_first_vpHeightF

n_by_cols_2_sixth_vpXF       = n_by_cols_2_fifth_vpXF + n_by_cols_2_fifth_vpWidthF + n_by_cols_2_hor_space
n_by_cols_2_sixth_vpYF       = n_by_cols_2_fifth_vpYF
n_by_cols_2_sixth_vpWidthF   = n_by_cols_2_first_vpWidthF
n_by_cols_2_sixth_vpHeightF  = n_by_cols_2_first_vpHeightF

n_by_cols_2_seventh_vpXF      = n_by_cols_2_fifth_vpXF
n_by_cols_2_seventh_vpYF      = n_by_cols_2_fifth_vpYF - n_by_cols_2_fifth_vpHeightF - n_by_cols_2_ver_space
n_by_cols_2_seventh_vpWidthF  = n_by_cols_2_first_vpWidthF
n_by_cols_2_seventh_vpHeightF = n_by_cols_2_first_vpHeightF

n_by_cols_2_eighth_vpXF       = n_by_cols_2_seventh_vpXF + n_by_cols_2_seventh_vpWidthF + n_by_cols_2_hor_space
n_by_cols_2_eighth_vpYF       = n_by_cols_2_seventh_vpYF
n_by_cols_2_eighth_vpWidthF   = n_by_cols_2_first_vpWidthF
n_by_cols_2_eighth_vpHeightF  = n_by_cols_2_first_vpHeightF

n_by_cols_2_ninth_vpXF      = n_by_cols_2_fifth_vpXF
n_by_cols_2_ninth_vpYF      = n_by_cols_2_seventh_vpYF - n_by_cols_2_seventh_vpHeightF - n_by_cols_2_ver_space
n_by_cols_2_ninth_vpWidthF  = n_by_cols_2_first_vpWidthF
n_by_cols_2_ninth_vpHeightF = n_by_cols_2_first_vpHeightF

n_by_cols_2_tenth_vpXF       = n_by_cols_2_ninth_vpXF + n_by_cols_2_ninth_vpWidthF + n_by_cols_2_hor_space
n_by_cols_2_tenth_vpYF       = n_by_cols_2_ninth_vpYF
n_by_cols_2_tenth_vpWidthF   = n_by_cols_2_first_vpWidthF
n_by_cols_2_tenth_vpHeightF  = n_by_cols_2_first_vpHeightF

;;########################################################################
;;######                                                            ######
;;######       RELATIVE PANEL POSITIONS FOR N BY THREE PLOTS        ######
;;######                                                            ######
;;########################################################################
n_by_cols_3_second_vpXF      = n_by_cols_3_first_vpXF + n_by_cols_3_first_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_second_vpYF      = n_by_cols_3_first_vpYF
n_by_cols_3_second_vpWidthF  = n_by_cols_3_first_vpWidthF
n_by_cols_3_second_vpHeightF = n_by_cols_3_first_vpHeightF

n_by_cols_3_third_vpXF       = n_by_cols_3_second_vpXF + n_by_cols_3_second_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_third_vpYF       = n_by_cols_3_first_vpYF
n_by_cols_3_third_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_third_vpHeightF  = n_by_cols_3_first_vpHeightF

n_by_cols_3_fourth_vpXF      = n_by_cols_3_first_vpXF
n_by_cols_3_fourth_vpYF      = n_by_cols_3_first_vpYF - n_by_cols_3_first_vpHeightF - n_by_cols_3_ver_space
n_by_cols_3_fourth_vpWidthF  = n_by_cols_3_first_vpWidthF
n_by_cols_3_fourth_vpHeightF = n_by_cols_3_first_vpHeightF

n_by_cols_3_fifth_vpXF       = n_by_cols_3_fourth_vpXF + n_by_cols_3_fourth_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_fifth_vpYF       = n_by_cols_3_fourth_vpYF
n_by_cols_3_fifth_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_fifth_vpHeightF  = n_by_cols_3_first_vpHeightF

n_by_cols_3_sixth_vpXF       = n_by_cols_3_fifth_vpXF + n_by_cols_3_fifth_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_sixth_vpYF       = n_by_cols_3_fourth_vpYF
n_by_cols_3_sixth_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_sixth_vpHeightF  = n_by_cols_3_first_vpHeightF

n_by_cols_3_seventh_vpXF      = n_by_cols_3_first_vpXF
n_by_cols_3_seventh_vpYF      = n_by_cols_3_fourth_vpYF - n_by_cols_3_fourth_vpHeightF - n_by_cols_3_ver_space
n_by_cols_3_seventh_vpWidthF  = n_by_cols_3_first_vpWidthF
n_by_cols_3_seventh_vpHeightF = n_by_cols_3_first_vpHeightF

n_by_cols_3_eigth_vpXF       = n_by_cols_3_seventh_vpXF + n_by_cols_3_seventh_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_eigth_vpYF       = n_by_cols_3_seventh_vpYF
n_by_cols_3_eigth_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_eigth_vpHeightF  = n_by_cols_3_first_vpHeightF

n_by_cols_3_ninth_vpXF       = n_by_cols_3_eigth_vpXF + n_by_cols_3_eigth_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_ninth_vpYF       = n_by_cols_3_seventh_vpYF
n_by_cols_3_ninth_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_ninth_vpHeightF  = n_by_cols_3_first_vpHeightF

n_by_cols_3_tenth_vpXF      = n_by_cols_3_first_vpXF
n_by_cols_3_tenth_vpYF      = n_by_cols_3_seventh_vpYF - n_by_cols_3_seventh_vpHeightF - n_by_cols_3_ver_space
n_by_cols_3_tenth_vpWidthF  = n_by_cols_3_first_vpWidthF
n_by_cols_3_tenth_vpHeightF = n_by_cols_3_first_vpHeightF

n_by_cols_3_eleventh_vpXF       = n_by_cols_3_tenth_vpXF + n_by_cols_3_tenth_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_eleventh_vpYF       = n_by_cols_3_tenth_vpYF
n_by_cols_3_eleventh_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_eleventh_vpHeightF  = n_by_cols_3_first_vpHeightF

n_by_cols_3_twelveth_vpXF       = n_by_cols_3_eleventh_vpXF + n_by_cols_3_eleventh_vpWidthF + n_by_cols_3_hor_space
n_by_cols_3_twelveth_vpYF       = n_by_cols_3_tenth_vpYF
n_by_cols_3_twelveth_vpWidthF   = n_by_cols_3_first_vpWidthF
n_by_cols_3_twelveth_vpHeightF  = n_by_cols_3_first_vpHeightF
