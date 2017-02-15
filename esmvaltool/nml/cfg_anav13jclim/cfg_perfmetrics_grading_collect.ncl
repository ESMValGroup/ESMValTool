
diag_script_info = True

;; Used by: perfmetrics_grading_collect.ncl

;; Set labelbar values
diag_script_info@label_bounds = (/0.0, 1.0/)
diag_script_info@label_scale = 0.1

;; Add upper/lower triangle for values outside range
diag_script_info@label_lo = False
diag_script_info@label_hi = False

;; Set colormap
diag_script_info@colormap = "matlab_jet"
diag_script_info@cm_interval = (/2, 63/)
diag_script_info@cm_reverse = True

;; Sort model in alphabetic order (excluding MMM)
diag_script_info@sort = True

;; Display metric values
diag_script_info@disp_values = False

;; Diplay model rankings in decreasing order (higher first)
diag_script_info@disp_rankings = True
diag_script_info@rank_order = -1
