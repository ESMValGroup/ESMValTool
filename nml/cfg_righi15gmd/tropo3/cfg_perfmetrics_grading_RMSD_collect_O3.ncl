diag_script_info = True

;; Used by: perfmetrics_grading_collect.ncl

;; Set labelbar values
diag_script_info@label_bounds = (/0, 2.75/)
diag_script_info@label_scale = 0.25

;; Add upper/lower triangle for values outside range
diag_script_info@label_lo = False

;; Set colormap
diag_script_info@colormap = "WhiteYellowOrangeRed"
diag_script_info@cm_interval = (/18, 253/)

;; Show grading values on the plot
diag_script_info@disp_values = False

;; Plot options
diag_script_info@title = "Ozone columns RMSD metrics"
diag_script_info@scale_font = 1.4
