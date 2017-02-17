diag_script_info = True

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest", "1x1")

diag_script_info@grid = "1x1"
diag_script_info@ref_model = "nan"  ; needed only when grid="ref"

diag_script_info@stddev = True
diag_script_info@lsmask = 0  ; 0 = no masking, 1 = mask land grid cells, 2 = mask ocean grid cells

