diag_script_info = True

;; Used by: clouds_interannual.ncl
diag_script_info@projection = "CylindricalEquidistant"   ; map projection, e.g., Mollweide, Mercator
diag_script_info@colormap = "WhiteBlueGreenYellowRed"  ; e.g., WhiteBlueGreenYellowRed, rainbow

; --- multi-model mean parameters -----------------------------
; calculate multi-model mean?
; valid options: True, False
diag_script_info@calcmm = True

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest", "1x1")
; Note: if "ref" is selected, the reference model has to be specified via
; diag_script_info@ref_model
diag_script_info@grid = "1x1"
; specify which model is the reference model in case of grid = "ref"
diag_script_info@ref_model = "n.a."

; write plots for individual models to separate files? (True, False)
diag_script_info@extrafiles = False

; data sets to be skipped when calculating the multi-model mean
; (= observations)
diag_script_info@models_to_skip = (/"UWisc", "CERES-EBAF", "MODIS", "MODIS-CFMIP", "GPCP-SG", \
             "ERA-Interim", "PATMOS", "CLARA-A2", "ESACCI-CLOUD"/)
; -------------------------------------------------------------

