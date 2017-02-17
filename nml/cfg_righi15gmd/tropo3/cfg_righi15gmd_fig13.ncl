diag_script_info = True

;; Used by: ancyc_lat.ncl.ncl
diag_script_info@styleset = "righi15gmd"
diag_script_info@font = 21
diag_script_info@range_option = 0  ; 0 = as in nml, 1 = overlapping time period
diag_script_info@lbLabelBarOn = True  ; Plot a Label Bar (True/False)?
diag_script_info@mpProjection = "CylindricalEquidistant"
diag_script_info@gsnZonalMean = False  ; plot zonal mean beside map 

; Options added to create fig. 24 of Eyring et al. 2016
; ----------------------------------------------------
; Calculate multi-model mean? (valid options: True, False)
diag_script_info@calcmm = True
; Calculate difference to reference dataset? (valid options: True, False)
; Note: the name of reference data set is given via namelist (ref_model = "..."
; in the <variable>-tag).
diag_script_info@calcdiff = True
