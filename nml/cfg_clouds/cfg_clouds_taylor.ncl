diag_script_info = True

; -----------------------------------------------------------------------------
; used by clouds_taylor.ncl
;diag_script_info@styleset = "CMIP5"         ; "CMIP5" (if not set, plotting routine will
                                             ; create symbols and select colors)
; -----------------------------------------------------------------------------
; calculate multi-model mean (True, False)?
diag_script_info@calcmm = True
; -----------------------------------------------------------------------------
; --- time averaging ---
; timemean = "annualclim": create 1 plot (annual mean)
;          = "seasonalclim": create 4 separate plots (1 for each season - DJF, MAM, JJA, SON)
diag_script_info@timemean = "annualclim"
; -----------------------------------------------------------------------------
; target grid used for regridding
; grid = "ref": use grid of reference model (defined below)
;      = "coarsest": use coarsest grid of all models/observations
;      = "finest": use highest resolution grid of all models/observations
;      = "1x1": regrid all datasets to 1 deg x 1 deg grid
diag_script_info@grid = "ref"
; -----------------------------------------------------------------------------
; minimum fraction of valid data points (for creating a mask)
diag_script_info@valid_fraction = 0.5
; -----------------------------------------------------------------------------
; embracelegend = True:  1) save legend to separate file
;                        2) plot alternative reference data set (if available) as a black or red star
;                           and use "altern. ref. dataset" as label
;               = False: 1) plot legend into taylor plot (same file)
;                        2) plot alternative reference data set (if available) using model symbols
;                           and use name of observations as label
diag_script_info@embracelegend = False
; -----------------------------------------------------------------------------
; definitions of observations used as reference
if (variables(0).eq."lwp") then
    diag_script_info@ref_model = (/"UWisc"/)
end if
if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut","SW_CRE","LW_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
end if
if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"MODIS"/)
end if
if (any(variables(0).eq.(/"clivi"/))) then
    diag_script_info@ref_model = (/"MODIS-CFMIP"/)
;    diag_script_info@ref_model = (/"MODIS-L3"/)
end if
if (any(variables(0).eq.(/"pr","pr-mmday"/))) then
    diag_script_info@ref_model = (/"GPCP-SG"/)
end if
if (any(variables(0).eq.(/"prw"/))) then
    diag_script_info@ref_model = (/"SSMI-MERIS"/)
end if

