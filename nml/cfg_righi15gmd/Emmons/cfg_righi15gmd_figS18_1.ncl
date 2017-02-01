diag_script_info = True
;; Used by: Emmons.ncl

;; ############################################################################
;; OBS data
;;   All files matching the following will be scanned, all suitable data will
;;   be used (mandatory). Wildcards are allowed.
;; ############################################################################
path = getenv("ESMValTool_OBSPATH") + "/Tier2/Emmons/"
diag_script_info@obsfiles = (/\
    path + "*PEM-Tropics-A-DC8-Hawaii_T0M_vmrc3h8.nc"/)
delete(path)

;; ############################################################################
;; OBS parameters
;;   The target variable (e.g. "vmrc3h8") + the following extension gives
;;   the variable name in the obs file.
;;   Deactivate lines corresponding to parameters not in the obs file,
;;   but either obsvar_mean or obsvar_median is mandatory.
;;   "_N" is only used for grading: Will use equal weights if disabled.
;;   Disabling any other quantity will just cause the corresponding whisker
;;   not to be drawn. This may be used to switch off whiskers.
;; ############################################################################
diag_script_info@obsvar_N = "_N"  ; number of observations
diag_script_info@obsvar_mean = "_mean" 	; mean
diag_script_info@obsvar_stddev = "_stddev"  ; standard deviation
diag_script_info@obsvar_min = "_min"  ; min
diag_script_info@obsvar_max = "_max"  ; max
diag_script_info@obsvar_5 = "_5"  ; 5% percentile
diag_script_info@obsvar_25 = "_25"  ; 25% percentile
diag_script_info@obsvar_median = "_median"  ; 50% percentile
diag_script_info@obsvar_75 = "_75"  ; 75% percentile
diag_script_info@obsvar_95 = "_95"  ; 95% percentile

;; ############################################################################
;; Campaign parameters
;;   Expected to be global parameters of each obs file (all mandatory)
;; ############################################################################
diag_script_info@obsatt_campaign = "campaign"  ; for plot annotation
diag_script_info@obsatt_location = "location"  ; for plot annotation
diag_script_info@obsatt_period = "period"  ; for matching time of year ...
diag_script_info@obsatt_latmin = "latmin"  ; for matching region ...
diag_script_info@obsatt_latmax = "latmax"
diag_script_info@obsatt_lonmin = "lonmin"  ; ... in the model data
diag_script_info@obsatt_lonmax = "lonmax"

;; ############################################################################
;; GRADING parameters
;;   (enable to activate grading -> not implemented yet)
;; ############################################################################
;diag_script_info@grading_file = "./work/Emmons_grading.ncl"
    ;; target file for grading values, Default: no grading

;; ############################################################################
;; FILTER parameters for selection of obs data
;;   (none of them mandatory, because there are defaults)
;; ############################################################################
diag_script_info@roi = (/-90., 90., 0., 360./)  ; region of interest:
    ;; (/latmin,latmax,lonmin,lonmax/), Default: global=(/-90.,90.,0.,360./)
diag_script_info@roi_match = "center"  ; region defined by obsatt vs roi:
    ;; select from "contained" / "center" / "overlap" / "encloses" / "outside",
    ;; Default: "center"
diag_script_info@poi = (/1992, 2000/)  ; period of interest for obs:
    ;; (/first, last/) years to be considered, Default: (/1900, 2100/)
diag_script_info@poi_match = "overlap"  ; period defined by obsatt vs poi:
    ;; select from "contained" / "center" / "overlap" / "encloses" / "outside",
    ;; Default: "overlap"

;; ############################################################################
;; GENERAL parameters
;;   (none of them mandatory, because there are defaults)
;; ############################################################################
diag_script_info@quantity = "mean"  ; determines quantities to be evaluated and
    ;; plotted for grading
    ;; Default: "mean" -> N, mean, stddev, min,max
    ;; "median" (not fully implemented yet) -> N, median, 5, 25, 75, 95
;diag_script_info@ncdf = "default"  ; enable to output to netCDF:
    ;; either use "default" or give a full file name, Default: no netCDF output
diag_script_info@interpolation = -2  ; flag determining regridding from
    ;; simulation to observations' pressure levels:
    ;; abs: 1 linear, !=1 log; sign: -/+ with/without extrapolation
    ;; Default: -1 or +1 will be used, depending on levels

;; ############################################################################
;; PANELLING parameters
;;   (none of them mandatory, because there are defaults)
;; ############################################################################
diag_script_info@max_vert = 1  ; plots per column, Default: 1
diag_script_info@max_hori = 1  ; plots per row, Default: 1
diag_script_info@aux_info = "righi15gmd_figS18_1" 
    ;; string to include in outfile name, Default: period of interest 

;; ############################################################################
;; STYLE parameters (used by function profile_plev)
;;   (none of them mandatory, because there are defaults)
;; ############################################################################
diag_script_info@styleset = "righi15gmd"  ; see diag_scripts/lib/ncl/style.ncl
    ;; for available styles, Default: "DEFAULT"
;diag_script_info@colors = (/"black"/)  ; The following entries override ...
;diag_script_info@dashes = (/"0"/)  ; ... styleset settings
diag_script_info@thicks = (/2, 1, 1, 1, 1, 1, 1/) * 6.
;diag_script_info@annots = (/"data"/)
;diag_script_info@avgstd = (/"0"/)
diag_script_info@font = 21  ; overrides default font
    ;; (see www.ncl.ucar.edu/Document/Graphics/Resources/tx.shtml#txFont)

;; ############################################################################
;; PLOTTING parameters (used by function profile_plev)
;;   (none of them mandatory, because there are defaults)
;; ############################################################################
diag_script_info@ptop = 193.3   ; top pressure:
    ;; "auto" = observational range per file
    ;; otherwise give p in [hPa], Default: minimum of obs pressure coordinate
    ;; Emmons plev [hPa]: 1013.25, 898.746, 794.952, 701.085, 616.402, 540.199,
    ;;                    471.81, 410.607, 355.998, 307.425, 264.363, 226.321,
    ;;                    226.321, 193.304, 165.104, 141.018, 120.446
;diag_script_info@pbot = "auto"   ; surface pressure:
    ;; "auto" = observational range per file; otherwise give p in [hPa]
    ;; Default: use maximum of obs pressure coordinate
diag_script_info@zoom = (/0, 150/)  ; zoom x axis range:
    ;; "yes" = consider only mean/median OR also min/max/5%/95% ("no"),
    ;; OR specify a range, e.g. (/0, 300/), Default: "yes"
;diag_script_info@inline_legend = True  ; Legend in the plot?, Default = False
