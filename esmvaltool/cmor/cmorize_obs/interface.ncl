; #############################################################################
; INTERFACE TO HANDLE THE COMMUNICATION BETWEEN THE PYTHON WORKFLOW AND NCL
; #############################################################################
; Load the automatically-generated settings.ncl for the current dataset,
; utility and logging functions.
;
; No functions/procedures shall be further added to this script.
; #############################################################################

; Load interface settings
loadscript("$settings")

; Load logging functions
loadscript(getenv("esmvaltool_root") + "/interface_scripts/logging.ncl")

; Load utility functions
loadscript(getenv("esmvaltool_root") + "/cmor/cmorize_obs/cmorize_obs_utils.ncl")

