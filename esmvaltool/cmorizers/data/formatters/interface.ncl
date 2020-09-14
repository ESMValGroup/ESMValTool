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
loadscript(getenv("esmvaltool_root") + \
           "/../interface_scripts/logging.ncl")

; Load utility functions
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/utilities.ncl")

; Check trailing slash
if (str_get_cols(input_dir_path, -1, -1).ne."/") then
  input_dir_path = input_dir_path + "/"
end if
if (str_get_cols(output_dir_path, -1, -1).ne."/") then
  output_dir_path = output_dir_path + "/"
end if
