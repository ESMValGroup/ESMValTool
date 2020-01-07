; #############################################################################
; GENERAL ROUTINES FOR MODEL ENSEMBLE OPERATIONS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function get_start_year
;    function get_end_year
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/data_handling.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

; #############################################################################
undef("get_start_year")
function get_start_year(dataset_index[1]:numeric)
;
; Arguments:
;    dataset_index: index of array "datasets"; this function is supposed to be
;                 called inside a dataset loop
;
; Return value
;    An integer containing the first year
;
; Description
;    Chooses first year of dataset input with respect to
;    diag_script_info@range_option
;
; Caveats
;
; References
;
; Modification history
;    20140128-senftleben_daniel: written.
;
local funcname, scriptname, dataset_index
begin

  funcname = "get_start_year"
  scriptname = "diag_scripts/shared/ensemble.ncl"
  enter_msg(scriptname, funcname)

  if (.not.isatt(diag_script_info, "range_option")) then
    error_msg("f", "ensemble.ncl", funcname, "diag_script_info must " + \
              "contain the range_option attribute")
  end if
  if (diag_script_info@range_option.eq.0) then
    start_year = toint(input_file_info@start_year(dataset_index))
  else if (diag_script_info@range_option.eq.1) then
    start_year = max(toint(input_file_info@start_year))
  else
    error_msg("f", "ensemble.ncl", funcname, "range_option not " + \
              "properly set, must be either 0 or 1")
  end if
  end if

  leave_msg(scriptname, funcname)
  return(start_year)

end

; #############################################################################
undef("get_end_year")
function get_end_year(dataset_index:numeric)
;
; Arguments
;    dataset_index: index of array "datasets"; this function is supposed to be
;                 called inside a dataset loop
;
; Return value
;    An integer containing the last year
;
; Description
;    Chooses last year of dataset input with respect to
;    diag_script_info@range_option
;
; Caveats
;
; References
;
; Modification history
;    20140128-senftleben_daniel: written.
;
local funcname, scriptname, dataset_index
begin

  funcname = "get_end_year"
  scriptname = "diag_scripts/shared/ensemble.ncl"
  enter_msg(scriptname, funcname)

  if (.not.isatt(diag_script_info, "range_option")) then
    error_msg("f", "ensemble.ncl", funcname, "diag_script_info must " + \
              "contain the range_option attribute")
  end if
  if (diag_script_info@range_option.eq.0) then
    end_year = toint(input_file_info@end_year(dataset_index))
  else if (diag_script_info@range_option.eq.1) then
    end_year = min(toint(input_file_info@end_year))
  else
    error_msg("f", "ensemble.ncl", funcname, "range_option not " + \
              "properly set, must be either 0 or 1")
  end if
  end if

  leave_msg(scriptname, funcname)
  return(end_year)

end
