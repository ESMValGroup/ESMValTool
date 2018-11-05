; #############################################################################
; INTERFACE TO HANDLE THE COMMUNICATION BETWEEN THE PYTHON WORKFLOW AND NCL
; #############################################################################
; Load the automatically-generated interface files <variable>_info.ncl and
; settings.ncl for the current diagnostic and variable with the current
; preprocessor set. These files contains the NCL logical _info varialbe with
; all the required information to run the diagnostic script
;
; No functions/procedures shall be further added to this script.
; #############################################################################

; Load config_user_info and diag_script_info from settings.ncl
begin
  print("INFO    Loading settings from " + getenv("settings"))
  loadscript("$settings")
end

; #############################################################################

; Load input_file_info, dataset_info and variable_info from <variable>_info.ncl
begin
  vardeffiles = diag_script_info@input_files
  nvdf = ListCount(vardeffiles)
  do i = 0, nvdf - 1
    if (.not. ismissing(vardeffiles[nvdf - 1 - i])) then
      vardeffile = vardeffiles[nvdf - 1 - i]
      print("INFO    Loading input data description from " + vardeffile)
      loadscript(vardeffile)
    end if
  end do
end

; #############################################################################

begin

  ; Add trailing slash to paths
  config_user_info@plot_dir = config_user_info@plot_dir + "/"
  config_user_info@run_dir = config_user_info@run_dir + "/"
  config_user_info@work_dir = config_user_info@work_dir + "/"

  ; Define current script name as global variables (following Python syntax)
  _DIAG_SCRIPT = diag_script_info@script
  _NVARIABLES = ListCount(variable_info)
  if (isvar("variable_info")) then
    _VARIABLES = new(_NVARIABLES, string)
    _FIELD_TYPES = new(_NVARIABLES, string)
    do i = 0, _NVARIABLES - 1
      _VARIABLES(i) = variable_info[i]@short_name
      _FIELD_TYPES(i) = variable_info[i]@field
    end do
  end if

end

; #############################################################################
; Load standard interface scripts
load "./interface_scripts/logging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"

; #############################################################################
; Convenience functions for selecting input_file_info

undef("select_metadata_by_var_name")
function select_metadata_by_var_name(info:list,
                                     variable:string,
                                     field:string)
;
; Arguments
;    info: the input_file_info logical containing all metadata
;    variable: variable's short name (according to CMOR standard)
;    field: variable's field type
;
; Return value
;    A list of logicals containing only the metadata matching the variable's
;    name and field.
;
; Description
;    Metadata information are stored in a list of input_file_info logicals.
;    This functions allows to extract the list elements matching the given
;    variable's short name and field type.
;
; Caveats
;
; References
;
; Modification history
;   20181105_A-righ_ma: written.
;
begin

  funcname = "select_metadata_by_var_name"
  scriptname = "interface_scripts/interface.ncl"
  enter_msg(scriptname, funcname)

  ninfo = ListCount(info)
  item = NewList("fifo")
  do ii = 0, ninfo - 1
    if (.not.isatt(info[ii], "short_name")) then
      error_msg("f", scriptname, funcname, "missing short_name attribute")
    end if
    if (.not.isatt(info[ii], "field")) then
      error_msg("f", scriptname, funcname, "missing field attribute")
    end if
    if (info[ii]@short_name.eq.variable .and. info[ii]@field.eq.field) then
      ListAppend(item, info[ii])
    end if
  end do

  leave_msg(scriptname, funcname)
  return(item)

end

; #############################################################################
undef("select_metadata_by_var_index")
function select_metadata_by_var_index(info:list,
                                      idx:integer)
;
; Arguments
;    info: the input_file_info logical containing all metadata
;    idx: the variable index as in the variable_info list (e.g, 0 for the
;         first variable, 1 for the second, etc.)
;
; Return value
;    A list of logicals containing only the metadata matching the variable's
;    index
;
; Description
;    Metadata information are stored in a list of input_file_info logicals.
;    This functions allows to extract the list elements matching the given
;    variable's index
;
; Caveats
;
; References
;
; Modification history
;   20181105_A-righ_ma: written.
;
begin

  funcname = "select_metadata_by_var_index"
  scriptname = "interface_scripts/interface.ncl"
  enter_msg(scriptname, funcname)

  if (idx.gt.ListCount(variable_info) - 1) then
    error_msg("f", funcname, scriptname, "variable's index out of range: " + \
              idx + " selected, only " + ListCount(variable_info) + \
              " available")
  end if

  var = variable_info[idx]
  item = select_metadata_by_var_name(info, var@short_name, var@field)

  leave_msg(scriptname, funcname)
  return(item)

end
