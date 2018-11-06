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

undef("select_metadata_by_atts")
function select_metadata_by_atts(info:list,
                                 atts:list)
;
; Arguments
;    info: a list of input_file_info logicals containing all metadata
;    atts: a list of arrays of two elements containing attribute=value pairs
;
; Return value
;    A list of logicals containing only the metadata matching the given
;    attribute value pairs
;
; Description
;    Metadata information are stored in a list of input_file_info logicals.
;    This functions allows to extract the list elements matching the attributes
;    value given as list in the second argument.
;
; Caveats
;
; References
;
; Modification history
;   20181106_A-righ_ma: written.
;
local funcname, scriptname, ninfo, natts, ii, jj, ca
begin

  funcname = "select_metadata_by_atts"
  scriptname = "interface_scripts/interface.ncl"
  enter_msg(scriptname, funcname)

  ninfo = ListCount(info)
  natts = ListCount(atts)
  item = NewList("fifo")
  do ii = 0, ninfo - 1  ; loop over input_file_info elements
    do jj = 0, natts - 1  ; loop over requested attribute=value pairs
      ca = atts[jj]
      if (.not.isatt(info[ii], ca[0])) then
        error_msg("f", scriptname, funcname, "missing " + ca[0] + " attribute")
      end if
      print(info[ii]@$ca[0]$ + " = " + ca[1])
      if (info[ii]@$ca[0]$.eq.ca[1]) then
        ListAppend(item, info[ii])
      end if
    end do
  end do

  leave_msg(scriptname, funcname)
  return(item)

end

; #############################################################################
undef("select_metadata_by_name")
function select_metadata_by_name(info:list,
                                 variable:string)

;
; Arguments
;    info: the input_file_info logical containing all metadata
;    variable: variable's short name (according to CMOR standard)
;
; Return value
;    A list of logicals containing only the metadata matching the variable's
;    short name
;
; Description
;    Metadata information are stored in a list of input_file_info logicals.
;    This functions allows to extract the list elements matching the given
;    variable's short name.
;
; Caveats
;
; References
;
; Modification history
;   20181105_A-righ_ma: written.
;
local funcname, scriptname
begin

  funcname = "select_metadata_by_name"
  scriptname = "interface_scripts/interface.ncl"
  enter_msg(scriptname, funcname)

  item = select_metadata_by_atts(info, [/[/"short_name", variable/]/])

  leave_msg(scriptname, funcname)
  return(item)

end

; #############################################################################
undef("select_metadata_by_index")
function select_metadata_by_index(info:list,
                                  idx:integer)
;
; Arguments
;    info: the input_file_info logical containing all metadata
;    idx: the variable index as in the variable_info list
;
; Return value
;    A list of logicals containing only the metadata matching the short name
;    of variable_info[idx]
;
; Description
;    Metadata information are stored in a list of input_file_info logicals.
;    This functions allows to extract the list elements matching the short name
;    of the variable with the given index
;
; Caveats
;
; References
;
; Modification history
;   20181105_A-righ_ma: written.
;
local funcname, scriptname
begin

  funcname = "select_metadata_by_index"
  scriptname = "interface_scripts/interface.ncl"
  enter_msg(scriptname, funcname)

  if (idx.gt.ListCount(variable_info) - 1) then
    error_msg("f", funcname, scriptname, "variable's index out of range")
  end if

  item = select_metadata_by_name(info, variable_info[idx]@short_name)

  leave_msg(scriptname, funcname)
  return(item)

end
