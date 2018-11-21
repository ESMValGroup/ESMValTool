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
  nvdf = dimsizes(vardeffiles)
  do i = 0, nvdf - 1
    isncl = .not.ismissing(str_match_ic_regex(vardeffiles(i), \
                                              ".*_info\" + ".ncl"))
    if (.not. ismissing(vardeffiles(i)) .and. isncl) then
      vardeffile = vardeffiles(i)
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

  ; Define current script name as global variable (following Python syntax)
  DIAG_SCRIPT = diag_script_info@script

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
                                 att[1]:logical)
;
; Arguments
;    info: a list of input_file_info logicals containing all metadata
;    atts: a scalar logical with the attributes to be matched in the input list
;
; Return value
;    A list of logicals containing only the metadata matching the given
;    attribute value pairs
;
; Description
;    Metadata information are stored in a list of input_file_info logicals.
;    This functions allows to extract the list elements whose attributes
;    match the requested ones.
;
; Caveats
;
; References
;
; Modification history
;   20181106_A-righ_ma: written.
;
local funcname, scriptname, ninfo, natts, l_flag, ii, jj, ca
begin

  funcname = "select_metadata_by_atts"
  scriptname = "interface_scripts/interface.ncl"
  enter_msg(scriptname, funcname)

  ninfo = ListCount(info)
  atts = getvaratts(att)
  natts = dimsizes(atts)
  item = NewList("fifo")
  do ii = 0, ninfo - 1  ; loop over input_file_info elements
    l_flag = True
    do jj = 0, natts - 1  ; loop over requested attributes
      if (.not.isatt(info[ii], atts(jj))) then
        error_msg("f", scriptname, funcname, \
                  "missing " + atts(jj) + " attribute")
      end if
      if (info[ii]@$atts(jj)$.ne.att@$atts(jj)$) then
        l_flag = False
      end if
    end do
    if (l_flag) then
      ListAppend(item, info[ii])
    end if
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

  atts = True
  atts@short_name = variable
  item = select_metadata_by_atts(info, atts)

  leave_msg(scriptname, funcname)
  return(item)

end
