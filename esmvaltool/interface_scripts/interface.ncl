; #############################################################################
; INTERFACE TO HANDLE THE COMMUNICATION BETWEEN THE PYTHON WORKFLOW AND NCL
; #############################################################################
; Load the automatically-generated interface files <variable>_info.ncl and
; settings.ncl for the current diagnostic and variable with the current
; preprocessor set. These files contains the NCL logical _info varialbe with
; all the required information to run the diagnostic script
;
; Contents:
;    function select_metadata_by_atts
;    function select_metadata_by_name
;    function metadata_att_as_array
;    procedure ESMValMD
;
; No further functions/procedures shall be added to this script.
;
; #############################################################################

; Load config_user_info and diag_script_info from settings.ncl
begin
  if .not. isvar("__interface_ncl__") then
    print("INFO    Loading settings from " + getenv("settings"))
    loadscript("$settings")
  end if
end

; #############################################################################

; Load input_file_info, dataset_info and variable_info from <variable>_info.ncl
begin
  if .not. isvar("__interface_ncl__") then
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
  end if
end

; #############################################################################

begin
  if .not. isvar("__interface_ncl__") then
    ; Add trailing slash to paths
    config_user_info@plot_dir = config_user_info@plot_dir + "/"
    config_user_info@run_dir = config_user_info@run_dir + "/"
    config_user_info@work_dir = config_user_info@work_dir + "/"

    ; Define current script name as global variable (following Python syntax)
    DIAG_SCRIPT = diag_script_info@script
  end if
end

; #############################################################################
; Load standard interface scripts
load "./interface_scripts/logging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"

; #############################################################################
; Convenience functions for selecting input_file_info and related attributes

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

; #############################################################################
undef("metadata_att_as_array")
function metadata_att_as_array(inlist:list,
                               att[1]:string)
;
; Arguments
;    inlist: a list
;    att: a string with the attribute to be searched for
;
; Description
;    This function retrieves the values of the attribute att from each element
;    of the list and returns them as an array.
;
; Caveats
;
; References
;
; Modification history
;    20181119-A_laue_ax: missing attributes are now allowed; error message is
;                        only issued if all attributes are missing
;    20180711-A_righ_ma: written
;
local funcname, scriptname
begin

  funcname = "metadata_att_as_array"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  do ii = 0, ListCount(inlist) - 1
    if (isatt(inlist[ii], att)) then
      if (.not.isdefined("out")) then
        set_type = typeof(inlist[ii]@$att$)
        out = new(ListCount(inlist), set_type)
      end if
      if (typeof(inlist[ii]@$att$).ne.set_type) then
        error_msg("f", scriptname, funcname, "attributes must be of the " + \
                  "same type")
      end if
      out(ii) = inlist[ii]@$att$
    end if
  end do

  if (all(ismissing(out))) then
    error_msg("f", scriptname, funcname, "attribute " + att + " not " + \
              "found in list element(s)")
  end if

  return(out)
  leave_msg(scriptname, funcname)

end

; Remember that the interface is now loaded.
begin
  __interface_ncl__ = True
end


; #############################################################################
undef("ESMValMD")
procedure ESMValMD(nc_file:string,
                   outfile:string,
                   caption:string,
                   statistics[*]:string,
                   domain:string,
                   plottype:string,
                   diag_authors[*]:string,
                   diag_refs[*]:string,
                   infiles[*]:string)
;
; Arguments
;  nc_file: netcdf file with data related to the plot
;  outfile: file name of the figure (including path)
;  caption: figure caption
;  statistics: list of metrics/statistics calculated
;  domain: of the data displayed
;  plottype: plot type of the figure
;  diag_authors: authors related to the diagnostic
;  diag_refs: references related to the diagnotic
;
; Description
;   Creates a yaml file with meta data from the diagnostic
;   script and save it in the plot directory
;
; Caveats
;
; References
;
local funcname, scriptname, verbosity
begin

  funcname = "ESMValMD"
  scriptname = "diag_scripts/shared/lib/meta_data.ncl"
  enter_msg(scriptname, funcname)

  str0 = "? " + nc_file
  n = dimsizes(infiles)
  str1 = ": ancestors: ["
  str1 = str1 + infiles(0)
  if(n.ge.1) then
    do i = 1, n-1
      str1 = str1 + ", " + infiles(i)
    end do
  end if
  str1 = str1 + "]"
  n = dimsizes(diag_authors)
  str2 = "  authors: ["
  str2 = str2 + diag_authors(0)
  if(n.ge.1) then
    do i = 1, n-1
      str2 = str2 + ", " + diag_authors(i)
    end do
  end if
  str2 = str2 + "]"
  str3 = "  caption: '" + caption + "'"
  str4 = "  domain: " + domain
  str5 = "  plotfile: " + outfile
  str6 = "  plot_type: " + plottype
  n = dimsizes(diag_refs)
  str7 = "  references: ["
  str7 = str7 + diag_refs(0)
  if(n.ge.1) then
    do i = 1, n-1
      str7 = str7 + ", " + diag_refs(i)
    end do
  end if
  str7 = str7 + "]"
  n = dimsizes(statistics)
  str8 = "  statistics: ["
  str8 = str8 + statistics(0)
  if(n.ge.1) then
    do i = 1, n-1
      str8 = str8 + ", " + statistics(i)
    end do
  end if
  str8 = str8 + "]"

  str := array_append_record(str7, str8, 0)
  str := array_append_record(str6, str, 0)
  str := array_append_record(str5, str, 0)
  str := array_append_record(str4, str, 0)
  str := array_append_record(str3, str, 0)
  str := array_append_record(str2, str, 0)
  str := array_append_record(str1, str, 0)
  str := array_append_record(str0, str, 0)

  str_list = [/str/]

  yaml_file = config_user_info@run_dir + "diagnostic_provenance.yml"
  write_table(yaml_file, "a", str_list, "%s %s %s %s %s %s %s %s %s")
  delete(str)

  log_info(" write meta data to " + yaml_file)

  leave_msg(scriptname, funcname)

end
