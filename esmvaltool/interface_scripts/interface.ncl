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

; Convert diag_script_info and config_user_info to scalar
begin
  tmp = diag_script_info[0]
  delete(diag_script_info)
  diag_script_info = tmp
  delete(tmp)
  delete(diag_script_info@_FillValue)

  tmp = config_user_info[0]
  delete(config_user_info)
  config_user_info = tmp
  delete(tmp)
  delete(config_user_info@_FillValue)
end


; Load input_file_info, dataset_info and variable_info from <variable>_info.ncl
begin
  vardeffiles = \
    str_match_ic_regex(diag_script_info@input_files, ".*_info\" + ".ncl")
  nvdf = dimsizes(vardeffiles)
  do i = 0, nvdf - 1
    if (.not. ismissing(vardeffiles(nvdf - 1 - i))) then
      vardeffile = diag_script_info@input_files(nvdf - 1 - i)
      print("INFO    Loading input data description from " + vardeffile)
      loadscript(vardeffile)
    end if
  end do
end

begin

  ; Add trailing slash to paths
  config_user_info@plot_dir = config_user_info@plot_dir + "/"
  config_user_info@run_dir = config_user_info@run_dir + "/"
  config_user_info@work_dir = config_user_info@work_dir + "/"

  ; Copy some info into legacy names (for backward compatibility with v1.0)
  diag_script = diag_script_info@script
  if (isvar("variable_info")) then
    variables = new(ListCount(variable_info), string)
    field_types = new(ListCount(variable_info), string)
    do i = 0, ListCount(variable_info) - 1
      variables(i) = variable_info[i]@short_name
      field_types(i) = variable_info[i]@field
    end do
  end if

  ; Convert lists to scalars if they contain only 1 element
  if (isvar("variable_info")) then
    if (ListCount(variable_info).eq.1) then
      tmp = variable_info[0]
      delete(variable_info)
      variable_info = tmp
      delete(tmp)
    end if
  end if

  if (isvar("dataset_info")) then
    if (ListCount(dataset_info).eq.1) then
      tmp = dataset_info[0]
      delete(dataset_info)
      dataset_info = tmp
      delete(tmp)
    end if
  end if

  if (isvar("input_file_info")) then
    if (ListCount(input_file_info).eq.1) then
      tmp = input_file_info[0]
      delete(input_file_info)
      input_file_info = tmp
      delete(tmp)
    end if
  end if

end

; Load other standard interface scripts
load "./interface_scripts/logging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
