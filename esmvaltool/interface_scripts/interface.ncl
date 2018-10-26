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

; Convert list to scalar (config_user_info and diag_script_info contain
; always only 1 element)
begin
  tmp = config_user_info[0]
  delete(config_user_info)
  config_user_info = tmp
  delete(tmp)
  tmp = diag_script_info[0]
  delete(diag_script_info)
  diag_script_info = tmp
  delete(tmp)
end

; Load input_file_info, dataset_info and variable_info from <variable>_info.ncl
begin
  vardeffiles = \
    str_match_ic_regex(diag_script_info@input_files, ".*_info\" + ".ncl")

  print(vardeffiles)
  do i = 0, dimsizes(vardeffiles) - 1
    if (.not. ismissing(vardeffiles(i))) then
      vardeffile = diag_script_info@input_files(i)
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

  ; Copy some info into legacy names
  diag_script = diag_script_info@script
  if (isvar("variable_info")) then
    variables = new(ListCount(variable_info), string)
    field_types = new(ListCount(variable_info), string)
    do i = 0, ListCount(variable_info) - 1
      variables(i) = variable_info[i]@short_name
      field_types(i) = variable_info[i]@field
    end do
  end if

end

; Load other interface scripts
load "./interface_scripts/logging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
