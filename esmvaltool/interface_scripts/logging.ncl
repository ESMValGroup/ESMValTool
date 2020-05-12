; #############################################################################
; GENERAL ROUTINES FOR HANDLING LOG OUTPUT
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    procedure log_info
;    procedure log_debug
;    procedure enter_msg
;    procedure leave_msg
;    procedure error_msg
;    procedure tool_stop
;    procedure exit_if_missing_atts
;    procedure log_provenance
;
; #############################################################################

; WARNING: no other interface_scripts shall be loaded here, to avoid nested
;          loops (since they will try to load logging.ncl themselves).

; #############################################################################
undef("log_info")
procedure log_info(output_string[*]:string)
;
; Arguments
;    output_string: the text to be output as message on screen
;
; Description
;    Write an info message to the log file. If the input is an array, each
;    element will be written on different lines.
;
; Caveats
;
; References
;
; Modification history
;    20180202-righi_mattia: written to replace info_output.
;
local nn, ii
begin

  if (any(config_user_info@log_level.eq.(/"info", "debug"/))) then
    nn = dimsizes(output_string)
    if (nn.eq.1) then
      print("INFO    " + output_string)
    else
      do ii = 0, nn - 1
        print("INFO    " + output_string(ii))
      end do
    end if
  end if

end

; #############################################################################
undef("log_debug")
procedure log_debug(output_string[*]:string)
;
; Arguments
;    output_string: the text to be output as message on screen
;
; Description
;    Write a debug message to the log file (only if log_level = debug in
;    config-user.yml). If the input is an array, each element will be
;    written on different lines.
;
; Caveats
;
; References
;
; Modification history
;    20180202-righi_mattia: written to replace info_output.
;
local nn, ii
begin

  if (config_user_info@log_level.eq."debug") then
    nn = dimsizes(output_string)
    if (nn.eq.1) then
      print("DEBUG   " + output_string)
    else
      do ii = 0, nn - 1
        print("DEBUG   " + output_string(ii))
      end do
    end if
  end if

end

; #############################################################################
undef("enter_msg")
procedure enter_msg(script[1]:string,
                    func[1]:string)
;
; Arguments
;    script: name of the ncl script file from where this function is called.
;    func: name of the function or procedure calling this function.
;
; Description
;    Print an enter message to the log file (debug mode only).
;
; Caveats
;
; References
;
; Modification history
;    20180202-righi_mattia: removed required verbosity.
;    20150319-righi_mattia: written.
;
begin

  if (func.eq."") then
    log_debug("<<< Entering " + script)
  else
    log_debug("<<< Entering " + func + " (" + script + ")")
  end if

end

; #############################################################################
undef("leave_msg")
procedure leave_msg(script[1]:string,
                    func[1]:string)
;
; Arguments
;    script: name of the ncl script file from where this function is called.
;    func: name of the function or procedure calling this function.
;
; Description
;    Print a leave message to the log file (debug mode only).
;
; Caveats
;
; References
;
; Modification history
;    20180202-righi_mattia: removed required verbosity.
;    20150319-righi_mattia: written.
;
begin

  if (func.eq."") then
    log_debug(">>> Leaving " + script)
  else
    log_debug(">>> Leaving " + func + " (" + script + ")")
  end if

end

; #############################################################################
undef("error_msg")
procedure error_msg(type[1]:string,
                    script[1]:string,
                    func[1]:string,
                    msg:string)
;
; Arguments
;    type: type of message, "f" for fatal, "w" for warning.
;    script: name of the ncl script file from where this function is called.
;    func: name of the function or procedure calling this function.
;    msg: actual error message.
;
; Description
;    Write an error/warning message to the log file.
;
; Caveats
;
; References
;
; Modification history
;    20150924-lauer_axel: added type "error" (treated as "fatal").
;    20140929-righi_mattia: written.
;
local nn, ii, msg_level, msg_origin
begin

  nn = dimsizes(msg)

  ; Determine which message level we're at (fatal/warning)
  if (type .eq. "f" .or. type .eq. "fatal") then
    msg_level = "fatal"
  else if (any(type.eq.(/"w", "warning", "error"/))) then
    msg_level = "warning"
  else
    print("fatal: in error_msg (logging.ncl), " + \
          "invalid value for variable type (=" + type + ")")
    status_exit(1)
  end if
  end if

  ; Determine origin of message
  if (func .ne. "") then
    msg_origin = func + " (" + script + "), "
  else
    msg_origin = script + ", "
  end if

  ; Write error/warning message
  do ii = 0, nn - 1
    log_info(msg_level + ": in " + msg_origin + msg(ii))
  end do

  ; Exit if fatal
  if (msg_level .eq. "fatal") then
    status_exit(1)
  end if

end

; #############################################################################
undef("tool_stop")
procedure tool_stop()
;
; Arguments
;
; Description
;    Force the tool execution to stop at any point of an ncl script.
;    For dubugging purposes only.
;
; Caveats
;    This should not be used for the standard exit-on-error. See the function
;    error_msg in logging.ncl instead.
;
; References
;
; Modification history
;    20150521-righi_mattia: written.
;
begin

  print("fatal: tool stop forced by the user")
  status_exit(1)

end

; #############################################################################
undef("exit_if_missing_atts")
procedure exit_if_missing_atts(in[1], \
                               reqatts[*]:string)
;
; Arguments
;    in: a logical or file variable to be checked for given attributes
;    reqatts: an array of strings with the required attributes of the variable
;
; Description
;    Check that the given variable contains the given list of attributes.
;    This function shall be used at the beginning of each diag_script to
;    check that diag_script_info contains the required attributes.
;    It can also be used for checking attributes in observations files, when
;    these are read directly in the diag_script (e.g. Emmons.ncl).
;
; Caveats
;
; References
;
; Modification history
;    20141002-righi_mattia: written.
;
local funcname, scriptname
begin

  funcname = "exit_if_missing_atts"
  scriptname = "interface_scripts/logging.ncl"
  enter_msg(scriptname, funcname)

  if (.not.all(isatt(in, reqatts))) then
    error_msg("f", DIAG_SCRIPT, "", "missing required diag_script_info " + \
              "attribute: " + reqatts(ind(.not.isatt(in, reqatts))))
  end if

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("log_provenance")
procedure log_provenance(nc_file:string,
                         outfile:string,
                         caption:string,
                         statistics[*]:string,
                         domains[*]:string,
                         plot_types[*]:string,
                         diag_authors[*]:string,
                         diag_refs[*]:string,
                         infiles[*]:string)
;
; Arguments
;    nc_file: netcdf file with data related to the plot
;    outfile: file name of the figure (including path)
;    caption: figure caption
;    statistics: list of metrics/statistics calculated
;    domains: of the data displayed
;    plot_types: plot type of the figure
;    diag_authors: authors related to the diagnostic
;    diag_refs: references related to the diagnotic
;    infiles: input files (from preproc) used to generate the plot
;
; Description
;    Create a yaml file with meta data from the diagnostic script and save it
;    in the run directory
;
; Caveats
;
; References
;
; Modification history
;    20190415-righi_mattia: turn domains into a list.
;    20190415-righi_mattia: extended to avoid overwriting previous call.
;    20190225-bock_lisa: written.
;
local funcname, scriptname, yaml_file, outstring, existing
begin

  funcname = "log_provenance"
  scriptname = "interface_scripts/logging.ncl"
  enter_msg(scriptname, funcname)

  ; Define output
  outstring = new(9, string)

  ; Set yaml file path
  yaml_file = config_user_info@run_dir + "diagnostic_provenance.yml"

  ; Save entries
  outstring(0) = "? " + nc_file
  outstring(1) = ": ancestors: [" + str_join(infiles, ", ") + "]"
  outstring(2) = "  authors: [" + str_join(diag_authors, ", ") + "]"
  outstring(3) = "  caption: '" + caption + "'"
  outstring(4) = "  domains: [" + str_join(domains, ", ") + "]"
  outstring(5) = "  plot_file: " + outfile
  outstring(6) = "  plot_types: [" + str_join(plot_types, ", ") + "]"
  outstring(7) = "  references: [" + str_join(diag_refs, ", ") + "]"
  outstring(8) = "  statistics: [" + str_join(statistics, ", ") + "]"

  ; Save existing information to avoid overwriting
  if (fileexists(yaml_file)) then
    existing = asciiread(yaml_file, -1, "string")
    outstring := array_append_record(existing, outstring, 0)
    delete(existing)
  end if

  ; Save provenance information
  asciiwrite(yaml_file, outstring)

  log_info("Provenance information saved in " + yaml_file)

  leave_msg(scriptname, funcname)

end
