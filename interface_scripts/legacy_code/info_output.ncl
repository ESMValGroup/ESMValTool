;;#############################################################################
;; GENERAL ROUTINES FOR HANDLING LOG OUTPUT
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure no_operation
;;    procedure info_output
;;    procedure error_msg
;;    procedure check_atts
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;#############################################################################

undef("no_operation")
procedure no_operation()
;; Description
;;   Defines an NCL no-op operator
begin
end

undef("info_output")
procedure info_output(output_string[1]:string,
                      verbosity[1]:integer,
                      required_verbosity[1]:integer)
;;
;; Argument
;;   output_string: the text to be output as message on screen
;;   verbosity: the current verbosity set as global attribute in the main nml
;;   required_verbosity: the minimum verbosity for the message to be printed
;;
;; Description
;;   Writes a log message on the stdout depending on the selected verbosity.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20130119-A_eval_ma: written.
;;
local indent, indent_lev, indent_tmp, indent_str, indent_step, cwd
begin

    cwd  = getenv("ESMValTool_cwd")
    if (ismissing(cwd)) then
        if(.not.ismissing(verbosity)) then
            print("info: " + output_string)
        end if 
    else  ; Full info output
        indent_step = 3
        if (verbosity .eq. 2) then
            padding_space = 13
        else
            padding_space = 1
        end if

        if (verbosity.ge.required_verbosity) then
            ;; Increase indentation if we're entering an NCL routine
            if (isStrSubset(output_string, "<<<<<<<< Entering")) then
                indent = asciiread(cwd + "/interface_data/curr_trace_indent.txt",\
                                         1, "integer")
                indent = indent + indent_step
                asciiwrite(cwd + "/interface_data/curr_trace_indent.txt",\
                           indent)
            end if

            ;; Add a prefix of ">>", "<<" or "  " to output_string
            indent = asciiread(cwd + "/interface_data/curr_trace_indent.txt",\
                                     1, "integer")
            if (indent .lt. 1) then
                indent = 1
            end if
            indent_lev = new(indent + padding_space, "string")
            indent_tmp = where(ismissing(indent_lev), " ", "-")
            if (isStrSubset(output_string, "<")) then
                delete(indent_lev)
                delete(indent_tmp)
                indent_lev = new(indent, "string")
                indent_tmp = where(ismissing(indent_lev), "<", "-")
            end if
            if (isStrSubset(output_string, ">")) then
                delete(indent_lev)
                indent_lev = new(indent, "string")
                delete(indent_tmp)
                indent_tmp = where(ismissing(indent_lev), ">", "-")
            end if

            indent_str = str_concat(indent_tmp)
            print("info: " + indent_str + output_string)

            ;; Decrease indentation if we're leaving an NCL routine
            if (isStrSubset(output_string, ">>>>>>>> Leaving")) then
                indent = asciiread(cwd + "/interface_data/curr_trace_indent.txt",\
                                         1, "integer")
                indent = indent - indent_step
                asciiwrite(cwd + "/interface_data/curr_trace_indent.txt",\
                           indent)
            end if
        end if
    end if
end

;;#############################################################################
undef("error_msg")
procedure error_msg(type[1]:string, 
                    script[1]:string,
                    func[1]:string,
                    msg:string)
;;
;; Arguments
;;   type: type of message, "f" for fatal, "w" for warning.
;;   script: name of the ncl script file from where this function is called.
;;   func: name of the function or procedure calling this function.
;;   msg: actual error message.
;;   
;; Description
;;   Writes a log message on the stdout depending on the selected verbosity.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20140929-A_righ_ma: written.
;;
local nn, ii, msg_level, msg_origin
begin
    nn = dimsizes(msg)

    ;; Determine which message level we're at (fatal/warning) 
    if (type .eq. "f" .or. type .eq. "fatal") then
        msg_level = "fatal"
    else if (type .eq. "w" .or. type .eq. "warning") then
             msg_level = "warning"
    else
        print("fatal: in error_msg (info_output.ncl), invalid value for variable type (=" + type + ")")
        status_exit(1)
    end if
    end if

    ;; Determine origin of message
    if (func .ne. "") then
        msg_origin = func + " (" + script + "), "
    else
        msg_origin = script + ", "
    end if

    ;; Write error/warning message
    do ii = 0, nn - 1
        info_output(msg_level + ": in " + msg_origin + msg(ii), 1, 1)
    end do

    ;; Exit if fatal
    if (msg_level .eq. "fatal") then
        status_exit(1)
    else
        no_operation
    end if
end

;;#############################################################################
undef("check_atts")
procedure check_atts(in, \
                     reqatts[*]:string)
;;
;; Arguments
;;   in: a logical or file variable to be checked for given attributes
;;   reqatts: an array of strings with the required attributes of the variable
;;   
;; Description
;;   Checks that the given variable contains the given list of attributes.
;;   This function shall be used at the beginning of each diag_script to
;;   check that diag_script_info contains the required attributes.
;;   It can also be used for checking attributes in observations files, when
;;   these are read directly in the diag_script (e.g. Emmons.ncl).
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20141002-A_righ_ma: written.
;;
begin

    if (.not.all(isatt(in, reqatts))) then
        error_msg("f", diag_script, "", "missing required " + \
                  "diag_script_info attribute: " + \
                  reqatts(ind(.not.isatt(in, reqatts))))
    end if

end

undef("exit_if_missing_attributes")
procedure exit_if_missing_attributes(in[1]:logical,\
                                     reqatts[*]:string)
;; Description
;;   Wrapper for function "check_atts(...)" which terminates execution if any
;;   of the required_attribues 'reqatts' are missing from the logical 'in'
begin
    check_atts(in, reqatts)
end
