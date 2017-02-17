;;#############################################################################
;; GENERAL ROUTINES FOR HANDLING LOG OUTPUT
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure no_operation
;;    procedure info_output
;;    procedure inquire_and_save_fileinfo
;;    procedure enter_msg
;;    procedure leave_msg
;;    procedure error_msg
;;    procedure exit_if_missing_atts
;;    procedure write_references
;;    procedure write_filelist
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;; Warning: no other interface_scripts shall be loaded here, to avoid nested
;;          loops (since they will try to load messaging.ncl themselves).

; a_laue_ax+
    gOldVar = ""
    gOldDiag = ""

    ;; definitions of functions write_header and write_diag_header
    load "./interface_scripts/write_header.ncl"
; a_laue_ax-

;;#############################################################################
undef("no_operation")
procedure no_operation()
;;
;; Arguments
;;
;; Description
;;    Define an NCL no-op operator.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;     20130119-A_eval_ma: written.
;;
begin

end

;;#############################################################################
undef("info_output")
procedure info_output(output_string[*]:string,
                      verbosity[1]:integer,
                      required_verbosity[1]:integer)
;;
;; Arguments
;;    output_string: the text to be output as message on screen
;;    verbosity: the current verbosity set as global attribute in the main nml
;;    required_verbosity: the minimum verbosity for the message to be printed
;;
;; Description
;;    Writes a log message on the stdout depending on the selected verbosity.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130119-A_eval_ma: written.
;;
local indent, indent_lev, indent_tmp, indent_str, indent_step, cwd
begin

    required_verbosity = min((/required_verbosity, 4/))

    cwd  = getenv("ESMValTool_cwd")
    if (ismissing(cwd)) then
        print("info: " + output_string)
    else  ; Full info output
        indent_step = 3
        if (verbosity .eq. 2) then
            padding_space = 13
        else
            padding_space = 1
        end if

        if (verbosity .ge. required_verbosity) then
            ;; Increase indentation if we're entering an NCL routine
            if (dimsizes(output_string) .eq. 1) then
                if (isStrSubset(output_string, "<<<<<<<< Entering")) then
                    indent = asciiread(cwd + \
                                       "/interface_data/curr_trace_indent.txt",\
                                       1, "integer")
                    indent = indent + indent_step
                    asciiwrite(cwd + "/interface_data/curr_trace_indent.txt",\
                               indent)
                end if
            end if

            ;; Add a prefix of ">>", "<<" or "  " to output_string
            indent = asciiread(cwd + "/interface_data/curr_trace_indent.txt",\
                                     1, "integer")
            if (indent .lt. 1) then
                indent = 1
            end if
            indent_lev = new(indent + padding_space, "string")
            indent_tmp = where(ismissing(indent_lev), " ", "-")
            if (dimsizes(output_string) .eq. 1) then
                if (isStrSubset(output_string, "<")) then
                    delete(indent_lev)
                    delete(indent_tmp)
                    indent_lev = new(indent, "string")
                    indent_tmp = where(ismissing(indent_lev), "<", "-")
                end if
            end if
            if (dimsizes(output_string) .eq. 1) then
                if (isStrSubset(output_string, ">")) then
                    delete(indent_lev)
                    indent_lev = new(indent, "string")
                    delete(indent_tmp)
                    indent_tmp = where(ismissing(indent_lev), ">", "-")
                end if
            end if

            indent_str = str_concat(indent_tmp)
            print("info: " + indent_str + output_string)

            ;; Decrease indentation if we're leaving an NCL routine
            if (dimsizes(output_string) .eq. 1) then
                if (isStrSubset(output_string, ">>>>>>>> Leaving")) then
                    indent = asciiread(cwd + \
                                       "/interface_data/curr_trace_indent.txt", \
                                       1, "integer")
                    indent = indent - indent_step
                    asciiwrite(cwd + "/interface_data/curr_trace_indent.txt", \
                               indent)
                end if
            end if
        end if
    end if
end

;;#############################################################################
undef("enter_msg")
procedure enter_msg(script[1]:string,
                    func[1]:string,
                    req_verbosity[1]:integer)
;;
;; Arguments
;;    script: name of the ncl script file from where this function is called.
;;    func: name of the function or procedure calling this function.
;;    required_verbosity: the minimum verbosity for the message to be printed.
;;
;; Description
;;    Prints an enter message on the stdout.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150319-A_righ_ma: written.
;;
local verbosity
begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    if (func.eq."") then
        info_output("<<<<<<<< Entering " + script, verbosity, req_verbosity)
    else
        info_output("<<<<<<<< Entering " + func + " (" + script + ")", \
                    verbosity, req_verbosity)
    end if

end

;;#############################################################################
undef("leave_msg")
procedure leave_msg(script[1]:string,
                    func[1]:string,
                    req_verbosity[1]:integer)
;;
;; Arguments
;;    script: name of the ncl script file from where this function is called.
;;    func: name of the function or procedure calling this function.
;;    required_verbosity: the minimum verbosity for the message to be printed.
;;
;; Description
;;    Prints a leave message on the stdout.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150319-A_righ_ma: written.
;;
local verbosity
begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    if (func.eq."") then
        info_output(">>>>>>>> Leaving " + script, verbosity, req_verbosity)
    else
        info_output(">>>>>>>> Leaving " + func + " (" + script + ")", \
                    verbosity, req_verbosity)
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
;;    type: type of message, "f" for fatal, "w" for warning.
;;    script: name of the ncl script file from where this function is called.
;;    func: name of the function or procedure calling this function.
;;    msg: actual error message.
;;
;; Description
;;    Writes a log message on the stdout depending on the selected verbosity.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150924-A_laue_ax: added type "error" (treated as "fatal")
;;    20140929-A_righ_ma: written.
;;
local nn, ii, msg_level, msg_origin
begin

    nn = dimsizes(msg)

    ;; Determine which message level we're at (fatal/warning)
    if (type .eq. "f" .or. type .eq. "fatal") then
        msg_level = "fatal"
; A-laue_ax+
;    else if (type .eq. "w" .or. type .eq. "warning") then
    else if (type .eq. "w" .or. type .eq. "warning" .or. type .eq. "error") then
; A-laue_ax-
        msg_level = "warning"
    else
        print("fatal: in error_msg (messaging.ncl), " + \
              "invalid value for variable type (=" + type + ")")
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
        no_operation()
    end if
end

;;#############################################################################
undef("exit_if_missing_atts")
procedure exit_if_missing_atts(in[1], \
                               reqatts[*]:string)
;;
;; Arguments
;;    in: a logical or file variable to be checked for given attributes
;;    reqatts: an array of strings with the required attributes of the variable
;;
;; Description
;;    Checks that the given variable contains the given list of attributes.
;;    This function shall be used at the beginning of each diag_script to
;;    check that diag_script_info contains the required attributes.
;;    It can also be used for checking attributes in observations files, when
;;    these are read directly in the diag_script (e.g. Emmons.ncl).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141002-A_righ_ma: written.
;;
local funcname, scriptname, verbosity
begin

    funcname = "exit_if_missing_atts"
    scriptname = "interface_scripts/messaging.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    if (.not.all(isatt(in, reqatts))) then
        error_msg("f", diag_script, "", "missing required " + \
                  "diag_script_info attribute: " + \
                  reqatts(ind(.not.isatt(in, reqatts))))
    end if

    leave_msg(scriptname, funcname, 4)

end

; a_laue_ax+

;;#############################################################################
undef("write_filelist")
procedure write_filelist(script:string,
                         variable:string,
                         model:string,
                         filelist:string,
                         sourcefilelist:string,
                         fixfile:string,
                         sourcefixfile:string,
                         sourcereflist:string,
                         sourcetidlist:string,
                         version:string)
;;
;; Arguments
;;    script: name of the diagnostic script for which the input files are
;;            written to the log file
;;    variable: name of variable
;;    model: name of model
;;    filelist: list of filenames processed to produce the file to be written
;;    sourcefilelist: list of filenames processed to produce the input files
;;                    of the file to be written (e.g., input files used to
;;                    produce climo files then used to process write a derived
;;                    variable)
;;    fixfile: filename of a FIXFILE applied to process the data to be written
;;             (if applicable)
;;    sourcefixfile: filename of a FIXFILE applied to the files in
;;                   sourcefilelist (if applicable)
;;    sourcereflist: references from original source files (if available)
;;    sourcetidlist: tracking_id(s) of original source files (if available)
;;    version: version of the ESMValTool
;;
;; Description
;;    Saving the following details to the log file
;;        - ESMValTtool version
;;        - list of input files processed
;;        - name of the "FIXFILE" applied (if any)
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150630-A_laue_ax: written.
;;
local funcname, scriptname, verbosity, output_refs, hline
begin

    funcname = "write_filelist"
    scriptname = "interface_scripts/messaging.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    if (filelist(0).eq."") then  ; derived variable or no information
        return
    end if

    hline = "-------------------------"
    hline = hline + hline + hline

    ;; Output refs file
    output_refs = getenv("ESMValTool_out_refs")

    s_open  = "echo " + str_get_dq
    s_close = str_get_dq + " >> " + output_refs

;    ;; If first time (empty output_refs) write header
    if (toint(systemfunc("wc -c < " + output_refs)) .eq. 0) then
        write_header(s_open, hline, s_close)
    end if

    if (script.ne.gOldDiag) then
        ;; Add script name
;        system("echo " + " >> " + output_refs)
;        system(s_open + "=== " + script + " ===" + s_close)
;        system("echo " + " >> " + output_refs)
        system(s_open + "PREPROCESSING/REFORMATTING (ESMValTool v" + \
               version + "):" + s_close)
        system("echo " + " >> " + output_refs)
        gOldDiag = script
    end if

    nfields = str_fields_count(filelist(0), "/")
    if (nfields.gt.1) then
        path = "/"
        do i = 1, nfields - 1
            path = path + str_get_field(filelist(0), i, "/") + "/"
        end do
    else
        path = ""
    end if

    if (gOldVar.ne.variable) then
        system(s_open + "  Variable: " + variable + s_close)
        system("echo " + " >> " + output_refs)
    end if
    gOldVar = variable
    system(s_open + "    Model: " + model + s_close)
    system(s_open + "    Input path: " + path + s_close)
    system(s_open + "    Input file(s):" + s_close)

    dsizes = dimsizes(filelist)
    do i = 0, dsizes(0) - 1
        nfields = str_fields_count(filelist(i), "/")
        fname = str_get_field(filelist(i), nfields, "/")
        system(s_open + "      " + sprinti("(%i) ",i + 1) + fname + s_close)
    end do

    dsizes = dimsizes(sourcefilelist)
    if (sourcefilelist(0).ne."") then
        system(s_open + "      Original source file(s) of all input file(s):"  \
               + s_close)
        do i = 0, dsizes(0) - 1
            if (.not.ismissing(sourcetidlist(i))) then
                system(s_open + "        -S- " + sprinti("(%i) ",i + 1)  \
                       + sourcefilelist(i) + " (tracking_id: " \
                       + sourcetidlist(i) + ") " + s_close)
            else
                system(s_open + "        -S- " + sprinti("(%i) ",i + 1)  \
                       + sourcefilelist(i) + s_close)
            end if
;            if (.not.ismissing(sourcetidlist(i))) then
;                system(s_open + "            --> tracking_id = " \
;                   + sourcetidlist(i) + s_close)
;            end if
        end do
    end if

    if ((sourcefixfile.ne."").and.(sourcefilelist(0).ne."")) then
        system(s_open + "      Fixes applied to original source file(s): "  \
               + sourcefixfile + s_close)
    else
        system(s_open + "    Fixes applied to input file(s): " + fixfile  \
               + s_close)
    end if

    dsizes = dimsizes(sourcereflist)
    if (sourcereflist(0).ne."") then
        system(s_open + "    Reference(s) of original source file(s):"  \
               + s_close)
        do i = 0, dsizes(0) - 1
            system(s_open + "      " + sprinti("(%i) ",i + 1)  \
                   + sourcereflist(i) + s_close)
        end do
    end if

;    system(s_open + hline + s_close)
    system("echo " + " >> " + output_refs)

    leave_msg(scriptname, funcname, 4)
end

; a_laue_ax-

;;#############################################################################
undef("write_references")
procedure write_references(script[1]:string,
                           auth:string,
                           contr:string,
                           diag:string,
                           obs:string,
                           proj:string)
;;
;; Arguments
;;    script: name of the diagnostic script for which reference are written
;;    auth: list of authors tags (A_YyXx)
;;    contr: list of cotributors tags (A_YyXx)
;;    diag: list of references for the diag_scripte diagnostics (DNNN)
;;    obs: list of references for the observations (ONNN)
;;    proj: list of projects to be acknowledged (PNNN)
;;
;; Description
;;    Writes references and acknowledgments in the log file
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151028-A_laue_ax: moved writing of diagnostic section header to
;;                        separate routine (write_diag_header) for easier
;;                        use with Python routines
;;    20150630-A_laue_ax: moved writing of header to separate routine
;;                        (write_header)
;;    20121128-A_righ_ma: written.
;;
; a_laue_ax+
;local funcname, scriptname, verbosity, namelist, master_refs, output_refs, \
local funcname, scriptname, verbosity, master_refs, output_refs, \
; a_laue_ax-
    hline, pos, str ;, stand_alone_entries
begin

    funcname = "write_references"
    scriptname = "interface_scripts/messaging.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    hline      = "-------------------------"
    hline_dble = "========================="
    hline = hline + hline + hline
    hline_dble = hline_dble + hline_dble + hline_dble

;    ;; D*** entries that shall not be extended by "Please cite "
;    stand_alone_entries = (/"D000", "D003", "D004"/)

; a_laue_ax+
;    ;; Namelist in use
;    namelist = getenv("ESMValTool_xml_name")
; a_laue_ax-

    ;; Master refs file
    master_refs = getenv("ESMValTool_in_refs")

    ;; Output refs file
    output_refs = getenv("ESMValTool_out_refs")

    s_open  = "echo " + str_get_dq
    s_close = str_get_dq + " >> " + output_refs

    ;; If first time (empty output_refs) write header
    if (toint(systemfunc("wc -c < " + output_refs)) .eq. 0) then
; a_laue_ax+
;        system(s_open + hline + s_close)
;        system(s_open + "+++++++++++++ " + \
;               "ESMVal TOOL REFERENCES and ACKNOWLEDGEMENTS LOG" + \
;               " +++++++++++++" + s_close)
;        system(s_open + hline + s_close)
;        system(s_open + "" + s_close)
;        system(s_open + "Namelist file: " + namelist + s_close)
;        system(s_open + "Creation date: `date`" + s_close)
;        system(s_open + "" + s_close)
;        system(s_open + "Please acknowledge the use of the ESMValTool." + \
;               s_close)
;        system(s_open + "Please cite Righi et al., Geosci. Model Dev., " + \
;               "8, 733-768 doi:10.5194/gmd-8-733-2015, 2015." + \
;               s_close)
;        system(s_open + "For the specific diagnostics, see below." + s_close)
;        system(s_open + "" + s_close)
;        system(s_open + hline + s_close)
        write_header(s_open, hline, s_close)
;        gOldDiag = script
; a_laue_ax-
    end if

    ;; write diagnostic header
    write_diag_header(s_open, hline_dble, s_close,  \
                      master_refs, output_refs,     \
                      auth, contr, diag, obs, proj, \
                      script)

;    ;; Add script name
;    system(s_open + hline_dble + s_close)
;    system(s_open + "=== " + script + " ===" + s_close)
;    system("echo " + " >> " + output_refs)
;
;    ;; Read master refs file
;    data = readAsciiTable(master_refs, 1, "string", 0)
;
;    ;; Author list
;    if (all(auth.ne."")) then
;        system(s_open + "AUTHOR(S): " + s_close)
;        do ii = 0, dimsizes(auth) - 1
;            pos = data(ind(str_squeeze(\
;                str_get_field(data(:, 0), 1, "=")) .eq. auth(ii)), 0)
;            str = str_squeeze(str_get_field(pos, 2, "="))
;            system(s_open + " -A- " + str + s_close)
;        end do
;        system("echo " + " >> " + output_refs)
;    end if
;
;    ;; Contributor list
;    if (all(contr.ne."")) then
;        system(s_open + "CONTRIBUTOR(S): " + s_close)
;        do ii = 0, dimsizes(contr) - 1
;            pos = data(ind(str_squeeze(\
;                       str_get_field(data(:, 0), 1, "=")) .eq. contr(ii)), 0)
;            str = str_squeeze(str_get_field(pos, 2, "="))
;            system(s_open + " -C- " + str + s_close)
;        end do
;        system("echo " + " >> " + output_refs)
;    end if
;
;    ;; Diagnostic reference list
;    if (all(diag.ne."")) then
;        system(s_open + "REFERENCE(S) FOR THIS DIAGNOSTIC: "\
;                      + s_close)
;        do ii = 0, dimsizes(diag) - 1
;            pos = data(ind(str_squeeze(\
;                       str_get_field(data(:, 0), 1, "=")) .eq. diag(ii)), 0)
;            str = str_squeeze(str_get_field(pos, 2, "="))
;            if (all(diag(ii).ne.stand_alone_entries)) then
;                str = "Please cite " + str
;            end if
;            system(s_open + " -R- " + str + s_close)
;        end do
;        system("echo " + " >> " + output_refs)
;    end if
;
;    ;; Observations reference list
;    if (all(obs.ne."")) then
;        system(s_open + "REFERENCE(S) FOR THE OBSERVATIONS: "\
;                      + s_close)
;        do ii = 0, dimsizes(obs) - 1
;            pos = data(ind(str_squeeze(\
;                       str_get_field(data(:, 0), 1, "=")) .eq. obs(ii)), 0)
;            str = str_squeeze(str_get_field(pos, 2, "="))
;            system(s_open + " -R- " + str + s_close)
;        end do
;        system("echo " + " >> " + output_refs)
;    end if
;
;    ;; Projects reference list
;    if (all(proj.ne."")) then
;        system(s_open + "ACKNOWLEDGEMENTS FOR THE PROJECTS: "\
;                      + s_close)
;        do ii = 0, dimsizes(proj) - 1
;            pos = data(ind(str_squeeze(\
;                       str_get_field(data(:, 0), 1, "=")) .eq. proj(ii)), 0)
;            str = str_squeeze(str_get_field(pos, 2, "="))
;            system(s_open + " -P- " + str + s_close)
;        end do
;        system("echo " + " >> " + output_refs)
;    end if
;
;    system(s_open + hline + s_close)

    leave_msg(scriptname, funcname, 4)

end


; A-laue_ax+
;;#############################################################################
undef("inquire_and_save_fileinfo")
procedure inquire_and_save_fileinfo(filename:string,
                                    model:string,
                                    variable:string)
;;
;; Arguments
;;    filename: file name of input file (netCDF)
;;    model: name of model
;;    variable: name of variable
;;
;; Description
;;    Read attributes from given file and write info to the log file
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151028-A_laue_ax: written.
;;
begin
    fid = addfile(filename, "r")

    filelist = new(10000, string)

    if (isatt(fid, "fixfile")) then
        fixfile = fid@fixfile
    else
        fixfile = "none"
    end if

    if (isatt(fid, "version")) then
        version = fid@version
    else
        version = "unknown"
    end if

    delete(fid)

    filelist(0) = filename
    n = 1
    ns = 0
    nr = 0

    if (n.eq.0) then
        filelist(0)    = "not specified"  ; list of files read
        n = 1
        sourcefilelist = ""  ; FIXFILE applied to original source files
        sourcefixfile  = ""  ; FIXFILE
        sourcereflist  = ""  ; list of references
        sourcetidlist  = ""  ; list of tracking_ids
    else
        ; now we check whether the files read (i.e., all names in "filelist")
        ; contain attribute(s) "infile_????". If so, we read these attributes
        ; and store them in "sourcefilelist". sourcefilelist then contains
        ; a list of original source files from which the files in "filelist"
        ; have originally been created, e.g., the original CMIP5 files or the
        ; original observational files.

        sourcefilelist    = new(dimsizes(filelist), string)
        sourcefilelist(0) = ""
        sourcereflist     = new(dimsizes(filelist), string)
        sourcereflist(0)  = ""
        sourcefixfile     = "none"
        sourcetidlist     = new(dimsizes(filelist), string)

        do m = 0, n - 1
            if (isfilepresent(filelist(m))) then
                sfile = addfile(filelist(m), "r")

                if (isatt(sfile, "fixfile")) then
                    sourcefixfile = sfile@fixfile
                end if

                ns = 0
                do while (ns.lt.dimsizes(sourcefilelist))
                    att = sprinti("infile_%04d", ns)
                    if (isatt(sfile, att)) then
                        sourcefilelist(ns) = sfile@$att$
                        ns = ns + 1
                    else
                        break
                    end if
                end do

                delete(sfile)
            end if
        end do

        ; Check for possible references and/or tracking_ids in the original
        ; source files.

        if (ns.gt.0) then
            ns = ns - 1

            do i = 0, ns
                if (isfilepresent(sourcefilelist(i))) then
                    sfile = addfile(sourcefilelist(i), "r")

                    ; check for tracking_id

                    if (isatt(sfile, "tracking_id")) then
                        sourcetidlist(i) = sfile@tracking_id
                    end if

                    ; check for references

                    if (isatt(sfile, "reference")) then
                        ; check if reference is already in list;
                        ; if reference is already in list, do not add again
                        newref = True
                        do j = 0, nr - 1
                            if (sourcereflist(j).eq.sfile@reference) then
                                newref = False
                                break
                            end if
                        end do
                        ; reference was not found in list --> add to list
                        if (newref) then
                            sourcereflist(nr) = sfile@reference
                            nr = nr + 1
                        end if
                    end if

                    ; if the source files contain derived variables they might contain
                    ; lists of source files as well
                    ; ---> check for source files of source files and (if any found)
                    ;      add them to the list of files to be written to the log-file

                    nss = 1
                    do while ((ns + nss).lt.dimsizes(sourcefilelist))
                        att = sprinti("infile_%04d", nss - 1)
                        if (isatt(sfile, att)) then
                             sourcefilelist(ns + nss) = sfile@$att$
                             nss = nss + 1
                        else
                            break
                        end if
                    end do

                    delete(sfile)
                    ns = ns + nss - 1

                end if
            end do
            if (nr.gt.0) then
               nr = nr - 1
            end if

        end if  ; if ns > 0
    end if

    ; write information on processed files to log-file

    write_filelist(diag_script, variable, model, filelist(0:n-1), \
                   sourcefilelist(0:ns), fixfile, sourcefixfile,  \
                   sourcereflist(0:nr), sourcetidlist(0:ns), version)

    delete(filelist)
    delete(sourcefilelist)
    delete(sourcereflist)
    delete(sourcetidlist)
end

