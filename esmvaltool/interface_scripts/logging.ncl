;;#############################################################################
;; GENERAL ROUTINES FOR HANDLING LOG OUTPUT
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure log_info
;;    procedure log_debug
;;    procedure enter_msg
;;    procedure leave_msg
;;    procedure error_msg
;;    procedure tool_stop
;;    procedure exit_if_missing_atts
;;    procedure write_filelist
;;    procedure write_references
;;    procedure inquire_and_save_fileinfo
;;
;;#############################################################################
MASTER_REFS = "../doc/MASTER_authors-refs-acknow.txt"

gOldVar = ""
gOldDiag = ""
load "./interface_scripts/write_header.ncl"

;; WARNING: no other interface_scripts shall be loaded here, to avoid nested
;;          loops (since they will try to load logging.ncl themselves).

;;#############################################################################
undef("log_info")
procedure log_info(output_string[*]:string)
;;
;; Arguments
;;    output_string: the text to be output as message on screen
;;
;; Description
;;    Write an info message to the log file. If the input is an array, each
;;    element will be written on different lines.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20180202-A_righ_ma: written to replace info_output.
;;
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

;;#############################################################################
undef("log_debug")
procedure log_debug(output_string[*]:string)
;;
;; Arguments
;;    output_string: the text to be output as message on screen
;;
;; Description
;;    Write a debug message to the log file (only if log_level = debug in
;;    config-user.yml). If the input is an array, each element will be
;;    written on different lines.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20180202-A_righ_ma: written to replace info_output.
;;
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

;;#############################################################################
undef("enter_msg")
procedure enter_msg(script[1]:string,
                    func[1]:string)
;;
;; Arguments
;;    script: name of the ncl script file from where this function is called.
;;    func: name of the function or procedure calling this function.
;;
;; Description
;;    Print an enter message to the log file (debug mode only).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20180202_A_righ_ma: removed required verbosity.
;;    20150319-A_righ_ma: written.
;;
begin

    if (func.eq."") then
        log_debug("<<< Entering " + script)
    else
        log_debug("<<< Entering " + func + " (" + script + ")")
    end if

end

;;#############################################################################
undef("leave_msg")
procedure leave_msg(script[1]:string,
                    func[1]:string)
;;
;; Arguments
;;    script: name of the ncl script file from where this function is called.
;;    func: name of the function or procedure calling this function.
;;
;; Description
;;    Print a leave message to the log file (debug mode only).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20180202_A_righ_ma: removed required verbosity.
;;    20150319-A_righ_ma: written.
;;
begin

    if (func.eq."") then
        log_debug(">>> Leaving " + script)
    else
        log_debug(">>> Leaving " + func + " (" + script + ")")
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
;;    Write an error/warning message to the log file.
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
    else if (any(type.eq.(/"w", "warning", "error"/))) then
        msg_level = "warning"
    else
        print("fatal: in error_msg (logging.ncl), " + \
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
        log_info(msg_level + ": in " + msg_origin + msg(ii))
    end do

    ;; Exit if fatal
    if (msg_level .eq. "fatal") then
        status_exit(1)
    end if

end

;;#############################################################################
undef("tool_stop")
procedure tool_stop()
;;
;; Arguments
;;
;; Description
;;    Forces the tool execution to stop at any point of an ncl script.
;;    For dubugging purposes only.
;;
;; Caveats
;;    This should not be used for the standard exit-on-error. See the function
;;    error_msg in logging.ncl instead.
;;
;; References
;;
;; Modification history
;;    20150521-A_righ_ma: written.
;;
begin

    print("fatal: tool stop forced by the user")
    status_exit(1)

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
;;    Check that the given variable contains the given list of attributes.
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
local funcname, scriptname
begin

    funcname = "exit_if_missing_atts"
    scriptname = "interface_scripts/logging.ncl"
    enter_msg(scriptname, funcname)

    if (.not.all(isatt(in, reqatts))) then
        error_msg("f", diag_script, "", "missing required " + \
                  "diag_script_info attribute: " + \
                  reqatts(ind(.not.isatt(in, reqatts))))
    end if

    leave_msg(scriptname, funcname)

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
;;                    produce preproc files then used to process a derived
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
;;    Save the following details to the log file
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
local funcname, scriptname, output_refs, hline
begin

    funcname = "write_filelist"
    scriptname = "interface_scripts/logging.ncl"
    enter_msg(scriptname, funcname)

    if (filelist(0).eq."") then  ; derived variable or no information
        return
    end if

    hline = "-------------------------"
    hline = hline + hline + hline

    ;; Output refs file
    output_refs = config_user_info@run_dir + "/references-acknowledgements.txt"

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

    leave_msg(scriptname, funcname)
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
;;    Write references and acknowledgments in the log file
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
local funcname, scriptname, master_refs, output_refs, hline, pos, str
begin

    funcname = "write_references"
    scriptname = "interface_scripts/logging.ncl"
    enter_msg(scriptname, funcname)

    hline      = "-------------------------"
    hline_dble = "========================="
    hline = hline + hline + hline
    hline_dble = hline_dble + hline_dble + hline_dble

;    ;; D*** entries that shall not be extended by "Please cite "
;    stand_alone_entries = (/"D000", "D003", "D004"/)

; a_laue_ax+
;    ;; Namelist in use
;    namelist = diag_script_info@namelist
; a_laue_ax-

    ;; Master refs file
    master_refs = MASTER_REFS

    ;; Output refs file
    output_refs = config_user_info@run_dir + "/references-acknowledgements.txt"

    s_open  = "echo " + str_get_dq
    s_close = str_get_dq + " >> " + output_refs

    ;; If first time (empty output_refs) write header
    if (.not. fileexists(output_refs)) then
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

    leave_msg(scriptname, funcname)

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

