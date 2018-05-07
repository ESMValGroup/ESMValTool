;;#############################################################################
undef("write_header")
procedure write_header(s_open:string,
                       hline:string,
                       s_close:string)
;;
;; Arguments
;;    s_open: 1st part of system command line to add text to log file
;;    hline: separator (horizontal line)
;;    s_close: 2nd part of system command line to add text to log file
;;
;; Description
;;    Writes header of the log file
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20160713-A_laue_ax: updated ESMValTool reference
;;    20150630-A_laue_ax: written.
;;
local funcname, scriptname, namelist, version
begin

    funcname = "write_header"
    scriptname = "interface_scripts/logging.ncl"

    ;; Namelist in use
    namelist = diag_script_info@namelist

    ;; ESMValTool version
    version = diag_script_info@version

    dq = "\" + str_get_dq()

    system(s_open + hline + s_close)
    system(s_open + "+++++++++++++ " + \
           "ESMVal TOOL REFERENCES and ACKNOWLEDGEMENTS LOG" + \
           " +++++++++++++" + s_close)
    system(s_open + hline + s_close)
    system(s_open + "" + s_close)
    system(s_open + "Namelist file: " + namelist + s_close)
    system(s_open + "Creation date: `date`" + s_close)
    system(s_open + "ESMValTool version: " + version + s_close)
    system(s_open + "Host name: " + getenv("HOST") + s_close)
    system(s_open + "User name: " + getenv("USER") + s_close)
    system(s_open + "" + s_close)

    system(s_open + "Licensed under the Apache License, Version 2.0 (the " \
           + dq + "License" + dq + "); you may" + s_close)
    system(s_open + "not use this file except in compliance with the " \
           + "License. You may obtain" + s_close)
    system(s_open + "a copy of the License at" + s_close)
    system(s_open + "" + s_close)
    system(s_open + "http://www.apache.org/licenses/LICENSE-2.0" + s_close)
    system(s_open + "" + s_close)
    system(s_open + "Unless required by applicable law or agreed to in " \
           + "writing, software" + s_close)
    system(s_open + "distributed under the License is distributed on an " \
           + dq + "AS IS" + dq + "BASIS," + s_close)
    system(s_open + "WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either " \
           + "express or implied." + s_close)
    system(s_open + "See the License for the specific language governing " \
           + "permissions and" + s_close)
    system(s_open + "limitations under the License." + s_close)

    system(s_open + "" + s_close)
    system(s_open + "Please acknowledge the use of the ESMValTool." + s_close)
;    system(s_open + "Please cite Righi et al., Geosci. Model Dev., " + \
;           "8, 733-768 doi:10.5194/gmd-8-733-2015, 2015." + s_close)
    system(s_open + "Please cite Eyring et al., Geosci. Model Dev., " + \
           "9, 1747-1802, doi:10.5194/gmd-9-1747-2016, 2016." + s_close)
    system(s_open + "For the specific diagnostics, see below." + s_close)
    system(s_open + "" + s_close)
;    system(s_open + hline + s_close)
end

;;#############################################################################
undef("write_diag_header")
procedure write_diag_header(s_open:string, hline:string, s_close:string,
                            master_refs:string, output_refs:string,
                            auth:string, contr:string, diag:string, obs:string,
                            proj:string, script:string)
;;
;; Arguments
;;    s_open:      1st part of system command line to add text to log file
;;    hline:       separator (horizontal line)
;;    s_close:     2nd part of system command line to add text to log file
;;    master_ref:  file name of the master reference file
;;                 (i.e. MASTER_authors-refs-acknow.txt)
;;    output_refs: filename of log file
;;    auth:        list of authors
;;    contr:       list of contributors
;;    diag:        list of references for diagnostic
;;    obs:         list of references for observations
;;    proj:        list of project acknowledgements
;;    script:      name of diagnostic script
;;
;; Description
;;    Writes diagnostic header section to the log file
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151028-A_laue_ax: written.
;;
begin
    ;; D*** entries that shall not be extended by "Please cite "
    stand_alone_entries = (/"D000","D003","D004"/)

    ;; Add script name
    system(s_open + hline + s_close)
    system(s_open + "=== " + script + " ===" + s_close)
    system("echo " + " >> " + output_refs)

    ;; Read master refs file
    data = readAsciiTable(master_refs, 1, "string", 0)

    ;; Author list
    if (all(auth.ne."")) then
        system(s_open + "AUTHOR(S): " + s_close)
        do ii = 0, dimsizes(auth) - 1
            pos = data(ind(str_squeeze(\
                str_get_field(data(:, 0), 1, "=")) .eq. auth(ii)), 0)
            str = str_squeeze(str_get_field(pos, 2, "="))
            system(s_open + " -A- " + str + s_close)
        end do
        system("echo " + " >> " + output_refs)
    end if

    ;; Contributor list
    if (all(contr.ne."")) then
        system(s_open + "CONTRIBUTOR(S): " + s_close)
        do ii = 0, dimsizes(contr) - 1
            pos = data(ind(str_squeeze(\
                       str_get_field(data(:, 0), 1, "=")) .eq. contr(ii)), 0)
            str = str_squeeze(str_get_field(pos, 2, "="))
            system(s_open + " -C- " + str + s_close)
        end do
        system("echo " + " >> " + output_refs)
    end if

    ;; Diagnostic reference list
    if (all(diag.ne."")) then
        system(s_open + "REFERENCE(S) FOR THIS DIAGNOSTIC: "\
                      + s_close)
        do ii = 0, dimsizes(diag) - 1
            pos = data(ind(str_squeeze(\
                       str_get_field(data(:, 0), 1, "=")) .eq. diag(ii)), 0)
            str = str_squeeze(str_get_field(pos, 2, "="))
            if (all(diag(ii).ne.stand_alone_entries)) then
                str = "Please cite " + str
            end if
            system(s_open + " -R- " + str + s_close)
        end do
        system("echo " + " >> " + output_refs)
    end if

    ;; Observations reference list
    if (all(obs.ne."")) then
        system(s_open + "REFERENCE(S) FOR THE OBSERVATIONS: "\
                      + s_close)
        do ii = 0, dimsizes(obs) - 1
            pos = data(ind(str_squeeze(\
                       str_get_field(data(:, 0), 1, "=")) .eq. obs(ii)), 0)
            str = str_squeeze(str_get_field(pos, 2, "="))
            system(s_open + " -R- " + str + s_close)
        end do
        system("echo " + " >> " + output_refs)
    end if

    ;; Projects reference list
    if (all(proj.ne."")) then
        system(s_open + "ACKNOWLEDGEMENTS FOR THE PROJECTS: "\
                      + s_close)
        do ii = 0, dimsizes(proj) - 1
            pos = data(ind(str_squeeze(\
                       str_get_field(data(:, 0), 1, "=")) .eq. proj(ii)), 0)
            str = str_squeeze(str_get_field(pos, 2, "="))
            system(s_open + " -P- " + str + s_close)
        end do
        system("echo " + " >> " + output_refs)
    end if
end
