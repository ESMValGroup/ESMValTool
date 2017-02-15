;#######################################################################
; This procedure reads references information from the master file and
; writes it in the reference/acknowledgement output file.
;#######################################################################
load "./interface_scripts/info_output.ncl"
load "./interface_scripts/auxiliary.ncl"

; Loading temporary files. These are updated on the fly by Python and
; used for Python->NCL communication
load "./interface_data/ncl.interface"

load "./diag_scripts/lib/ncl/set_operators.ncl"

undef("write_references")
procedure write_references(script:string, \  ; name of the script
                           auth:string,   \  ; author info code
                           contr:string,  \  ; contributor info code
                           diag:string,   \  ; diagnostic info code
                           obs:string,    \  ; observations info code
                           proj:string    \  ; projects info code
                           )
local namelist, master_refs, output_refs, hline, pos, str, stand_alone_entries
begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<<<<<<<<<<<<<<<<< Entering write_references.ncl", \
                verbosity, 2)
    hline = "-------------------------"
    hline = hline + hline + hline
  
    ;; Namelist in use
    namelist = xml  ; From ncl.interface

    ;; Master refs file
    master_refs = in_refs  ; From ncl.interface

    ;; Output refs file
    output_refs = out_refs  ; From ncl.interface

    s_open  = "echo " + str_get_dq
    s_close = str_get_dq + " >> " + output_refs

    ;; If first time (empty output_refs) write header
    if (toint(systemfunc("wc -c < " + output_refs)) .eq. 0) then
        system(s_open + hline + s_close)
        system(s_open + "+++++++++++++ "\
                      + "ESMValTool REFERENCES and ACKNOWLEDGEMENTS LOG"\
                      + " +++++++++++++" \
                      + s_close)
        system(s_open + hline + s_close)
        system(s_open + "" + s_close)
        system(s_open + "Namelist file: " + namelist + s_close)
        system(s_open + "Creation date: `date`" + s_close)
        system(s_open + "" + s_close)
        system(s_open + "Please acknowledge the use of the ESMValTool." + s_close)
        system(s_open + "Please cite Righi et al., Geosci. Model Dev., " \
                      + "8, 733-768 doi:10.5194/gmd-8-733-2015, 2015." \
                      + s_close)
        system(s_open + "For the specific diagnostics, see below." + s_close)
        system(s_open + "" + s_close)
        system(s_open + hline + s_close)
    end if

    ;; Read master refs file
    data = readAsciiTable(master_refs, 1, "string", 0)

    ;; Add script name
    system(s_open + "=== " + script + " ===" + s_close)
    system("echo " + " >> " + output_refs)

    ;; Author list
    if (all(auth.ne."")) then
        system(s_open + "AUTHOR(S): " + s_close)
        do ii = 0, dimsizes(auth) - 1
            pos = data(ind(str_squeeze(str_get_field(data(:, 0), 1, "=")) .eq. auth(ii)), 0)
            str = str_squeeze(str_get_field(pos, 2, "="))
            system(s_open + "  " + str + s_close)
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
            system(s_open + "  " + str + s_close)
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
            system(s_open + "  " + str + s_close)
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
            system(s_open + "  " + str + s_close)
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
            system(s_open + "  " + str + s_close)
        end do
        system("echo " + " >> " + output_refs)
    end if

    system(s_open + hline + s_close)
    info_output(">>>>>>>>>>>>>>>>>>>>>>> Leaving write_references.ncl", \
                verbosity, 2)

end
