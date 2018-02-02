;#######################################################################
; This code reads references information from the master file and
; writes it to the reference/acknowledgement output file.
; This is a stand-alone version of the function "write_references"
; (interface_scripts/logging.ncl) called by Python to write the
; reference/acknowledgement output file.
;#######################################################################

;; definition of functions write_header and write_diag_header
load "./interface_scripts/write_header.ncl"


begin

    hline      = "-------------------------"
    hline_dble = "========================="
    hline = hline + hline + hline
    hline_dble = hline_dble + hline_dble + hline_dble

    ;; Master refs file "master_refs" and output refs file "output_refs" as
    ;; well as "auth", "contr", "diag", "obs", "proj" are defined in
    ;; "$ESMValTool_interface_data/ncl.ref.info"

    s_open  = "echo " + str_get_dq
    s_close = str_get_dq + " >> " + out_refs

    ;; If first time (empty output_refs) write header
    if (toint(systemfunc("wc -c < " + out_refs)) .eq. 0) then
        write_header(s_open, hline, s_close)
    end if

    ;; write diagnostic header section
    write_diag_header(s_open, hline_dble, s_close, in_refs, out_refs, \
                      ref_auth, ref_contr, ref_diag, ref_obs, ref_proj, \
                      ref_script)

end
