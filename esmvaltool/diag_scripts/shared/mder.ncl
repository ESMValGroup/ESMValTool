; #############################################################################
; mder.ncl
; #############################################################################
; AUXILIARY FUNCTIONS FOR MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; #############################################################################
; Please consider using of extending existing routines before adding new
; ones. Check the header of each routine for documentation.
;
; Contents:
;
;   * function diag_file_name
;
; #############################################################################

load "interface_scripts/interface.ncl"


undef("diag_file_name")
function diag_file_name(dataset_idx[1]: integer, \
                        opt[1]: logical)
;
; Arguments:
;   dataset_idx : index of the current model
;   opt : options, needs the attributes
;                            @wdiag        (diagnostic name)
;                            @period_avg   (year-, month- or season name)
;                            @level        (specify plev range)
;                            @region       (specify lat and/or lon range)

; Return value:
;   file_name : returns the file name of the calculated climatological
;               diagnostic inclusive directory
; Description:
;   * provides a file name for the calculation results of diagsnostics.
;
; Modification history:
;   * 20180717-A_schl_ma: ported to v2.0
;   * 20140723-A_wenz_sa: written
;

local outname

begin

  enter_msg("mder.ncl", "")

  ; File name
  outname = basename(systemfunc("basename " + \
                                input_file_info@filename(dataset_idx)))

  ; Add extra info to the file
  if (isatt(opt, "wdiag")) then
    outname = outname + "_" + opt@wdiag
  end if
  if (isatt(opt, "period_avg")) then
    if (dimsizes(opt@period_avg) .eq. 2) then
      outname = outname + "_" + str_sub_str(opt@period_avg(0) + \
                                            opt@period_avg(1), " ", "")
    else
      outname = outname + "_" + opt@period_avg
    end if
  end if
  if (isatt(opt, "level")) then
    if (dimsizes(opt@level) .eq. 2) then
      outname = outname + "_" + str_sub_str(opt@level(0) + \
                                            opt@level(1), " ", "")
    else
      outname = outname + "_" + str_sub_str(opt@level, " ", "")
    end if
  end if
  if (isatt(opt, "region")) then
    if (dimsizes(opt@region) .eq. 4) then
      outname = outname + "_" + str_sub_str(opt@region(0) + opt@region(1) + \
                                            opt@region(2) + opt@region(3), \
                                            " ", "")
    else if (dimsizes(opt@region) .eq. 2) then
      outname = outname + "_" + str_sub_str(opt@region(0) + opt@region(1), \
                                            " ", "")
    else
      outname = outname + "_" + str_sub_str(opt@region, " ", "")
    end if
    end if
  end if
  outname = outname + ".nc"

  ; Return file name
  leave_msg("mder.ncl", "")
  return(outname)

end
