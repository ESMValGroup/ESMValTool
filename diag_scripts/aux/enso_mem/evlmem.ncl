;;#############################################################################
;; GENERAL ROUTINES FOR MAXIMUM ENTROPY METHOD (MEM)
;;#############################################################################
;;
;; Contents:
;;    function evlmem
;;
;;#############################################################################

load "./interface_scripts/messaging.ncl"

;;#############################################################################

undef("evlmem")
function evlmem(fdt:numeric, d[*]:numeric, m:integer, xms:numeric)
;;
;; Arguments
;;    fdt: frequency times the sampling interval
;;    d: linear prediction coefficient (dimension = m)
;;    m: number of linear prediction coefficients
;;    xms: mean square discrepancy
;;
;; Return value
;;    estimated power spectrum (maximum entropy method)
;;
;; Description
;;    estimate power spectrum as a fucntion of frequency times the sampling interval
;;    code adopted from "Numerical Recipes in C" (sect. 13.7)
;;
;; Caveats
;;
;; References
;;     Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery:
;;     Numerical recipes in C: The art of scientific computing (sect. 13.7),
;;     second edition, Camebridge University Press, ISBN 0-521-43108-5.
;;
;; Modification history
;;    20170315-A_laue_ax: written.
;;

local sumr, sumi, wtemp, wr, wi, theta, wpr, wpi, i

begin

    funcname = "evlmem"
    scriptname = "diag_scripts/aux/enso_mem/evlmem.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    sumr = 1.0d
    sumi = 0.0d

    wr = 1.0d
    wi = 0.0d

    theta = 6.28318530717959d * fdt
    wpr = cos(theta)
    wpi = sin(theta)

    do i = 0, m - 1
        wtemp = wr
        wr = wr * wpr - wi * wpi
        wi = wi * wpr + wtemp * wpi
        sumr = sumr - d(i) * wr
        sumi = sumi - d(i) * wi
    end do

    leave_msg(scriptname, funcname, 6)

    return(xms / (sumr * sumr + sumi * sumi))

end
