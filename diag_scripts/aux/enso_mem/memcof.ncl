;;#############################################################################
;; GENERAL ROUTINES FOR MAXIMUM ENTROPY METHOD (MEM)
;;#############################################################################
;;
;; Contents:
;;    procedure memcof
;;
;;#############################################################################

load "./interface_scripts/messaging.ncl"

;;#############################################################################

undef("memcof")
procedure memcof(data:numeric, m:integer, xms:numeric, d:numeric)
;;
;; Arguments
;;    data: time series to analyze (input)
;;    m:    number of requested linear prediction coefficients (input)
;;    xms:  mean square discrepancy (output)
;;    d:    linear prediction coefficients, dimension m (output)
;;
;; Description
;;    estimation of the linear prediction (LP) coefficients d(j) of a data set
;;    code adopted from "Numerical Recipes in C" (sect. 13.6)
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

local i, k, n, d1, wk1, wk2, wkm, p, denom, xnum

begin

    funcname = "memcof"
    scriptname = "diag_scripts/aux/enso_mem/memcof.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    n = dimsizes(data)

    if (m .gt. n) then
        m = n
    end if

    ; note: arrays d1, wk1, wk2, wkm will be subscripted with indices 1...m (1...n)
    ;       while array data is subscripted with 0...n-1 (NCL standard)

    d1 = new((/m + 1/), float)

    wk1 = new((/n + 1/), float)
    wk2 = new((/n + 1/), float)
    wkm = new((/m + 1/), float)

    wkm = 0.0

    p = sum(data(:)^2)

    xms = p / tofloat(n)

    wk1(1) = data(0)
    wk2(n-1) = data(n - 1)

    wk1(2:n-1) = data(1:n-2)
    wk2(1:n-2) = data(1:n-2)

    do k = 1, m
        xnum = sum(wk1(1:n-k) * wk2(1:n-k))
        denom = sum(wk1(1:n-k)^2 + wk2(1:n-k)^2)
        d1(k) = 2.0 * xnum / denom
        xms = xms * (1.0 - d1(k)^2)

        do i = 1, k - 1
            d1(i) = wkm(i) - d1(k) * wkm(k-i)
        end do

        if (k .eq. m) then
            break
        end if

        wkm(1:k) = d1(1:k)

        do j = 1, n - k - 1
            wk1(j) = wk1(j) - wkm(k) * wk2(j)
            wk2(j) = wk2(j+1) - wkm(k) * wk1(j+1)
        end do

    end do  ; k-loop

    delete(wk1)
    delete(wk2)
    delete(wkm)

    d(:) = d1(1:m)
    delete(d1)

    leave_msg(scriptname, funcname, 6)

end
