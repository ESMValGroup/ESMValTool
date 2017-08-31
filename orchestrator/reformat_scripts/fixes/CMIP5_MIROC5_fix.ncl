;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Axel Lauer (DLR, Germany)                                   #
;#######################################################################
;# PROJECT:     CMIP5_ETHZ                                             #
;# ENSEMBLE:    r1i1p1                                                 #
;# MODEL:       MIROC5C                                                #
;# VARIABLES:   snc, snw                                               #
;# FIELDS:      T2Ds                                                   #
;# DESCRIPTION: missing entries in time coordinate                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local opt,y1,y2,yy,mm,dd,ii
begin

    ret = 1

    ; calendar = "noleap"
    dayspermonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    if ((name.eq."snc".or.name.eq."snw").and.FIELD.eq."T2Ds".and. \
        ENSEMBLE.eq."r1i1p1") then
        opt = 0
        opt@calendar = var&time@calendar
        ; get start date from time attribute "days_since_xxxx"
        t = 0.0
        t@calendar = var&time@calendar
        t@units = var&time@units
        res = cd_calendar(t, -5)
        yy = res(0, 0)
        mm = res(0, 1)
        dd = res(0, 2)
        do ii = 0, dimsizes(var&time) - 1
            var&time(ii) = tofloat(cd_inv_calendar(yy, mm, dd, 12, 0, 0, \
                                   var&time@units, opt))
            dd = dd + 1
            if (dd.gt.dayspermonth(mm-1)) then
                mm = mm + 1
                dd = 1
            end if
            if (mm.gt.12) then
                mm = 1
                yy = yy + 1
            end if
        end do
        ret = 0
    end if

    if (name.eq."sftof") then
        var = var * 100.
        ret = 0
    end if

    return(ret)

end
