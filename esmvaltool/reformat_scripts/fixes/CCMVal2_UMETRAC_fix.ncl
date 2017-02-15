;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR, Germany)                                 #
;#######################################################################
;# PROJECT:     CCMVal2                                                #
;# ENSEMBLE:    REF-B1                                                 #
;# MODEL:       UMETRAC                                                #
;# VARIABLES:   ua                                                     #
;# FIELDS:      T3M                                                    #
;# DESCRIPTION: missing entries in time coordinate                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local opt,y1,y2,yy,mm,ii
begin

  ret = 1

  if (name.eq."ua".and.FIELD.eq."T3M") then
    opt = 0
    opt@calendar = var&time@calendar
    yy = 1960
    mm = 2
    do ii=0,dimsizes(var&time)-1
      var&time(ii) = tofloat(cd_inv_calendar(yy, mm, 15, 0, 0, 0,\
                                 var&time@units, opt))
      mm = mm+1
      if (mm.eq.13) then
        mm = 1
        yy = yy+1
      end if
    end do
    var@units = "m s-1"
    ret = 0
  end if

  return(ret)

end
