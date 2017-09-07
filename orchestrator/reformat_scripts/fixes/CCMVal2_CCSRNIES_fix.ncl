;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR, Germany)                                 #
;#######################################################################
;# PROJECT:     CCMVal2                                                #
;# ENSEMBLE:    REF-B1                                                 #
;# MODEL:       CCSRNIES                                               #
;# VARIABLES:   va, zg                                                 #
;# FIELDS:      T2Mz                                                   #
;# DESCRIPTION: time coordinate not properly defined                   #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local opt,y1,y2,yy,mm,ii
begin

  ret = 1

  if (any(name.eq.(/"va","zg"/)).and.FIELD.eq."T2Mz") then
    opt = 0
    opt@calendar = var&time@calendar
    y1 = 1960
    y2 = 1960 + dimsizes(var&time)/12 - 1
    do yy = y1,y2
      do mm = 1,12
        ii = 12*(yy-y1) + mm - 1
        var&time(ii) = cd_inv_calendar(yy, mm, 15, 0, 0, 0,\
                               var&time@units, opt)
      end do
    end do
    ret = 0
  end if

  return(ret)

end
