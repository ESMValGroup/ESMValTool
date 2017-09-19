;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR, Germany)                                 #
;#######################################################################
;# PROJECT:     CCMVal2                                                #
;# ENSEMBLE:    REF-B1                                                 #
;# MODEL:       CAM3.5                                                 #
;# VARIABLES:   [all]                                                  #
;# FIELDS:      [all]
;# DESCRIPTION: time coordinate must be shifted 2 days back            #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local opt,y1,y2,yy,mm,ii
begin

  ret = 1
  var&time = var&time - 2.
  ret = 0

  if (name.eq."va".and.FIELD.eq."T3M") then
    opt = 0
    opt@calendar = var&time@calendar
    yy = 1960
    mm = 1
    do ii=0,dimsizes(var&time)-1
      var&time(ii) = tofloat(cd_inv_calendar(yy, mm, 15, 0, 0, 0,\
                                 var&time@units, opt))
      mm = mm+1
      if (mm.eq.13) then
        mm = 1
        yy = yy+1
      end if
    end do
  end if

  return(ret)

end

