;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR, Germany)                                 #
;#######################################################################
;# PROJECT:     CCMVal2                                                #
;# ENSEMBLE:    REF-B1                                                 #
;# MODEL:       WACCM                                                  #
;# VARIABLES:   [all] except toz                                       #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: time coordinate must be shifted 2 days back            #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

  ret = 1

  if (name.ne."toz") then
    var&time = var&time - 2.

    ret = 0
  end if

  return(ret)

end

