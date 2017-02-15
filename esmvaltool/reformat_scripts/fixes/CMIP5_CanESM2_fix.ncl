;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR-Germany)                                  #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       CanESM2                                                #
;# VARIABLES:   fgco2                                                  #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: Inconsistent units                                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

  ret = 1

  if (name.eq."fgco2") then
      var = var * 12. / 44.
      ret = 0
  end if

  return(ret)

end
