;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR-Germany)                                  #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       EC-EARTH                                               #
;# VARIABLES:   sic                                                    #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: Inconsistent units                                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

  ret = 1

  if (name.eq."sic") then
      var = 100.*var
      ret = 0
  end if

  if (name.eq."sftlf") then
      var = 100. * var
      ret = 0
  end if

  return(ret)

end
