;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR-Germany)                                  #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       inmcm4                                                 #
;# VARIABLES:   gpp                                                    #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: Inconsistent units                                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

  ret = 1

  if (name.eq."gpp") then
      var = -1 * var
      ret = 0
  end if

  if (name.eq."lai") then
      var = var / 100.
      ret = 0
  end if

  return(ret)

end
