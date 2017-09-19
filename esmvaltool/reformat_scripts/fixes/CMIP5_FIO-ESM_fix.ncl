;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR-Germany)                                  #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       FIO-ESM                                                #
;# VARIABLES:   ch4                                                    #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: inconsistent units (ch4,co2)                           #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local f, files, tmp
begin

  ret = 1

  ;; Converting ch4 units from kg/kg to ppmb
  if (name.eq."ch4") then
      var =  var * 29/16 * 1.e9
      ret = 0
  end if

  ;; Converting co2 units from kg/kg to ppmb
  if (name.eq."co2") then
      var =  var * 29/44 * 1.e6
      ret = 0
  end if

  return(ret)

end
