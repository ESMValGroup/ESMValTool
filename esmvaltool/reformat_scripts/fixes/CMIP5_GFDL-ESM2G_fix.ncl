;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    esmHistorical                                          #
;# MODEL:       GFDL-ESM2M                                             #
;# VARIABLES:   time                                                   #
;# FIELDS:      all                                                    #
;# DESCRIPTION: inconsistent units (co2)                               #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

    ret = 1

  ;; Converting co2 units to ppmb
  if (name.eq."co2") then
      var =  var * 1.e6
      ret = 0
  end if

    return(ret)

end

