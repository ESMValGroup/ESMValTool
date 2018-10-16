;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR-Germany)                                  #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    1pctCO2, esmFixClim1                                   #
;# MODEL:       GFDL-CM2p1                                             #
;# VARIABLES:   sftof                                                  #
;# FIELDS:      all                                                    #
;# DESCRIPTION: inconsistent units                                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

    ret = 1

    if (name.eq."sftof") then
        var = var * 100.
        ret = 0
    end if

    return(ret)

end

