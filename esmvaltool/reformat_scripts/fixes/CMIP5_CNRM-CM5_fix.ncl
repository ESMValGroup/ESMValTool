;########################################################################
;# Model-specific fixes for reformatDEFAULT_main                        #
;# Author: Axel Lauer (DLR, Germany)                                    #
;########################################################################
;# PROJECT:     CMIP5_ETHZ                                              #
;# MODEL:       CNRM-CM5                                                #
;# VARIABLES:   msftmyz, msftmyzba                                      #
;# DESCRIPTION: wrong scaling factor                                    #
;#              (see http://www.cnrm.meteo.fr/cmip5/spip.php?article24) #
;########################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

    ret = 1

    if (name.eq."msftmyz".or.name.eq."msftmyzba") then
        var = var * 1.e6
        ret = 0
    end if

    return(ret)

end
