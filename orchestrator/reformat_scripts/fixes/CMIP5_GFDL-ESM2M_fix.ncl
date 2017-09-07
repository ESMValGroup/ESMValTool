;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    1pctCO2, esmFixClim1                                   #
;# MODEL:       GFDL-ESM2M                                             #
;# VARIABLES:   time                                                   #
;# FIELDS:      all                                                    #
;# DESCRIPTION: change time units                                      #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

    ret = 1

    if (name.eq."co2") then
        var = var * 1.e6
        ret = 0
    end if

    if (name.eq."sftof") then
        var = var * 100.
        ret = 0
    end if

    if (iscoord(var, "time")) then
        if (isatt(var&time,"units"))then 
            if (var&time@units.eq."days since 0001-01-01 00:00:00" .or. \
                var&time@units.eq."days since 1-01-01 00:00:00") then 
                var&time@units = "days since 1850-01-01 00:00:00"
                ret = 0
            end if
        end if
    end if

    return(ret)

end

