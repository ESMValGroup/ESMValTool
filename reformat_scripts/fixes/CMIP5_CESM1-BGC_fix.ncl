;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    1pctCO2, esmFixClim1                                   #
;# MODEL:       CESM1-BGC                                              #
;# VARIABLES:   time, co2, nbp                                         #
;# FIELDS:      all                                                    #
;# DESCRIPTION: various                                                #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

local dNames, lat_name, lon_name, lat, lon, lat_index, fout, areafile, \
    g, fx, area
begin

    ret = 1

    if (name.eq."co2") then
        if (var_attr_exp(0).eq."esmFixClim1") then
            var = var * (28.966 / 44.)
            ret = 0
        end if
    end if

    if (name.eq."nbp") then
        var@_FillValue = 1.e+33
        var@missing_value = var@_FillValue
        ret = 0
    end if
    
    if (iscoord(var, "time")) then
        if (isatt(var&time,"units"))then 
            if (any(var&time@units.eq.(/"days since 0001-01-01 00:00:00", \
                                        "days since 1-01-01 00:00:00"/)))
                var&time@units ="days since 1850-01-01 00:00:00"
                ret = 0
            end if
        end if
    end if

    return(ret)

end
