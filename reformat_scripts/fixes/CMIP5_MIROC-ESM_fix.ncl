;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    1pctCO2, esmFixClim1,esmHistorical                     #
;# MODEL:       MIROC-ESM                                              #
;# VARIABLES:   time, co2                                               #
;# FIELDS:      all                                                    #
;# DESCRIPTION: change time units and calendar                         #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

local dNames, lat_name, lon_name, lat, lon, lat_index, fout, areafile,g,fx,area
begin

    ret = 1

    if (iscoord(var, "time")) then
        if (isatt(var&time,"units"))then 
            if (var&time@units.eq."days since 0000-01-01 00:00:00") then
                var&time@units ="days since 1849-01-01 00:00:00"
                ret = 0
            end if
            if (var&time@units.eq."days since 1-1-1")then
                var&time@units ="days since 1850-01-01 00:00:00"
                ret = 0
            end if
        end if
    end if

    if (iscoord(var, "AR5PL35")) then
        var!1 = "plev"
        ret = 0
    end if

    if (name .eq. "co2") then
        if (isatt(var,"units")) then
            var@units = "1.e-6"
            ret = 0
        end if
    end if

    if (name .eq. "gpp") then
        if (isatt(var,"units")) then  ; [gC/m2/day]
            do y=0,dimsizes(var&time) - 1
                var(y, :, :) = var(y, :, :) / (10 * 3600.)
            end do   
            delete(y)
            var@units = "kg m-2 s-1"
            ret = 0
        end if
    end if

    return(ret)

end

