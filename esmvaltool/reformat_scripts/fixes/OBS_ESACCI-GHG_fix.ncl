;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     OBS                                                    #
;# ENSEMBLE:                                                           #
;# MODEL:       ESACCI-GHG                                             #
;# VARIABLES:   xco2Std*, xch4Std*                                     #
;# FIELDS:      all                                                    #
;# DESCRIPTION: change units                                           #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

local dNames, lat_name, lon_name, lat, lon, lat_index, fout, areafile,g,fx,area
begin

    ret = 1
    
    ;if (name .eq. (/"xco2"/)) then
    ;    var = where(var .gt. 0.000402, var@_FillVallue, var)
        ;do tt = 0, dimsizes(var&time)-1
        ;    do la = 1, dimsizes(var&lat)-1
        ;        do lo = 1, dimsizes(var&lon)-1
        ;            if (.not.ismissing(var(tt,la,lo))) then
        ;                var(tt,la,lo) = where(var(tt,la,lo) .gt. \
        ;                                      var(tt,la-1,lo-1)*1.1, /
        ;                                      var@_FillVallue, \
        ;                                      var(tt,la,lo))
        ;            end if
        ;        end do
        ;    end do
        ;end do
    ;        ret = 0
    ;end if    

    if (any(name .eq. (/"xco2Stderr","xco2Stddev"/))) then
        if (isatt(var,"units")) then
            var = var *1.e6
            var@units = "1e-6"
            ret = 0
        end if
    end if

    if (any(name .eq. (/"xch4Stderr","xch4Stddev"/))) then
        if (isatt(var,"units")) then
            var = var *1.e9
            var@units = "1e-9"
            ret = 0
        end if
    end if

    return(ret)

end

