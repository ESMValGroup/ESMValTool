;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    1pctCO2, esmFixClim1                                   #
;# MODEL:       MPI-ESM-LR                                             #
;# VARIABLES:   time, pctisccp                                         #
;# FIELDS:      all, T2Ds                                              #
;# DESCRIPTION: change time units and calendar                         #
;#              pctisccp: convert from hPa to Pa                       #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
begin

    ret = 1

    if (iscoord(var, "time")) then
        if (isatt(var&time,"units")) then 
            if (var&time@units.eq."days since 0001-01-01 00:00:00".or. \
                var&time@units.eq."days since 1-01-01 00:00:00")
                var&time@units ="days since 1850-01-01 00:00:00"
                ret = 0
            end if
        end if
        if (isatt(var&time,"calendar")) then
            if (var&time@calendar.eq."proleptic_gregorian") then
                var&time@calendar="standard"
            end if
            ret = 0
        end if
    end if

    if (name.eq."pctisccp".and.FIELD.eq."T2Ds") then
        var = 100.*var
        ret = 0
    end if

    return(ret)

end
