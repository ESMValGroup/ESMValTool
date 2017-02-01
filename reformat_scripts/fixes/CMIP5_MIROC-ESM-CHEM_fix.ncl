;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Sabrina Wenzel (DLR-Germany)                                #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    historical                                             #
;# MODEL:       MIROC-ESM-CHEM                                         #
;# VARIABLES:   time, tro3                                             #
;# FIELDS:      all                                                    #
;# DESCRIPTION: change time units and calendar                         #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

local dNames, lat_name, lon_name, lat, lon, lat_index, fout, areafile,g,fx,area
begin

    ret = 1

    if (name .eq. "tro3") then
        ; files say unit of ozone is "1e-9" ut unit is actually "1e-6"
        var = var * 1.0e3
        if (iscoord(var, "time")) then
;            tt = cd_calendar(var&time,0)
;            var&time@_FillValue = default_fillvalue(typeof(var&time))
;            do it = 1,dimsizes(tt(:,0))-1
;                if (.not.ismissing(tt(it-1,0)).and.tt(it-1,0).gt.tt(it,0)+1) then
;                    tmp_it = tt(it,0)
;                    var&time(it) = var&time@_FillValue
;                else if (ismissing(tt(it-1,0)).and.tt(it,0)+1.gt.tmp_it) then
;                    var&time(it) = var&time@_FillValue
;                end if
;                end if
;            end do
            do it = 1, dimsizes(var&time) - 1
                if (var&time(it).eq.0) then
                    tt = tointeger(cd_calendar(var&time(it-1), 0))
                    tt(0, 1) = tt(0, 1) + 1  ; month
                    if (tt(0, 1).gt.12) then
                        tt(0, 1) = 1
                        tt(0, 0) = tt(0, 0) + 1  ; year
                    end if
                    var&time(it) = cd_inv_calendar(tt(0, 0), tt(0, 1), tt(0, 2), tt(0, 3), \
                                                   tt(0, 4), tt(0, 5), var&time@units, 0)
                end if
            end do
        ret = 0   
        end if
    end if

    return(ret)

end

