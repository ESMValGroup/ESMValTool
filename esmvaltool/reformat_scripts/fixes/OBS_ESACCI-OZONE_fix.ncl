;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Name Surname (Affiliation)                                  #
;#######################################################################
;# PROJECT:                                                            #
;# ENSEMBLE:                                                           #
;# MODEL:       ESA CCI ozone profiles                                 #
;# VARIABLES:   tro3prof                                               #
;# FIELDS:                                                             #
;# DESCRIPTION: vertical coordinate says "Pa" but is actually "hPa"    #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

local ret, tmp, ti, tt
begin
    verbosity  = diag_script_info@verbosity

    ret = 1

    if (name .eq. "tro3prof") then
        if (isdim(var,"plev")) then
            var&plev = var&plev * 100.0  ; hPa --> Pa
            ret = 0
        end if
    end if

    ; The return value is 0 if the fix has been applied, and 1 otherwise
    return(ret)
end

