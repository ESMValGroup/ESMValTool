;_______________________________________________________________________________
;      function depth_conform( depth:numeric)
;_______________________________________________________________________________
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
undef("depth_conform")
function depth_conform( depth : numeric)
;;               return_val : [*] float
;; Arguments:
;;    depth  - Pressure level data in Pa or hPa
;; Return value:
;;    depth   - the same pressure level data but in units hPa
;; Description:
;;    Converts the input pressure level data from Pa (or hPa) to hPa
local no_operation, depth_return
begin
    if (depth@units .eq. "Km" .or. depth@units .eq. "km" ) then
        depth_return = depth_return * 1000.
    else 
        if ( depth@units .eq. "m" .or. depth@units .eq. "meters" .or. depth@units .eq. "metres" ) then
            depth_return = depth
        else
            print("fatal: depth_conform.ncl: Only m or Km are valid units for depths")
            status_exit(1)
        end if
    end if
    return depth_return
end
