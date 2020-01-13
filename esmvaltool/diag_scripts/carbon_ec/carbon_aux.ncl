; #############################################################################
; # EXTERNAL SEA ICE FUNCTIONS
; #############################################################################
;  Please consider using of extending existing routines before adding new ones.
;  Check the header of each routine for documentation.
;
;  Contents:
;     function sea_ice_area
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/latlon.ncl"

; #############################################################################
undef("cumulate_carbon")
function cumulate_carbon(NBP[*][*][*][*]:numeric)

;  Arguments:
;     NBP            : carbon flux variable
;
;  Return value:
;     3-d array (time) containing monthly cumulated carbon fluxes
;
;  Description:
;     Calculatesn cumulated carbon fluxes from T3M fields
;
;  Caveats
;
;  Modification history:
;     20191118-zechlau_sabrina: written.
;
local funcname, scriptname, verbosity, idx_new, NBP_tmp, NBP
begin

  funcname = "cumulate_carbon"
  scriptname = "diag_scripts/carbon_ec/carbon_aux.ncl"
  enter_msg(scriptname, funcname)

    ;; Integrate
    idx_new = dimsizes(NBP)
    NBP_tmp = new((/idx_new(0), idx_new(1), idx_new(2)/), typeof(NBP))
    NBP_tmp(0, :, :) = 0.0
    do y = 0, dimsizes(NBP(:, 0, 0)) - 2
        NBP_tmp(y + 1, :, :) = NBP_tmp(y, :, :) + NBP(y, :, :)
    end do
    delete(y)
    NBP = (/NBP_tmp/)

    NBP!0 = "time"
    NBP&time = tmp&time
    NBP!1 = "lat"
    NBP&lat = tmp&lat
    NBP!2 = "lon"
    NBP&lon = tmp&lon
    NBP@units = variable_info@units
    NBP@_FillValue = tmp@_FillValue

  return(NBP)
end
