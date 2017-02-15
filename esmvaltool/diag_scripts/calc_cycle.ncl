load "./diag_scripts/calc_seasonal.ncl"

; calc_cycle.ncl
; data(time[0:11],lat, lon) => data(time[0:12],lat)
; data(time[0:11],lat)     => data(time[0:12],lat)

undef("calc_cycle")
function calc_cycle(data)
local ds, rank, data_zonal, data_cycle
begin
    ds = dimsizes(data)
    rank = dimsizes(ds)

    if ((rank .lt. 2) .or. (rank .gt. 3)) then
        print("calc_cycle: data should be 2D/3D")
        exit
    end if

    if ((rank .eq. 3) .and. (isdim(data, "lon"))) then
        data_zonal = dim_avg_Wrap(data)
    else
        data_zonal = data
    end if

    if (ds(0) .gt. 12) then
        data_mod = clmMonTL(data_zonal)
    else
        data_mod = data_zonal
    end if

    data_cycle = new((/13, ds(1)/), float)
    data_cycle!0 = "time"
    data_cycle&time = ispan(1, 13, 1)
    data_cycle!1 = "lat"
    data_cycle&lat = data&lat

    copy_VarAtts(data_zonal, data_cycle)

    data_cycle(0:11, :) = (/data_mod/)
    data_cycle(12, :) = (/data_cycle(0, :)/)

    return(data_cycle)
end
