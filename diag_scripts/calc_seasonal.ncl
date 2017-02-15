load "./interface_scripts/messaging.ncl"

; ******************************************************************
; Calculate long term monthly means (monthly climatology)
; requires named dimensions
;
; The time dimension must be a multiple of 12
;
;    x(time, lat)  <==== INPUT DIMENSION ORDER
;    x!0 = "time"
;    x!1 = "lat"
;
undef("clmMonTL")
function clmMonTL(x[*][*]:numeric)
local dimx, nlat, ntim, i, nmo, nmos, xLT, aveMonth
begin

    funcname = "clmMonTL"

    dimx = dimsizes(x)
    ntim = dimx(0)
    nlat = dimx(1)

    nmos  = 12
    modCheck("clmMonTL", ntim, nmos)   ; error check

    do i = 0, 1
        if (ismissing(x!i)) then
            error_msg("f", "calc_seasonal.ncl", funcname, "clmMonTL: all " + \
                      "dimensions must be named")
        end if
    end do

    xLT = x($x!1$|:, $x!0$|:)

    ;; Compute all 12 monthly averages first.
    aveMonth = new((/nlat, nmos/), typeof(x), getFillValue(x))

    do nmo = 0, nmos - 1
        aveMonth(:, nmo) = dim_avg(xLT(:, nmo:ntim - 1:nmos))
    end do

    ;; copy attributes and add a couple of informational attributes
    ;; only need to do this if plotting or writing to netCDF file
    copy_VarAtts(x, aveMonth)
    aveMonth@time_op_ncl = "Climatology: " + (ntim/nmos) + " years"

    ;; copy spatial (eg, lat) coordinate variables
    aveMonth!0 = xLT!0
    if (iscoord(xLT, xLT!0)) then
        aveMonth&$aveMonth!0$ = xLT&$xLT!0$
    end if

    aveMonth!1     = "month"           ; create a "month" named dim
    aveMonth&month = ispan(0, nmos - 1, 1)  ; create a month coord var

    return (aveMonth(month|:, $x!1$|:))   ; return (month, lat)
end
; **************************************************


; ==================================================
undef("calc_SEASON")
function calc_SEASON(data,
                     season : string)
;;       return_val dimsizes(data) - 1 : typeof(data)
;;
;; Arguments:
;;    data  --  one of 2D (time, lat), 3D (time, lat, lon) or 4D (time, lev, lat, lon)
;;     season  --  requested season in upper case (any subset of JFMAMJJASOND)
;; Return value:
;;    data_new -- average of climatology for given season
;; Description:
;;    Compute the climatology and then the average for the requested season
local data_new, ds, rank, dimNames, data_season, data_mod, idx_mm
begin
    ds = dimsizes(data)
    rank = dimsizes(ds)
    ;; Compute climatology
    if (ds(0) .gt. 12) then
        if (rank .eq. 2) then
            data_mod = clmMonTL(data)
        end if
        if (rank .eq. 3) then
            data_mod = clmMonTLL(data)
        end if
        if (rank .eq. 4) then
            data_mod = clmMonTLLL(data)
        end if
    else
        data_mod = data
    end if

    ;; Create return value array
    if (rank .eq. 2) then
        data_new = new((/ds(1)/), typeof(data))
    end if
    if (rank .eq. 3) then
        data_new = new((/ds(1), ds(2)/), typeof(data))
    end if
    if (rank .eq. 4) then
        data_new = new((/ds(1), ds(2), ds(3)/), typeof(data))
    end if
    copy_VarAtts(data, data_new)
    dimNames = getvardims(data)

    do i = 1, rank - 1, 1
        ip1 = i - 1
        data_new!ip1 = dimNames(i)
        data_new&$dimNames(i)$ = data&$dimNames(i)$
    end do

    ;; Compute average for the requested season
    if (season .eq. "ANN" .or. season .eq. "ANNUAL") then
        data_season = month_to_annual(data_mod, 1)
    else
        idx_mm = calc_season_index(season)

        if (rank .eq. 2) then
            data_season = data_mod(idx_mm(0), :)
            do idx = 1, dimsizes(idx_mm) - 1
                data_season = data_season + data_mod(idx_mm(idx), :)
            end do
        end if

        if (rank .eq. 3) then
            data_season = data_mod(idx_mm(0), :, :)
            do idx = 1, dimsizes(idx_mm) - 1
                data_season = data_season + data_mod(idx_mm(idx), :, :)
            end do
        end if

        if (rank .eq. 4) then
            data_season = data_mod(idx_mm(0), :, :, :)
            do idx = 1, dimsizes(idx_mm) - 1
                data_season = data_season + data_mod(idx_mm(idx), :, :, :)
            end do
        end if
        idx_len = int2flt(dimsizes(idx_mm))
        data_season = data_season / idx_len
    end if

    data_new = (/data_season/)
    copy_VarAtts(data, data_new)

    return(data_new)

end
; **************************************************

; ==================================================
undef("day_to_month")
function day_to_month(data, fno)
;; undocumented function
local data_new, t_new, tim, dt, ds, t_utc, nt, t_st_yr, t_ed_yr, dimNames, \
      year, month, idx_i, idx_ym, rank
begin
    ;; data: T2Ds(time, lat, lon), T2Dz(time, lev, lat), T1Ds(time, lat)
    ds = dimsizes(data)
    rank = dimsizes(ds)

    tim = data&time
    dt = dimsizes(tim)

    t_utc = cd_calendar(tim, 0)
    t_st_yr = floattoint(t_utc(0, 0))
    t_ed_yr = floattoint(t_utc(dimsizes(tim) - 1, 0))
    nt = (t_ed_yr - t_st_yr + 1) * 12

    if (rank .eq. 2) then
        data_new = new((/nt, ds(1)/), typeof(data))
    end if
    if (rank .eq. 3) then
        data_new = new((/nt, ds(1), ds(2)/), typeof(data))
    end if

    copy_VarAtts(data, data_new)
    dimNames = getvardims(data)
    t_new = new((/nt/), typeof(tim))
    t_new@units = tim@units

    do i = 1, rank - 1
        data_new!i = dimNames(i)
        data_new&$dimNames(i)$ = data&$dimNames(i)$
    end do

    do year = t_st_yr, t_ed_yr
        do month = 1, 12
            idx_i = (year - t_st_yr) * 12 + (month - 1)
            idx_ym = ind((t_utc(:, 0) .eq. year) .and. (t_utc(:, 1) .eq. month))
            if (.not. all(ismissing(idx_ym))) then
                if (rank .eq. 3) then
                    if (fno .eq. "T2Ds") then
                        data_new(idx_i, :, :) = (/dim_avg_Wrap(data(lat|:, lon|:, time|min(idx_ym):max(idx_ym)))/)
                    else
                        data_new(idx_i, :, :) = (/dim_avg_Wrap(data(plev|:, lat|:, time|min(idx_ym):max(idx_ym)))/)
                    end if
                else
                    data_new(idx_i, :) = (/dim_avg_Wrap(data(lat|:, time|min(idx_ym):max(idx_ym)))/)
                end if
                t_new(idx_i) = avg(tim(min(idx_ym):max(idx_ym)))
            else
                t_new(idx_i) = (/cd_inv_calendar(year, month, 15, 0, 0, 0, tim@units, 0)/)
            end if
            delete(idx_ym)
        end do
    end do

    data_new!0 = "time"
    data_new&time = t_new

    return(data_new)
end
; **************************************************

; ==================================================
undef("clim_season_average")
function clim_season_average(indata [*][*][*]:float,
                             wgt_lat [*]:float,
                             season  [1]:string)
local climdata, tmp
begin
    climdata = calc_SEASON(indata, season)

    ;; compute the weighted average
    wgt_lon = conform_dims((/dimsizes(indata&lon)/), 1.0, -1)
    tmp = wgt_areaave(climdata, wgt_lat, wgt_lon, 0 )
    return climdata
end
