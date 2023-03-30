; #############################################################################
; ROUTINES FOR CALCULATING GROWTH RATE AND SEASONAL CYCLE AMPLITUDE
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function calc_gr
;    function calc_sca
;    function var_detrend
;
; #############################################################################
undef("calc_gr")
function calc_gr(data[*]: numeric, \
                 opt[1]:string, \
                 min_nmonth: integer)
;
; Description:
;   Calculates growth rate of a 1D array using the Buchwitz et al 2018 method,
;   with a maximum of min_nmonth permitted missing values per year.
;
; Return value:
;
;
; Arguments:
;     data: requires time coordinate
;     opt:  "monthly"       - monthly GR
;           "yearly"        - yearly avg GR
;           "monthlyfilled" - monthly, missing values filled with yearly
;                             average
;     min_nmonth: maximum amount of missing values per year
;
;
; Modification history:
;    20180813_gier_bettina: written
;
local result, ntime, grmon, gryr, data1D, gryrf, mvnum
begin
    data1D = data

    ; Calc 1D case
    ntime = dimsizes(data1D)
    grmon = new(ntime, float)
    grmon(6:ntime-7) = data1D(12:) - data1D(:ntime-13)

    if opt.ne."monthly" then
        gryr = new(ntime/12, float)

        do itim = 0, ntime/12 - 1
            mvnum = num(ismissing(grmon(itim*12:itim*12+11)))
            if mvnum .le. min_nmonth then
                gryr(itim) = avg(grmon(itim*12:itim*12+11))
            end if
        end do

        if opt.eq."yearly" then
            result = gryr
        else if opt.eq."monthlyfilled" then
            gryrf = new(ntime, float)
            do itim = 0, ntime/12 - 1
                gryrf(itim*12:itim*12+11) = \
                              where(ismissing(grmon(itim*12:itim*12+11)), \
                                    gryr(itim),\
                                    grmon(itim*12:itim*12+11))
            end do
            result = gryrf
        end if
        end if

    else
        result = grmon
    end if

    return result
end

; #############################################################################
undef("calc_sca")
function calc_sca(data[*]: numeric, \
                  min_nmonth: integer)
;
; Description:
;   Calculates seasonal cycle amplitude by max - min of a 1D array, with a
;   maximum of min_nmonth permitted missing values per year.
;
; Return value:
;
;
; Arguments:
;     data:       1D array with data timeseries
;     min_nmonth: maximum amount of missing values per year
;
;
; Modification history:
;    20180813_gier_bettina: written
;
local sca, detrended, gryrf, ntime, mvnum, data1D
begin
    ; Get GR to detrend
    data1D = data
    ntime = dimsizes(data)
    sca = new(ntime/12, float)
    gryrf = calc_gr(data, "monthlyfilled", min_nmonth)
    do itim = 0, ntime/12 - 1
        detrended = data(itim*12:itim*12+11) \
                    - dim_cumsum(gryrf(itim*12:itim*12+11)/12., 1)
        mvnum = num(ismissing(detrended))
        if mvnum .le. min_nmonth then
            sca(itim) = max(detrended) - min(detrended)
        end if
    end do

    return sca
end

; #############################################################################
undef("var_detrend")
function var_detrend(x[*]: numeric,
                     y[*]: numeric)
;
; Description:
;  detrends a variable y with respect to x by linear regression
;
; Return value:
;
;
; Arguments:
;     x: var to detrend against (usually time)
;     y: var to detrend
;
;
; Modification history:
;    20181207_gier_bettina: written
;
local rc, rc_line, y_detr
begin
    rc = regline(x, y)
    rc_line = rc * x + rc@yintercept
    y_detr = y - rc_line
    copy_VarMeta(y, y_detr)
    return y_detr
end
