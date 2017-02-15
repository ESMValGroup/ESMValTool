;;#############################################################################
;; GENERAL ROUTINES FOR STATISTICS
;;#############################################################################
;; Please consider using or extending existing routines before adding new
;; ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;    function calc_wave_freq
;;    function calc_pr_ua_lag
;;    function calc_pr_ua_mean
;;    function calc_pr_ua_vari
;;    function calc_olr_ua_spectra
;;    function calc_daily_univari_eof
;;    function calc_daily_multi_vari_ceof
;;    function calc_daily_olr_uv_life_cycle
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
;-------------------------------
undef("calc_wave_freq")
function calc_wave_freq(nml_data_vault_pr[1]: logical,\
                                       di[1]: logical,\
                               idx_mod_pr[1]: numeric,\
                                   season[1]: string,\
                                 Ref_Read[1]: logical)
local mean_pr,       nFill,  no_of_harmonics,\
    optWavFrq,    pr,     prClmDay,\
    prClmDay_sm,  pr_sm,  pr_tmp,\
    transforms,   var0,   wfO,\
    ymd,          yyyyddd
begin
    transforms = True
    transforms@apply_these = (/"latrange"/)
    transforms@latrange = di@isv_latrange

    var0 = variables(0)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; entire file
    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if
    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3
    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    pr_sm = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    mean_pr = dim_avg_n_Wrap(pr_sm(:, :, :), 1)

    ;; replacing missing value 1e+20 with 0
    if (isatt(mean_pr, "_FillValue")) then
        nFill = num(ismissing(mean_pr))
        if (nFill.gt.0) then
            missv = mean_pr@_FillValue
            delete(mean_pr@_FillValue)
            mean_pr = where(mean_pr.eq.missv, 0, mean_pr)
            mean_pr@_FillValue = missv
            delete(missv)
        end if
    end if

    optWavFrq = False

    wfO = mjo_wavenum_freq_season(mean_pr, ymd, season, optWavFrq)

    delete(prClmDay)
    delete(prClmDay_sm)
    delete(pr_sm)
    delete(mean_pr)

    return(wfO)
end
;----------------------------------------
undef("calc_pr_ua_lag")
function calc_pr_ua_lag(nml_data_vault_pr[1]: logical,\
                        nml_data_vault_ua[1]: logical,\
                                       di[1]: logical,\
                               idx_mod_pr[1]: numeric,\
                               idx_mod_ua[1]: numeric,\
                                   season[1]: string)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local aver,         date_u,            fca,\
    fcb,          ihp,               lat,\
    latn,         latN_globe,        latN_IO,\
    lats,         latS_globe,        latS_IO,\
    lonl,         lonL_IO,           lonr,\
    lonR_IO,      missv,             mxlag,\
    nFill,        no_of_harmonics,   nWgt,\
    optXcor,      P,                 PIO,\
    pr,           prClmDay,          prClmDay_sm,\
    pr_sm,        pr_tmp,            P_timeLat,\
    P_timeLon,    rp_timelat,        rp_timelon,\
    ru_timelat,   ru_timelon,        sigma,\
    transforms,   U,                 ua,\
    uaClmDay,     uaClmDay_sm,       ua_sm,\
    ua_tmp,       U_timeLat,         U_timeLon,\
    var0,         var1,              wgt,\
    wyIO,         ymd,               yyyyddd
begin
    ;; Compute lancoz weights
    ihp = di@ihp     ; bpf = >band pass filter
    nWgt = di@nWgt
    sigma = di@sigma   ; Lanczos sigma
    fca = di@fca     ; start freq
    fcb = di@fcb     ; end freq
    wgt = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

    ;;
    ;; Region definition
    ;;

    latS_IO = di@eio_latrange(0)
    latN_IO = di@eio_latrange(1)
    lonL_IO = di@eio_lonrange(0)
    lonR_IO = di@eio_lonrange(1)

    latS_globe = di@global_latrange(0)        ; global subset [Fig 6]
    latN_globe = di@global_latrange(1)

    lats = di@prop_lag_latrange(0)
    latn = di@prop_lag_latrange(1)        ; lat band for (lag, lon) Fig 5

    lonl = di@prop_lag_lonrange(0)        ; lon band for (lag, lat) Fig 6
    lonr = di@prop_lag_lonrange(1)

    transforms = True
    transforms@apply_these = (/"latrange"/)
    transforms@latrange = di@global_latrange

    var0 = variables(0)
    var1 = variables(1)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)

    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; entire file

    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if

    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3
    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    pr_sm = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    ;; IO in 'wyIO'/pr_smIO refers to the reference area 'Indian Ocean'
    lat = pr&lat
    wyIO = lat(ind(lat .ge. latS_IO .and.\
                   lat .le. latN_IO))
    wyIO = cos(0.017459 * wyIO)         ; spatial wgt
    delete(lat)

    ua_tmp = extract_data_extended(idx_mod_ua,\
                                   var1,\
                                   nml_data_vault_ua,\
                                   -1,\
                                   0,\
                                   0)
    ua = apply_transforms(transforms, ua_tmp)
    uaClmDay = clmDayTLL_calendar(ua, yyyyddd)

    ;; replacing missing value 1e+20 with 0
    if (isatt(uaClmDay, "_FillValue")) then
        nFill = num(ismissing(uaClmDay))
        if (nFill.gt.0) then
            missv = uaClmDay@_FillValue
            aver = avg(uaClmDay)
            delete(uaClmDay@_FillValue)
            uaClmDay = where(uaClmDay.eq.missv, aver, uaClmDay)
            uaClmDay@_FillValue = missv
            delete(missv)
            delete(aver)
        end if
    end if

    uaClmDay_sm = smthClmDayTLL_calendar(uaClmDay, no_of_harmonics)

    ua_sm = calcDayAnomTLL_calendar(ua, yyyyddd, uaClmDay_sm)

    delete(ua_tmp)
    delete(pr_tmp)
    delete(pr)
    delete(ua)
    delete(prClmDay)
    delete(prClmDay_sm)
    delete(uaClmDay)
    delete(uaClmDay_sm)

    ;;************************************************
    ;; Create wgted area average of the base IO precip series (time)
    ;; Really, no need to area weight here .... area is very small.
    ;;************************************************
    P = pr_sm
    U = ua_sm

    PIO = wgt_areaave_Wrap(P(:, {latS_IO:latN_IO}, {lonL_IO:lonR_IO}), wyIO, 1., 0)
    PIO = dtrend(PIO, False)                              ; rmv overall trend
    PIO = wgt_runave_leftdim(PIO, wgt, 0)                 ; apply filter

    ;;************************************************
    ;; Create LAT average of the global Precip and U series (time, lon)
    ;; Really, no need to area weight here
    ;;************************************************

    ;; (time, lon)
    P_timeLon = dim_avg_Wrap(P(time|:, {lon|:}, {lat|lats:latn}))
    P_timeLon = dtrend_leftdim(P_timeLon, False)           ; rmv overall trend
    P_timeLon = wgt_runave_leftdim(P_timeLon, wgt, 0)     ; apply filter

    ;; (time, lon)
    U_timeLon = dim_avg_Wrap(U(time|:, {lon|:}, {lat|lats:latn}))
    U_timeLon = dtrend_leftdim(U_timeLon, False)           ; rmv overall trend
    U_timeLon = wgt_runave_leftdim(U_timeLon, wgt, 0)

    ;;************************************************
    ;; Create LON average of the global Precip and U series (time, lat)
    ;; Really, no need to area weight here
    ;;************************************************
    ;; (time, lat)
    P_timeLat = dim_avg_Wrap(P(time|:, lat|:, {lon|lonl:lonr}))
    P_timeLat = dtrend_leftdim(P_timeLat, False)            ; rmv overall trend
    P_timeLat = wgt_runave_leftdim(P_timeLat, wgt, 0)      ; apply filter

    ;; (time, lat)
    U_timeLat = dim_avg_Wrap(U(time|:, lat|:, {lon|lonl:lonr}))
    U_timeLat = dtrend_leftdim(U_timeLat, False)            ; rmv overall trend
    U_timeLat = wgt_runave_leftdim(U_timeLat, wgt, 0)

    optXcor = False
    mxlag = di@prop_filter_mxlag

    date_u = ymd
    rp_timelon = mjo_xcor_lag(PIO, P_timeLon, date_u, mxlag, season, optXcor)
    ru_timelon = mjo_xcor_lag(PIO, U_timeLon, date_u, mxlag, season, optXcor)

    rp_timelat = mjo_xcor_lag(PIO, P_timeLat, date_u, mxlag, season, optXcor)
    ru_timelat = mjo_xcor_lag(PIO, U_timeLat, date_u, mxlag, season, optXcor)

    return([/ rp_timelon, ru_timelon, rp_timelat, ru_timelat /])
end

;----------------------------------
undef("calc_pr_ua_mean")
function calc_pr_ua_mean(nml_data_vault_pr[1]: logical,\
                         nml_data_vault_ua[1]: logical,\
                                        di[1]: logical,\
                                idx_mod_pr[1]: numeric,\
                                idx_mod_ua[1]: numeric,\
                                  i_season[*]: numeric)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local mean_pr,   mean_ua,   pr,\
    pr_tmp,    pr_tmp2,   transforms,\
    ua,        ua_tmp,    ua_tmp2,\
    var0,      var1
begin
    transforms = True
    transforms@apply_these = (/"latrange", "lonrange"/)
    transforms@latrange = di@latrange_basic
    transforms@lonrange = di@lonrange_basic

    var0 = variables(0)
    var1 = variables(1)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    ua_tmp = extract_data_extended(idx_mod_ua,\
                                    var1,\
                                    nml_data_vault_ua,\
                                    -1,\
                                    0,\
                                    0)
    pr_tmp2 = apply_transforms(transforms, pr_tmp)
    ua_tmp2 = apply_transforms(transforms, ua_tmp)
    pr = clmMonTLL(pr_tmp2)
    ua = clmMonTLL(ua_tmp2)
    delete(ua_tmp)
    delete(pr_tmp)
    delete(ua_tmp2)
    delete(pr_tmp2)

    mean_ua = dim_avg_Wrap(ua(lat|:, lon|:, month|i_season(:)))
    mean_pr = dim_avg_Wrap(pr(lat|:, lon|:, month|i_season(:)))
    delete(ua)
    delete(pr)

    return([/ mean_pr, mean_ua /])
end

;-----------------------------------
undef("calc_pr_ua_vari")
function calc_pr_ua_vari(nml_data_vault_pr[1]: logical,\
                         nml_data_vault_ua[1]: logical,\
                                        di[1]: logical,\
                                idx_mod_pr[1]: numeric,\
                                idx_mod_ua[1]: numeric,\
                                   seaName[1]: string)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local iSeaS,         mms,          no_of_harmonics,\
    pr,            prBPF,        prClmDay,\
    prClmDay_sm,   pr_sm,        pr_tmp,\
    prVarBPF,      transforms,   ua,\
    uaBPF,         uaClmDay,     uaClmDay_sm,\
    ua_sm,         ua_tmp,       uaVarBPF,\
    var0,          var1,         wgt,\
    ymd,           ymdhms,       yyyyddd
begin
    ;; Compute lancoz weights
    ihp = di@ihp     ; bpf = >band pass filter
    nWgt = di@nWgt
    sigma = di@sigma   ; Lanczos sigma
    fca = di@fca     ; start freq
    fcb = di@fcb     ; end freq
    wgt = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

    transforms = True
    transforms@apply_these = (/"latrange", "lonrange", "fill_missing_w_poisson"/)
    transforms@latrange = di@latrange_seasonal
    transforms@lonrange = di@lonrange_seasonal
    transforms@is_cyclic = di@is_cyclic

    var0 = variables(0)
    var1 = variables(1)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; entire file

    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if

    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3
    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    pr_sm = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    ymdhms = cd_calendar(pr&time, 0)
    mms = floattoint(ymdhms(:, 1))        ; just months

    if (seaName.eq."summer") then
        iSeaS = ind(mms.ge.5 .and. mms.le.10)
    else if (seaName.eq."winter") then
        iSeaS = ind((mms.ge.1 .and. mms.le. 4) .or.\
                    (mms.ge.11.and. mms.le.12))
    end if
    end if

    prBPF = wgt_runave_Wrap(pr_sm(lat|:, lon|:, time|:), wgt, 0)
    prVarBPF = dim_variance_Wrap(prBPF(:, :, iSeaS))

    ua_tmp = extract_data_extended(idx_mod_ua,\
                                   var1,\
                                   nml_data_vault_ua,\
                                   -1,\
                                   0,\
                                   0)
    ua = apply_transforms(transforms, ua_tmp)
    uaClmDay = clmDayTLL_calendar(ua, yyyyddd)
    uaClmDay_sm = smthClmDayTLL_calendar(uaClmDay, no_of_harmonics)
    ua_sm = calcDayAnomTLL_calendar(ua, yyyyddd, uaClmDay_sm)

    uaBPF = wgt_runave_Wrap(ua_sm(lat|:, lon|:, time|:), wgt, 0)
    uaVarBPF = dim_variance_Wrap(uaBPF(:, :, iSeaS))

    delete(ua_tmp)
    delete(pr_tmp)
    delete(pr)
    delete(ua)
    delete(prClmDay)
    delete(prClmDay_sm)
    delete(uaClmDay)
    delete(uaClmDay_sm)
    delete(prBPF)
    delete(uaBPF)

    return([/ prVarBPF, uaVarBPF /])
end

undef("calc_olr_ua_spectra")
function calc_olr_ua_spectra(nml_data_vault_pr[1]: logical,\
                             nml_data_vault_ua[1]: logical,\
                                            di[1]: logical,\
                                    idx_mod_pr[1]: numeric,\
                                    idx_mod_ua[1]: numeric)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local aver,              interp_dir,    latd,\
    LATD1,             lond,          LOND1,\
    MLON5,             nFill,         NLAT5,\
    no_of_harmonics,   Opt,           pr,\
    prClmDay,          prClmDay_sm,   pr_sm,\
    pr_sm_tmp,         pr_tmp,        STC,\
    transforms,        ua,            uaClmDay,\
    uaClmDay_sm,       ua_sm,         ua_sm_tmp,\
    ua_tmp,            var0,          var1,\
    ymd,               yyyyddd
begin
    transforms = True
    transforms@apply_these = (/"latrange", "fill_missing_w_poisson"/)
    transforms@latrange = di@latrange_cross
    transforms@is_cyclic = di@is_cyclic

    var0 = variables(0)
    var1 = variables(1)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; entire file
    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if
    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3

    ;; replacing missing value 1e+20 with 0
    if (isatt(prClmDay, "_FillValue")) then
        nFill = num(ismissing(prClmDay))
        if (nFill.gt.0) then
            missv = prClmDay@_FillValue
            delete(prClmDay@_FillValue)
            prClmDay = where(prClmDay.eq.missv, 0, prClmDay)
            prClmDay@_FillValue = missv
            delete(missv)
        end if
    end if

    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    pr_sm = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    ua_tmp = extract_data_extended(idx_mod_ua,\
                                   var1,\
                                   nml_data_vault_ua,\
                                   -1,\
                                   0,\
                                   0)
    ua = apply_transforms(transforms, ua_tmp)
    uaClmDay = clmDayTLL_calendar(ua, yyyyddd)

    ;; Replacing missing value 1e+20 with 0
    if (isatt(uaClmDay, "_FillValue")) then
        nFill = num(ismissing(uaClmDay))
        if (nFill.gt.0) then
            missv = uaClmDay@_FillValue
            aver = avg(uaClmDay)
            delete(uaClmDay@_FillValue)
            uaClmDay = where(uaClmDay.eq.missv, aver, uaClmDay)
            uaClmDay@_FillValue = missv
            delete(missv)
            delete(aver)
        end if
    end if

    uaClmDay_sm = smthClmDayTLL_calendar(uaClmDay, no_of_harmonics)
    ua_sm = calcDayAnomTLL_calendar(ua, yyyyddd, uaClmDay_sm)

    ;; Interpolation the data to the 0.5 degree if
    ;; pr and ua from different datasets
    if (idx_mod_pr.ne. idx_mod_ua) then
        Opt = True
        Opt@ForceOverwrite = True
        Opt@PrintTimings = True
        Opt@Title = "0.5degree"

        NLAT5 = 361
        MLON5 = 720
        LATD1 = latGlobeF(NLAT5, "lat", "latitude", "degrees_north")
        LOND1 = lonGlobeF(MLON5, "lon", "longitude", "degrees_east")
        latS = di@latrange_cross(0)
        latN = di@latrange_cross(1)
        latd = LATD1({latS:latN})
        lond = LOND1(:)

        interp_dir = plot_dir + "/interpolation/"
        system("mkdir -p " + interp_dir)

        dstGridName_pr = interp_dir\
                         + models@name(idx_mod_pr)\
                         + "_" + "Scrip_05.nc"
        weight_pr = interp_dir\
                    + models@name(idx_mod_pr) + "_"\
                    + "weight_" + "conserve.nc"

        Opt@DstGridLat = latd
        Opt@DstGridLon = lond
        Opt@InterpMethod = "conserve"
        Opt@DstFileName = dstGridName_pr
        Opt@WgtFileName = weight_pr

        ;; Name of output source file
        Opt@SrcFileName = interp_dir\
                          + models@name(idx_mod_pr)\
                          + "_pr_src_SCRIP.nc"
        Opt@SrcTitle = "src Grid"

        pr_sm_tmp = pr_sm
        delete(pr_sm)
        pr_sm = ESMF_regrid(pr_sm_tmp, Opt)

        ; ua interpolation
        dstGridName_ua = interp_dir\
                         + models@name(idx_mod_ua) + "_"\
                         + "Scrip_05.nc"
        weight_ua = interp_dir\
                    + models@name(idx_mod_ua) + "_"\
                    + "weight_" + "conserve.nc"

        Opt@DstFileName = dstGridName_ua
        Opt@WgtFileName = weight_ua

        ;; Name of output source file
        Opt@SrcFileName = interp_dir\
                          + models@name(idx_mod_ua)\
                          + "_ua_src_SCRIP.nc"

        ua_sm_tmp = ua_sm
        delete(ua_sm)
        ua_sm = ESMF_regrid(ua_sm_tmp, Opt)

        delete(pr_sm_tmp)
        delete(ua_sm_tmp)
    end if

    ;;*****************************************************
    ;; Calculate the cross-spectra
    ;;*****************************************************
    segLen = di@segLen
    segOverLap = di@segOverLap
    STC = mjo_cross(pr_sm, ua_sm, segLen, segOverLap, False)

    delete(ua_tmp)
    delete(pr_tmp)
    delete(pr)
    delete(ua)
    delete(prClmDay)
    delete(prClmDay_sm)
    delete(uaClmDay)
    delete(uaClmDay_sm)
    delete(pr_sm)
    delete(ua_sm)

    return(STC)
end

;------------------------

undef("calc_daily_univari_eof")
function calc_daily_univari_eof(nml_data_vault_pr[1]: logical,\
                                               di[1]: logical,\
                                       idx_mod_pr[1]: numeric)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local eof,     eof_ts,     no_of_harmonics,\
    pr,      prClmDay,   prClmDay_sm,\
    pr_sm,   pr_tmp,     transforms,\
    var0,    wgt,        xbp,\
    ymd,     yyyyddd
begin
    transforms = True
    transforms@apply_these = (/"latrange"/)
    transforms@latrange = di@latrange_eof

    var0 = variables(0)
    neof = di@num_eof

    ;; Compute lancoz weights
    ihp = di@ihp      ; bpf = >band pass filter
    nWgt = di@nWgt
    sigma = di@sigma  ; Lanczos sigma
    fca = di@fca      ; start freq
    fcb = di@fcb      ; end freq
    wgt = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; Entire file

    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if

    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3
    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    pr_sm = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    ;; Do bandpass filter
    xbp = wgt_runave_Wrap(pr_sm(lat|:, lon|:, time|:), wgt, 0)

    ;; Remove means of band pass series: *not* necessary
    xbp = dim_rmvmean(xbp)

    ;; Calculate EOF
    eof = eofunc_Wrap(xbp, neof, False)
    eof_ts = eofunc_ts_Wrap(xbp, eof, False)

    delete(yyyyddd)
    delete(pr_tmp)
    delete(pr)
    delete(ymd)
    delete(prClmDay)
    delete(prClmDay_sm)
    delete(pr_sm)
    delete(xbp)

    return(eof)
end

;------------------------

undef("calc_daily_multi_vari_ceof")
function calc_daily_multi_vari_ceof(nml_data_vault_pr[1]: logical,\
                                    nml_data_vault_ua[1]: logical,\
                                                   di[1]: logical,\
                                           idx_mod_pr[1]: numeric,\
                                           idx_mod_ua[1]: numeric)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local ncdf,              ccr_12,            cdata,\
    ceof,               ceof_ts,           dimx,\
    eof_cdata,          eof_ts_cdata,      imax_200,\
    imax_850,           imin_200,          imin_850,\
    interp_dir,         latd,              LATD1,\
    lon,                lond,              LOND1,\
    lonmax_200,         lonmax_850,        mjo_ts_index,\
    mjo_ts_index_smt,   mlon,              MLON5,\
    mxlag,              netcdf_dir,        nGood,\
    NLAT5,              no_of_harmonics,   nStrong,\
    ntim,               nvar,              olr,\
    OLR,                Opt,               PC1,\
    PC2,                pcv_eof_olr,       pcv_eof_u200,\
    pcv_eof_u850,       pr,                prClmDay,\
    prClmDay_sm,        pr_sm,             pr_sm_tmp,\
    pr_tmp,             rlag_01,           rlag_10,\
    r_olr_u200,         r_olr_u850,        r_u850_u200,\
    time,               transforms,        u200,\
    U200,               u200ClmDay,        u200ClmDay_sm,\
    u200_sm,            u200_sm_tmp,       u200_tmp,\
    u200_tmp2,          u850,              U850,\
    u850ClmDay,         u850ClmDay_sm,     u850_sm,\
    u850_sm_tmp,        u850_tmp,          u850_tmp2,\
    ua_full,            var0,              var1,\
    var_olr,            var_u200,          var_u850,\
    wgt,                ymd,               yyyyddd,\
    zavg_var_olr,       zavg_var_u200,     zavg_var_u850
begin
    transforms = True
    transforms@apply_these = (/"latrange", "fill_missing_w_poisson"/)
    transforms@latrange = di@latrange_ceof
    transforms@is_cyclic = di@is_cyclic

    var0 = variables(0)
    var1 = variables(1)

    field_type0 = field_types(0)
    field_type1 = field_types(1)

    neof = di@num_ceof

    ;; Compute lancoz weights
    ihp = di@ihp      ; bpf = >band pass filter
    nWgt = di@nWgt
    sigma = di@sigma  ; Lanczos sigma
    fca = di@fca      ; start freq
    fcb = di@fcb      ; end freq
    wgt = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; entire file
    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if
    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3
    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    pr_sm = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    ua_full = extract_data_extended(idx_mod_ua,\
                                    var1,\
                                    nml_data_vault_ua,\
                                    -1,\
                                    0,\
                                    0)
    u200_tmp = ua_full(:, 0, :, :)
    u850_tmp = ua_full(:, 1, :, :)

    u200_tmp2 = apply_transforms(transforms, u200_tmp)
    u850_tmp2 = apply_transforms(transforms, u850_tmp)

    u200ClmDay = clmDayTLL_calendar(u200_tmp2, yyyyddd)
    u200ClmDay_sm = smthClmDayTLL_calendar(u200ClmDay, no_of_harmonics)
    u200_sm = calcDayAnomTLL_calendar(u200_tmp2, yyyyddd, u200ClmDay_sm)

    u850ClmDay = clmDayTLL_calendar(u850_tmp2, yyyyddd)
    u850ClmDay_sm = smthClmDayTLL_calendar(u850ClmDay, no_of_harmonics)
    u850_sm = calcDayAnomTLL_calendar(u850_tmp2, yyyyddd, u850ClmDay_sm)

    ;; interpolation the data to the 0.5 degree if
    ;; pr and ua from different datasets
    if (idx_mod_pr.ne. idx_mod_ua) then
        Opt = True
        Opt@ForceOverwrite = True
        Opt@PrintTimings = True
        Opt@Title = "0.5degree"

        NLAT5 = 361
        MLON5 = 720
        LATD1 = latGlobeF(NLAT5, "lat", "latitude", "degrees_north")
        LOND1 = lonGlobeF(MLON5, "lon", "longitude", "degrees_east")
        latS = di@latrange_cross(0)
        latN = di@latrange_cross(1)
        latd = LATD1({latS:latN})
        lond = LOND1(:)

        interp_dir = plot_dir + "/interpolation/"
        system("mkdir -p " + interp_dir)

        dstGridName_pr = interp_dir\
                         + models@name(idx_mod_pr) + "_"\
                         + "Scrip_05.nc"
        weight_pr = interp_dir\
                    + models@name(idx_mod_pr) + "_"\
                    + "weight_" + "conserve.nc"

        Opt@DstGridLat = latd
        Opt@DstGridLon = lond
        Opt@InterpMethod = "conserve"

        Opt@DstFileName = dstGridName_pr
        Opt@WgtFileName = weight_pr

        ;; Name of output source file
        Opt@SrcFileName = interp_dir\
                          + models@name(idx_mod_pr)\
                          + "_pr_src_SCRIP.nc"
        Opt@SrcTitle = "src Grid"
        pr_sm_tmp = pr_sm
        delete(pr_sm)

        pr_sm = ESMF_regrid(pr_sm_tmp, Opt)

        ;; u850 interpolation
        dstGridName_u850 = interp_dir\
                           + models@name(idx_mod_ua) + "_"\
                           + "Scrip_05.nc"
        weight_u850 = interp_dir\
                      + models@name(idx_mod_ua) + "_"\
                      + "weight_" + "conserve.nc"
        Opt@DstFileName = dstGridName_u850
        Opt@WgtFileName = weight_u850

        ;; Name of output source file
        Opt@SrcFileName = interp_dir\
                          + models@name(idx_mod_ua)\
                          + "_u850_src_SCRIP.nc"
        u850_sm_tmp = u850_sm
        delete(u850_sm)

        u850_sm = ESMF_regrid(u850_sm_tmp, Opt)

        ;; u200 interpolation
        dstGridName_u200 = interp_dir\
                           + models@name(idx_mod_ua) + "_"\
                           + "Scrip_05.nc"
        weight_u200 = interp_dir\
                      + models@name(idx_mod_ua) + "_"\
                      + "weight_" + "conserve.nc"
        Opt@DstFileName = dstGridName_u200
        Opt@WgtFileName = weight_u200

        ;; Name of output source file
        Opt@SrcFileName = interp_dir\
                          + models@name(idx_mod_ua)\
                          + "_u200_src_SCRIP.nc"
        u200_sm_tmp = u200_sm
        delete(u200_sm)

        u200_sm = ESMF_regrid(u200_sm_tmp, Opt)

        delete(pr_sm_tmp)
        delete(u850_sm_tmp)
        delete(u200_sm_tmp)
    end if

    OLR = dim_avg_n_Wrap(pr_sm, 1)
    U850 = dim_avg_n_Wrap(u850_sm, 1)
    U200 = dim_avg_n_Wrap(u200_sm, 1)

    delete(pr_tmp)
    delete(pr)
    delete(prClmDay)
    delete(prClmDay_sm)
    delete(pr_sm)
    delete(ua_full)
    delete(u200_tmp)
    delete(u850_tmp)
    delete(u200_tmp2)
    delete(u850_tmp2)
    delete(u200ClmDay)
    delete(u200ClmDay_sm)
    delete(u200_sm)
    delete(u850ClmDay)
    delete(u850ClmDay_sm)
    delete(u850_sm)

    dimx = dimsizes(OLR)
    ntim = dimx(0)
    mlon = dimx(1)
    time = OLR&time
    lon = OLR&lon

    ;;************************************************
    ;; Apply the band pass filter to the original anomalies
    ;;************************************************
    olr = wgt_runave_n_Wrap(OLR, wgt, 0, 0)  ; (time, lon)
    u850 = wgt_runave_n_Wrap(U850, wgt, 0, 0)
    u200 = wgt_runave_n_Wrap(U200, wgt, 0, 0)

    ;;************************************************
    ;; remove temporal means of band pass series: *not* necessary
    ;;************************************************
    olr = dim_rmvmean_n(olr, 0)              ; (time, lon)
    u850 = dim_rmvmean_n(u850, 0)
    u200 = dim_rmvmean_n(u200, 0)

    ;;************************************************
    ;; Compute the temporal variance at each lon
    ;;************************************************
    var_olr = dim_variance_n_Wrap(olr, 0)     ; (lon)
    var_u850 = dim_variance_n_Wrap(u850, 0)
    var_u200 = dim_variance_n_Wrap(u200, 0)

    ;;************************************************
    ;; Compute the zonal mean of the temporal variance
    ;;************************************************
    zavg_var_olr = dim_avg_n_Wrap(var_olr, 0)
    zavg_var_u850 = dim_avg_n_Wrap(var_u850, 0)
    zavg_var_u200 = dim_avg_n_Wrap(var_u200, 0)

    ;;************************************************
    ;; Normalize by sqrt(avg_var*)
    ;;************************************************
    olr = olr / sqrt(zavg_var_olr)          ; (time, lon)
    u850 = u850 / sqrt(zavg_var_u850)
    u200 = u200 / sqrt(zavg_var_u200)

    ;;************************************************
    ;; Combine the normalized data into one variable
    ;;************************************************
    cdata = new((/ 3 * mlon, ntim/), typeof(olr), getFillValue(olr))
    do ml = 0, mlon - 1
        cdata(ml, :) = (/olr(:, ml)/)
        cdata(ml + mlon, :) = (/u850(:, ml)/)
        cdata(ml + 2 * mlon, :) = (/u200(:, ml)/)
    end do

    ;;************************************************
    ;; Compute **combined** EOF; Sign of EOF is arbitrary
    ;;************************************************
    eof_cdata = eofunc(cdata, neof, False)      ; (neof, 3*mlon)
    eof_ts_cdata = eofunc_ts(cdata, eof_cdata, False)   ; (neof, 3*ntim)

    ;;************************************************
    ;; For clarity, explicitly extract each variable. Create time series
    ;;************************************************
    nvar = 3  ; "olr", "u850", "u200"
    ceof = new((/nvar, neof, mlon/), typeof(cdata), getFillValue(cdata))

    do n = 0, neof - 1
        ceof(0, n, :) = eof_cdata(n, 0:mlon - 1)         ; olr
        ceof(1, n, :) = eof_cdata(n, mlon:2 * mlon - 1)  ; u850
        ceof(2, n, :) = eof_cdata(n, 2 * mlon:)          ; u200
    end do

    ceof!0 = "var"
    ceof!1 = "eof"
    ceof!2 = "lon"
    ceof&lon = olr&lon
    ceof_ts = new((/nvar, neof, ntim/), typeof(cdata), getFillValue(cdata))

    ;; (0, neof, ntim)
    ceof_ts(0, :, :) = eofunc_ts_Wrap(olr(lon|:, time|:), ceof(0, :, :), False)

    ;; (1, neof, ntim)
    ceof_ts(1, :, :) = eofunc_ts_Wrap(u850(lon|:, time|:), ceof(1, :, :), False)

    ;; (2, neof, ntim)
    ceof_ts(2, :, :) = eofunc_ts_Wrap(u200(lon|:, time|:), ceof(2, :, :), False)

    ;;**********************************************t*
    ;; Add code contributed by Marcus N. Morgan, Florida Institute of
    ;; Technology; Feb 2015 Calculate % variance (pcv_)accounted ford
    ;; by OLR, U850 and U200
    ;;************************************************
    pcv_eof_olr = new(neof, typeof(ceof))
    pcv_eof_u850 = new(neof, typeof(ceof))
    pcv_eof_u200 = new(neof, typeof(ceof))

    do n = 0, neof - 1
        pcv_eof_olr(n) = avg((ceof(0, n, :) * sqrt(ceof@eval(n))) ^ 2) * 100
        pcv_eof_u850(n) = avg((ceof(1, n, :) * sqrt(ceof@eval(n))) ^ 2) * 100
        pcv_eof_u200(n) = avg((ceof(2, n, :) * sqrt(ceof@eval(n))) ^ 2) * 100
        print("pcv: neof = " + (n + 1) + ":  "\
                             + pcv_eof_olr(n) + "  "\
                             + pcv_eof_u850(n) + "  "\
                             + pcv_eof_u200(n))
    end do

    ;;**********************************************t*
    ;; Change sign? Ad hoc approach[ Eun-Pa Lim]
    ;;************************************************
                                        ;; index (subscript)
    imax_850 = maxind(ceof(1, 0, :))     ; max u850 eof 1
    imin_850 = minind(ceof(1, 0, :))     ; min

    imax_200 = maxind(ceof(1, 1, :))     ; max u200 eof 2
    imin_200 = minind(ceof(1, 1, :))     ; min

    lonmax_850 = ceof&lon(imax_850)      ; longitude of max value
    lonmax_200 = ceof&lon(imax_200)

    if (lonmax_850.lt.50 .or. lonmax_850.gt.120) then
        ceof(:, 0, :) = -ceof(:, 0, :)
        ceof_ts(:, 0, :) = -ceof_ts(:, 0, :)
        eof_cdata(0, :) = -eof_cdata(0, :)
        eof_ts_cdata(0, :) = -eof_ts_cdata(0, :)
    end if

    if (lonmax_200.lt.90 .or. lonmax_200.gt.150) then
        ceof(:, 1, :) = -ceof(:, 1, :)
        ceof_ts(:, 1, :) = -ceof_ts(:, 1, :)
        eof_cdata(1, :) = -eof_cdata(1, :)
        eof_ts_cdata(1, :) = -eof_ts_cdata(1, :)
    end if

    ;;************************************************
    ;; Compute cross correlation of each variable's EOF time series at zero-lag
    ;;************************************************
    r_olr_u850 = escorc(ceof_ts(0, :, :), ceof_ts(1, :, :))   ; (neof)
    r_olr_u200 = escorc(ceof_ts(0, :, :), ceof_ts(2, :, :))
    r_u850_u200 = escorc(ceof_ts(1, :, :), ceof_ts(2, :, :))

    ;;************************************************
    ;; Compute cross correlation of the multivariate EOF; EOF 1 vs EOF 2
    ;;************************************************
    mxlag = 25

    ;; (N, mxlag+1)
    rlag_01 = esccr(eof_ts_cdata(0, :), eof_ts_cdata(1, :), mxlag)

    ;; (N, mxlag + 1)
    rlag_10 = esccr(eof_ts_cdata(1, :), eof_ts_cdata(0, :), mxlag)
    ccr_12 = new((/2 * mxlag + 1/), float)

    ccr_12(mxlag:) = rlag_10(0:mxlag)
    ccr_12(0:mxlag) = rlag_01(:: - 1)       ; reverse order

    ;;************************************************
    ;; Normalize the multivariate EOF 1&2 component time series
    ;; Compute (PC1^2+PC2^2): values > 1 indicate "strong" periods
    ;;************************************************
    eof_ts_cdata(0, :) = eof_ts_cdata(0, :) / stddev(eof_ts_cdata(0, :))
    eof_ts_cdata(1, :) = eof_ts_cdata(1, :) / stddev(eof_ts_cdata(1, :))

    mjo_ts_index = eof_ts_cdata(0, :) ^ 2 + eof_ts_cdata(1, :) ^ 2
    mjo_ts_index_smt = runave(mjo_ts_index, 91, 0)  ; 91-day running mean

    nGood = num(.not.ismissing(mjo_ts_index))     ; # non-missing
    nStrong = num(mjo_ts_index .ge. 1.0)
    print("nGood = " + nGood\
                     + "   nStrong = " + nStrong\
                     + "   nOther = " + (nGood - nStrong))

    ;;************************************************
    ;; Write PC results to netCDF for use in another example.
    ;;************************************************
    mjo_ts_index!0 = "time"
    mjo_ts_index&time = time
    mjo_ts_index@long_name = "MJO PC INDEX"
    mjo_ts_index@info = "(PC1^2 + PC2^2)"

    PC1 = eof_ts_cdata(0, :)
    PC1!0 = "time"
    PC1&time = time
    PC1@long_name = "PC1"
    PC1@info = "PC1/stddev(PC1)"

    PC2 = eof_ts_cdata(1, :)
    PC2!0 = "time"
    PC2&time = time
    PC2@long_name = "PC2"
    PC2@info = "PC2/stddev(PC2)"

    ;;  diro = "./"
    netcdf_dir = plot_dir + "/netcdf/"
    system("mkdir -p " + netcdf_dir)

    aux_title_info = "MJO_PC_INDEX"
    output_netcdf = interface_get_figure_filename("mjo_life_cycle_plot",\
                                                  var0,\
                                                  field_type0,\
                                                  aux_title_info,\
                                                  idx_mod_pr)

    ;; Remove any pre - existing file
    system("/bin/rm -f " + netcdf_dir + output_netcdf + ".nc")
    ncdf = addfile(netcdf_dir + output_netcdf + ".nc", "c")

    ;; Make time an UNLIMITED dimension
    filedimdef(ncdf, "time", -1, True)      ; recommended  for most applications
                                            ;; output variables directly
    ncdf->MJO_INDEX = mjo_ts_index
    ncdf->PC1 = PC1
    ncdf->PC2 = PC2
    ;; ------------------------------------------

    return([/ ceof, pcv_eof_olr, pcv_eof_u850, pcv_eof_u200 /])
end

;------------------------
undef("calc_daily_olr_uv_life_cycle")
function calc_daily_olr_uv_life_cycle(nml_data_vault_pr[1]: logical,\
                                      nml_data_vault_ua[1]: logical,\
                                      nml_data_vault_va[1]: logical,\
                                                     di[1]: logical,\
                                             idx_mod_pr[1]: numeric,\
                                             idx_mod_ua[1]: numeric)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local Anom_sm_u,         Anom_sm_v,     Anom_sm_x,\
    dimx,              mlon,          nlat,\
    no_of_harmonics,   ntim,          pr,\
    prClmDay,          prClmDay_sm,   pr_tmp,\
    transforms,        u,             u850ClmDay,\
    u850ClmDay_sm,     u850_tmp,      u850_tmp2,\
    v,                 v850ClmDay,    v850ClmDay_sm,\
    v850_tmp,          v850_tmp2,     var0,\
    var1,              var2,          wgt,\
    x,                 ymd,           yyyyddd
begin
    transforms = True
    transforms@apply_these = (/"latrange", "fill_missing_w_poisson"/)
    transforms@latrange = di@latrange_life
    transforms@is_cyclic = di@is_cyclic

    var0 = variables(0)
    var1 = variables(1)
    var2 = variables(2)

    ;; Compute lancoz weights
    ihp = di@ihp      ; bpf = >band pass filter
    nWgt = di@nWgt
    sigma = di@sigma  ; Lanczos sigma
    fca = di@fca      ; start freq
    fcb = di@fcb      ; end freq
    wgt = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

    pr_tmp = extract_data_extended(idx_mod_pr,\
                                   var0,\
                                   nml_data_vault_pr,\
                                   -1,\
                                   0,\
                                   0)
    pr = apply_transforms(transforms, pr_tmp)
    ymd = cd_calendar(pr&time, -2)   ; entire file
    if (isatt(ymd, "calendar")) then
        if (ymd@calendar .eq. "proleptic_gregorian") then
            ymd@calendar = "gregorian"
        end if
    end if
    yyyyddd = yyyymmdd_to_yyyyddd(ymd)
    prClmDay = clmDayTLL_calendar(pr, yyyyddd)
    no_of_harmonics = 3
    prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
    Anom_sm_x = calcDayAnomTLL_calendar(pr, yyyyddd, prClmDay_sm)

    u850_tmp = extract_data_extended(idx_mod_ua,\
                                     var1,\
                                     nml_data_vault_ua,\
                                     -1,\
                                     0,\
                                     0)
    v850_tmp = extract_data_extended(idx_mod_ua,\
                                     var2,\
                                     nml_data_vault_va,\
                                     -1,\
                                     0,\
                                     0)

    u850_tmp2 = apply_transforms(transforms, u850_tmp)
    v850_tmp2 = apply_transforms(transforms, v850_tmp)

    u850ClmDay = clmDayTLL_calendar(u850_tmp2, yyyyddd)
    u850ClmDay_sm = smthClmDayTLL_calendar(u850ClmDay, no_of_harmonics)
    Anom_sm_u = calcDayAnomTLL_calendar(u850_tmp2, yyyyddd, u850ClmDay_sm)

    v850ClmDay = clmDayTLL_calendar(v850_tmp2, yyyyddd)
    v850ClmDay_sm = smthClmDayTLL_calendar(v850ClmDay, no_of_harmonics)
    Anom_sm_v = calcDayAnomTLL_calendar(v850_tmp2, yyyyddd, v850ClmDay_sm)

    delete(pr_tmp)
    delete(pr)
    delete(prClmDay)
    delete(prClmDay_sm)
    delete(v850_tmp)
    delete(u850_tmp)
    delete(v850_tmp2)
    delete(u850_tmp2)
    delete(v850ClmDay)
    delete(v850ClmDay_sm)
    delete(u850ClmDay)
    delete(u850ClmDay_sm)

    ;;---------------------------------------

    dimx = dimsizes(Anom_sm_x)
    ntim = dimx(0)
    nlat = dimx(1)
    mlon = dimx(2)

    ;;************************************************
    ;; Apply the band pass filter to the original anomalies
    ;;************************************************

    x = wgt_runave_n(Anom_sm_x, wgt, 0, 0)      ; better
    u = wgt_runave_n(Anom_sm_u, wgt, 0, 0)
    v = wgt_runave_n(Anom_sm_v, wgt, 0, 0)

    copy_VarCoords(Anom_sm_x, x)
    copy_VarCoords(Anom_sm_u, u)
    copy_VarCoords(Anom_sm_v, v)

    delete(Anom_sm_x)
    delete(Anom_sm_u)
    delete(Anom_sm_v)

    return([/ x, u, v /])
end
