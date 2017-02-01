;
; Requires: va:T3M
; Monsoon Hadley Index compared to precipitation
info = True
info@derived = True
info@long_name="MHI"
info@units="Standard deviation"
info@vcp_cnlvl=(/0,1,2,3,4,5,6,7,8,9/)*10.
info@vcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@vcp_ptop=6e-2
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

info@scp_cnlvl = (/1.0625,2.125,3.1875,4.25,5.3125,6.375,7.4375,8.5,9.5625,10.625,11.6875,12.75,13.8125,14.875,15.9375/)*1e-5
info@pcp_cnlvl = (/1.0625,2.125,3.1875,4.25,5.3125,6.375,7.4375,8.5,9.5625,10.625,11.6875,12.75,13.8125,14.875,15.9375/)*1e-5

;info@acp_cnlvl = (/1.06,2.13,3.19,4.3,5.3,6.38,7.44,8.50,9.56,10.63,11.69,12.75,13.81,14.88,15.94/)*1e-5
info@acp_cnlvl = (/1.06,2.13,3.19,4.3,5.3,6.38,7.44,8.50,9.56,10.63,11.69,12.75,13.81,14.88,15.94/)
info@acp_dcnlvl = (/-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30/)*1e-6

info@scp_dcnlvl = (/-20.0,-17.5,-15.0,-12.5,-10.0,-7.5,-5.0,-2.5,0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5/)*1e-6

info@cf_latrange = (/9,31/)
info@cf_lonrange = (/69,111/)
info@ts_anom = False
;info@masks=(/"landsea.nc"/)

info@title = "Monsoon Hadley Index"

info@ac_latrange = (/-20,40/)
;info@ac_lonrange = (/77,100/)
info@refModel = (/"ERAINT"/)
info@season="JJAS"

load "./plot_type/calc_seasonal.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;for lsm tmp
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "ncl_code/apply_mask.ncl"


undef("calculate")
function calculate(index [1] : integer, 
                   variable [1] : string, 
                   field_number [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'temp_ncl_code/filesinfo'-file
;;    variable - logical with releveant variable as string attribute
;;    field_number  - string with field number classification
;; Return value:
;;    data_new - logical
;; Description:
;;    Computes the Monsoon Hadley Index for the designated models. See
;;    Goswami B.N., Krishnamurthy V., Annamalai H., 1999, A broad-scale circu-
;;    lation index for the interannual variability of the Indian summer monsoon,
;;    Quarterly Journal of the Royal Meteorological Society, vol. 125, issue 554,
;;    pp. 611-633
local result, i, T, lsdata, a, lsm, va850, va200, climva850_time, climva200_time
begin
    data_new=True
    result_va = read_data(index, "va", "T3M")
    data_va = extract_data(index, result_va, -1, 0, 0)

    ;; compute the climatology plev 850
    va850=data_va(:, {850}, :, :)
    wgt_lat=conform_dims((/dimsizes(va850&lat)/), 1.0, -1)
    climva850=clim_season_average(va850, wgt_lat, info@season)

    ;; compute the climatology plev 200
    va200=data_va(:, {200}, :, :)
    wgt_lat=conform_dims((/dimsizes(va200&lat)/), 1.0, -1)
    climva200=clim_season_average(va200, wgt_lat, info@season)

    ;; compute monsoon-hadley-index (mhi) 
    climva850_time=conform( va850, climva850, (/1,2/) )
    climva200_time=conform( va200, climva200, (/1,2/) )
    mhi = ( va850 - climva850_time ) - ( va200 - climva200_time ) 

    copy_VarMeta(va850, mhi)
    mhi@standard_name="mhi"
    mhi@long_name="Monsoon Hadely Index"
    mhi@units="Standard deviation"

    ;; Rewrite field specification to T2Ms (2D timeseries on the surface)
    c=stringtochar(field_number)
    c(1)=integertochar(50) ; "2"
    c(3)=integertochar(115) ; "s"
    field_number=chartostring(c)
    delete(c)

    add_data(index, data_new, mhi)
    delete(climva200)
    delete(climva850)
    delete(climva200_time)
    delete(climva850_time)
    delete(mhi)

    return(data_new)
end


