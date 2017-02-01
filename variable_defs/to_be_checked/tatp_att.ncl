;
; Requires: ta:*3*
;
info = True
info@derived = True   
;OLD Requires: ta:T3*
;info@derived = False 
info@long_name="Tropopause Air Temperature"
info@units="K"
info@scp_cnlvl=(/180,185,190,195,200,205,210,215,220,225/) ; ptp
info@pcp_cnlvl=(/180,185,190,195,200,205,210,215,220,225/) 
info@acp_cnlvl=(/180,185,190,195,200,205,210,215,220,225/)
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl = (/-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,12/)

;surfcontrend: surface contour trend plots
info@sct_cnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*1.e-2
info@sct_dcnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*1.e-2

info@refModel = "WACCM"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

info@ts_latrange = (/-25,25/)
info@ts_lonrange = (/0,360/)
info@ts_anom = False

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
local result, i, T, pr, Twork, ptp, ptropWork, Pwork, twork, ttrop
begin
  data_new=True
  result = read_data(index,"ta","*3*")
;  missing=-999.0
    T = extract_data(index,result,-1,0,0)
    printVarSummary(T)
;    dims=dimsizes(T)
;    ptp=new((/dims(0),dims(2),dims(3)/),"float")

    pr=T&plev

    Twork=T({time | :},{lat | :},{lon | :},{plev | ::-1})

    ptp = trop_wmo(pr(::-1), Twork, 0, False) 

;************************************************
; Interpolate Temperature to the 'ptrop' level
;************************************************

    ptropWork = conform(Twork, ptp, (/0,1,2/))
;    printVarSummary(ptropWork)
;    printVarSummary(Twork)
;    printVarSummary(ptp)

    Pwork=conform(Twork,pr(::-1),(/3/))

;    printVarSummary(Pwork)



    twork     = int2p(Pwork,Twork, ptropWork, 0)

    printVarSummary(twork)

;    ttrop     = Twork(:,:,:,0)             ; trick for meta data
    ttrop     = (/ twork(:,:,:,0) /)
    ttrop@long_name  = "temperature at tropopause"
    ttrop@units      = T@units
    ttrop!0 = T!0
    ttrop!1 = T!2
    ttrop!2 = T!3
    ttrop&time=T&time
    ttrop&lat=T&lat
    ttrop&lon=T&lon

    printVarSummary(ttrop)

;   for outputing altitude (km) instead
;    ztp = -7.500*log(ptp/1013.)

    add_data(index,data_new,ttrop)
    delete(T)
    delete(ptp)
    delete(pr)
    delete(ttrop)
    delete(Twork)
    delete(Pwork)
    delete(ptropWork)
    delete(twork)

  c=stringtochar(field_number)
  c(1)=integertochar(50) ; "2"
  c(3)=integertochar(115) ; "s"
  field_number=chartostring(c)
  print(field_number)
  delete(c)
  return(data_new)
end

