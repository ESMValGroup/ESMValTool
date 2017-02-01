;
; Requires: ta:*3*
;
info = True
info@derived = True   
;OLD Requires: ta:T3*
;info@derived = False 
info@long_name="Tropopause Air Pressure"
info@units="hPa"

info@scp_cnlvl=(/50,150,200,250,300,350,400,450,500,550,600/) ; ptp
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_cnlvl=(/50,150,200,250,300,350,400,450,500,550,600/) ; ptp
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl=(/50,150,200,250,300,350,400,450,500,550,600/) ; ptp
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@refModel = "ERA40"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

info@ts_latrange = (/-20,20/)
info@ts_lonrange = (/0,360/)
info@ts_anom = False
info@ts_multimean = True

info@ml_latrange = (/-20,20/)
info@ml_lonrange = (/0,360/)
info@ml_anom = False
info@ml_stdv_plot = "ERA40"
info@ml_multimean = True

;surfcontrend: surface contour trend plots
info@sct_cnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e-2
info@sct_dcnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e-2

info@fig07_lat_min = (/-10./)
info@fig07_lat_max = (/10./)
info@fig07_lev_sel = (/80./)


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
local tc, pr, T, i, pin, tin, ptp
begin
  data_new=True
  result = read_data(index,"ta","*3*")
  missing=-999.0
  T = extract_data(index,result,-1,0,0)


;check ordering....needs to be top to bottom (lowest pressure first)

  pr=T&plev
  if (pr(0).eq.max(pr)) then
    pin =  pr(::-1)
    tin =  T({time | :},{lat | :},{lon | :},{plev | ::-1})
  else
    pin = pr
    tin = T({time | :},{lat | :},{lon | :},{plev | :})
  end if	 

  if(max(pin).gt.2000.) then
    pin=pin/100.
  end if

  ptp = trop_wmo(pin, tin, 0, False)

; original that did not check for ordering
;    ptp = trop_wmo(pr(::-1), T({time | :},{lat | :},{lon | :},{plev | ::-1}), 0, False) 
 
  delete(pin)
  delete(tin)

;   for outputing altitude (km) instead
;    ztp = -7.500*log(ptp/1013.)

  ptp!0=T!0
  ptp!1=T!2
  ptp!2=T!3
  ptp&time=T&time
  ptp&lat=T&lat
  ptp&lon=T&lon
  ptp@long_name = info@long_name
  ptp@units = info@units
    
;    ptp@_FillValue = T@_FillValue
    ptp@_FillValue = -999.

;    printVarSummary(ptp)
;    print(max(ptp))
;    print(min(ptp))

  add_data(index,data_new,ptp)
  delete(T)
  delete(ptp)
  delete(pr)

  c=stringtochar(field_number)
  c(1)=integertochar(50) ; "2"
  c(3)=integertochar(115) ; "s"
  field_number=chartostring(c)
  print(field_number)
  delete(c)
  return(data_new)
end

