;
; Requires: ta:*3*
;
info = True
info@derived = True   
;info@derived = False 
; old: ta:T3*
info@long_name="Tropopause Air Altitude"
info@units="km"
info@scp_cnlvl=(/9,10,11,12,13,14,15,16,17,18,19,20/) ;
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_cnlvl=(/9,10,11,12,13,14,15,16,17,18,19,20/)
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@acp_cnlvl=(/9,10,11,12,13,14,15,16,17,18,19,20/)
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)

;surfcontrend: surface contour trend plots
info@sct_cnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*1.e-3
info@sct_dcnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*1.e-3

info@refModel = "CMAM"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

;tsline 
info@ts_latrange = (/-25,25/)
info@ts_lonrange = (/0,360/)
info@ts_anom = True

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
local result, i, pr, T, pin, tin, ptp, ztp
begin
  data_new=True
  result = read_data(index,"ta","*3*")
;  missing=-999.0
    T = extract_data(index,result,-1,0,0)
    printVarSummary(T)
;    dims=dimsizes(T)
;    ptp=new((/dims(0),dims(2),dims(3)/),"float")

    pr=T&plev
    
    ptp = trop_wmo(pr(::-1), T({time | :},{lat | :},{lon | :},{plev | ::-1}), 0, False) 
    

;   for outputing altitude (km) instead
    ztp = -7.5*log(ptp/1013.)
    delete(ptp)

    ztp!0=T!0
    ztp!1=T!2
    ztp!2=T!3
    ztp&time=T&time
    ztp&lat=T&lat
    ztp&lon=T&lon
    
    ztp@_FillValue = T@_FillValue
    ztp@cell_methods = "time:mean (of calendar months)"
    ztp@long_name = info@long_name
    ztp@units = info@units
    ztp@standard_name = info@long_name

    add_data(index,data_new,ztp)
    delete(T)
    delete(ztp)
    delete(pr)


  c=stringtochar(field_number)
  c(1)=integertochar(50) ; "2"
  c(3)=integertochar(115) ; "s"
  field_number=chartostring(c)
  print(field_number)
  delete(c)
  return(data_new)
end

