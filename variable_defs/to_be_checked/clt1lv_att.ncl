;
;  Requires: clt:*3*
;
info = True
info@derived = True
info@long_name="Avg Cloud Fraction"
info@units="fraction"
info@scp_cnlvl=(/-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80/)*0.01 
info@pcp_cnlvl=(/-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80/)*0.01
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*0.01
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*0.01
info@sdp_cnlvl = (/-30,-25,-20,-15,-10,-5,-2,0,2,5,10,15,20,25,30/)*0.01
info@acp_cnlvl=(/-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80/)*0.01
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*0.001
info@sct_cnlvl=(/-10,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,10/)*0.001
info@sct_dcnlvl=(/-10,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,10/)*0.001
info@seasons = (/"ANN","DJF","MAM","JJA","SON"/)
info@regions = (/"GLBL"/)
;info@refModel = "WACCM"

info@vcp_ptop=6e-2
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

info@lvrange=(/130.,90./)

info@ts_latrange = (/-15,-15/)
info@ts_lonrange = (/60,180/)
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
local tmp, dum, data, i
begin
  data_new = True
  tmp = read_data(index,"clt","*3*")
  dum = extract_data(index,tmp,1,info@lvrange(0),info@lvrange(1))
  data = dim_avg_n_Wrap(dum,1)
  data@long_name=info@long_name
  data@units = info@units
  add_data(index,data_new,data)

  delete(dum)
  delete(data)
  delete(tmp)

  return(data_new)
end
