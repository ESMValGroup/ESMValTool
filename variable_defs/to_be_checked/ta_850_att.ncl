;
; Requires: ta:*3*
;
info = True
info@derived = True
info@long_name="850 mb Temp"
info@units="K"
info@cell_method = "time: mean"
info@sdp_cnlvl = (/-30,-25,-20,-15,-10,-5,-2,0,2,5,10,15,20,25,30/)
info@scp_cnlvl = (/230,235,240,245,250,255,260,265,270,275,280,285,290,295,300/)
info@scp_dcnlvl = (/-8,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,8/)
info@acp_cnlvl = (/230,235,240,245,250,255,260,265,270,275,280,285,290,295,300/)
info@acp_dcnlvl = (/-8,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,8/)
info@rgb_file = "amwg.rgb"
info@pcp_cnlvl = (/210,220,230,240,250,255,260,265,270,275,280,285,290,295,300/)
info@pcp_dcnlvl = (/-8,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,8/)
info@refModel = "WACCM"
;info@table_one_type = "DIFF"
info@seasons = (/"ANN","DJF","MAM","JJA","SON"/)
info@regions = (/"GLBL"/)

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
local tmp, dum, i
begin
  data_new = True
  tmp = read_data(index,"ta","*3*")
    dum = extract_data(index,tmp,1,850.,850.)   ;extract_data(models,tmp,i,"plev",850.,850.)

    dum@long_name=info@long_name
    dum@units=info@units
    add_data(index,data_new,dum)

    delete(dum)
  delete(tmp)

  return(data_new)
end
