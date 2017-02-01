;
; Requires: none
;
info = True
info@derived = True
info@long_name="Surface pressure"
info@units= "millibars"
info@pcp_cnlvl = (/650,700,750,800,820,840,860,880,900,920,940,960,980,1000,1020/)
info@pcp_dcnlvl = (/-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15/)
info@scp_cnlvl = (/600,630,660,690,720,750,780,810,840,870,900,930,960,990,1020/)
info@scp_dcnlvl = (/-30,-25,-20,-15,-10,-5,-2,0,2,5,10,15,20,25,30/)
info@sdp_cnlvl = (/-30,-25,-20,-15,-10,-5,-2,0,2,5,10,15,20,25,30/)
info@acp_cnlvl = (/600,630,660,690,720,750,780,810,840,870,900,930,960,990,1020/)
info@acp_dcnlvl = (/-30,-25,-20,-15,-10,-5,-2,0,2,5,10,15,20,25,30/)
info@rgb_file = "amwg.rgb"
info@refModel = "UMETRAC"
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
local ps
; convert PS from Pa to mb and return
begin
  result = read_data(index,variable,field_number)
  ps = extract_data(index,result,-1,0,0)
  ps = (/ps/100./)       ; convert to mb
  ps@units = info@units
  if (typeof(ps).eq."double") then
    ps_flt = dble2flt(ps)
  else
    ps_flt = ps
  end if
  modify_data(index,result,ps_flt)
  delete(ps)
  delete(ps_flt)

  return(result)
end
