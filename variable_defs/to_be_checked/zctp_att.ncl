;
; Requires: ta:*3*
;
info = True
info@derived = True
info@long_name="Cold Point Tropopause Air Altitude"
info@units="km"
info@scp_cnlvl=(/9,10,11,12,13,14,15,16,17,18,19,20/) ;
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_cnlvl=(/9,10,11,12,13,14,15,16,17,18,19,20/)
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@acp_cnlvl=(/9,10,11,12,13,14,15,16,17,18,19,20/)
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@refModel = "WACCM"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"


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
local result, missing, i, dims, zctp, z, itime, ilat, ilon, T, tc, \
      xx, ng
begin
  data_new=True
  result = read_data(index,"ta","*3*")
  missing=-999.0
    print(" ")
    T = extract_data(index,result,-1,0,0)
;caf    printVarSummary(T)
    dims=dimsizes(T)
    zctp=new((/dims(0),dims(2),dims(3)/),"float")

    z = -7500*log(T&plev/1013)
    do itime = 0,dimsizes(T&time)-1
      do ilat  = 0,dimsizes(T&lat)-1
        do ilon  = 0,dimsizes(T&lon)-1
         tc=T(itime,:,ilat,ilon)
    
          pr=tc&plev
          if (typeof(tc).eq."double") then
            tc_flt = dble2flt(tc)
            delete(tc)
            tc=tc_flt
            delete(tc_flt)
          end if
;-----------------------------------------------------------------------------
;find cold point trop
;-----------------------------------------------------------------------------

          xx=ind(pr.lt.500..and.pr.gt.60.)
          ng=dimsizes(xx)
          if (ng .gt. 0) then 
              cpt=min(tc(xx))
              pp=ind(tc(xx) .eq. cpt)
;              printVarSummary(z(xx(pp(0))))
;              print(z(xx(pp(0))))
;              printVarSummary(zctp(itime,ilat,ilon))
;              printVarSummary(zctp)
              zctp(itime,ilat,ilon)=doubletofloat(z(xx(pp(0))) / 1000.)


              delete(pp)
              delete(cpt)
          else 
              zctp(itime,ilat,ilon)=missing
          end if

          delete(xx)
       end do  ;lon
     end do    ;lat
   end do      ;time


    zctp!0=T!0
    zctp!1=T!2
    zctp!2=T!3
    zctp&time=T&time
    zctp&lat=T&lat
    zctp&lon=T&lon
    zctp@long_name = info@long_name
    zctp@units = info@units
    add_data(index,data_new,zctp)
    delete(T)
    delete(zctp)
    delete(tc)
    delete(z)
    delete(pr)



  c=stringtochar(field_number)
  c(1)=integertochar(50) ; "2"
  c(3)=integertochar(115) ; "s"
  field_number=chartostring(c)
  print(field_number)
  delete(c)
  return(data_new)
end

