;
; Requires: ta:*3*
;
;OLD: Requires: ta:T3*
info = True
info@derived = True   
;info@derived = False
info@long_name="Cold Point Tropopause Air Temperature"
info@units="K"
info@scp_cnlvl=(/180,185,190,195,200,205,210,215,220,225/) ; ptp
info@pcp_cnlvl=(/180,185,190,195,200,205,210,215,220,225/) 
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl = (/-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,12/)
info@refModel = "ERA40"
info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

;surfcontrend: surface contour trend plots
info@sct_cnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*1.e-2
info@sct_dcnlvl= (/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*1.e-2

info@ts_latrange = (/-20,20/)
info@ts_lonrange = (/0,360/)
;info@ts_ymax = 220
;info@ts_ymin = 180
;info@ts_xmax_year = 2010
;info@ts_xmin_year = 1950
info@ts_anom = False
info@ts_multimean = True

info@ml_latrange = (/-20,20/)
info@ml_lonrange = (/0,360/)
info@ml_anom = False
info@ml_stdv_plot = "ERA40"
info@ml_multimean = True

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
local result, missing, i, T, dims, tactp, z, itime, ilat, ilon, tc, \
      xx, ng
begin
  data_new=True
  result = read_data(index,"ta","*3*")
  missing=-999.0
    print(" ")
    T = extract_data(index,result,-1,0,0)
    if (max(T&plev).gt.2000.) then
      T&plev = T&plev / 100.
    end if
    dims=dimsizes(T)
    tactp=new((/dims(0),dims(2),dims(3)/),"float")

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
              tactp(itime,ilat,ilon)=min(tc(xx))
          else 
              tactp(itime,ilat,ilon)=missing
          end if

       end do  ;lon
     end do    ;lat
   end do      ;time


    tactp!0=T!0
    tactp!1=T!2
    tactp!2=T!3
    tactp&time=T&time
    tactp&lat=T&lat
    tactp&lon=T&lon
    tactp@long_name=info@long_name
    tactp@units=info@units

    add_data(index,data_new,tactp)
    delete(T)
    delete(tactp)
    delete(tc)
    delete(z)
    delete(pr)
    delete(xx)


  c=stringtochar(field_number)
  c(1)=integertochar(50) ; "2"
  c(3)=integertochar(115) ; "s"
  field_number=chartostring(c)
  print(field_number)
  delete(c)
  return(data_new)
end

