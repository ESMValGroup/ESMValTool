;This file calculate the equivalent latitude from a potential vorticity (PV) field on a prescribed isentropic level.


;***********************************************************
load "./variable_defs/variable_defs_aux/function_eq_lat.ncl"

load "./diag_scripts/uniform_coordinates.ncl"

;************************************************************
undef("CalcEqLat")

function CalcEqLat(data,field)

local PV,time,lat,lon,nlat,nlov,nla,nlo,num_eq_lat,array_lat_lon,countour_S,countour_N,\
Equivalent_latitude_N,Equivalent_latitude_S,nl,nln,nlon,nlat,diff

begin
   PV=data  
  ; printVarSummary(PV)
   time=PV&time
   lat=PV&lat
   lon=PV&lon
   PV&lat=lat
   PV&lon=lon
   PV&time=time

   nlat=dimsizes(lat)
   nlon=dimsizes(lon)
        
   uniform_coordinates(PV(time|:,lat|:,lon|:),field)
   ;printVarSummary(PV)

   nla=73
   nlo=144
   num_eq_lat=100
   array_lat_lon=new((/dimsizes(time),nlat,nlon/),"double")
   countour_S=new((/dimsizes(time),num_eq_lat/),"double")
   countour_N=new((/dimsizes(time),num_eq_lat/),"double")

   do ntime=0,dimsizes(time)-1       
      countour_S(ntime,:)=  equivalent_latitude(PV(lat|:,lon|:,time|ntime),nla,nlo,num_eq_lat,"SH")
      countour_N(ntime,:)=  equivalent_latitude(PV(lat|:,lon|:,time|ntime),nla,nlo,num_eq_lat,"NH")

      Equivalent_latitude_S=countour_S&Equivalent_latitude
      Equivalent_latitude_N=countour_N&Equivalent_latitude
      do nl=0,nlat-1
         if(lat(nl).le.0.) then
            do nln=0,nlon-1               
               diff=countour_S          
               do k=0,num_eq_lat-1       

                  if(.not.(ismissing(countour_S(ntime,k)))) then                    
                     diff(ntime,k)=countour_S(ntime,k)-PV(lat|nl,lon|nln,time|ntime)                       
                  else
                     diff(ntime,k)=countour_S@_FillValue                     
                  end if
               end do

               if(.not.all(ismissing(diff))) then
                  array_lat_lon(ntime,nl,nln)=Equivalent_latitude_S(ind(min(diff)))
               else
                  array_lat_lon=countour_S@_FillValue
               end if                               
               delete(diff)
            end do            
         end if
         if(lat(nl).gt.0.) then
            do nln=0,nlon-1               
               diff=countour_N          
               do k=0,num_eq_lat-1       	            
                  if(.not.(ismissing(countour_N(ntime,k)))) then                    
                     diff(ntime,k)=countour_N(ntime,k)-PV(lat|nl,lon|nln,time|ntime)                       
                  else
                     diff(ntime,k)=countour_N@_FillValue
                  end if
               end do

            if(.not.all(ismissing(diff))) then
               array_lat_lon(ntime,nl,nln)=Equivalent_latitude_N(ind(min(diff)))
            else
               array_lat_lon=countour_N@_FillValue
            end if   
            delete(diff)
           end do
        end if
      end do    
   end do


   array_lat_lon!0="time"
   array_lat_lon!1="lat"
   array_lat_lon!2="lon"

   array_lat_lon&time=time
   array_lat_lon&lat=lat
   array_lat_lon&lon=lon
   return(array_lat_lon)

end
       
