;This file calculate the equivalent latitude from a potential vorticity (PV) field on a prescribed isentropic level.

;author Irene Cionni

undef("equivalent_latitude")

function equivalent_latitude(PV:numeric,jlat:integer,ilon:integer,num_eq_lat:integer,he:string)

local extract_PV,lat_new_S,lat_new_N,nlat,nlon,pi_p,degreetoradians,lat,lon,\
countour_S,countour_N,max_PV_S,min_PV_S,max_PV_N,min_PV_N,AJ_S,AT_S,AJ_N,AT_N,\
Num_N,Num_S,Tot_N,Tot_S,AT_S,AT_N,e_lat_S,dsizes_a_S,dsizes_a_N,indices_S,\
indices_N,ind_res_S,ind_res_N,dim1_S,dim1_N,A1_S,A1_N,lat_p_S,lat_p_N,lat_m_S,\
lat_m_N,lat_new_S,lat_new_N,AJ_S,AJ_N,AT_S,AT_N,e_lat_S,e_lat_N

begin
   
   PV=PV
   lat=PV&lat
   lon=PV&lon

   nlat=dimsizes(lat)
   nlon=dimsizes(lon)
   pi_p=new(1,"double")
   pi_p=-1.0       
   pi=acos(pi_p)
   degreetoradians= 2*acos(pi_p)/360

   ;The function g2fsh read latitude from -90 to 90 and does'nt want ciclic point in longitude values
   ;extract_PV =g2fsh_Wrap(PV(z0_level|:,time0|:,{latitude0|-90.0:90.0},{longitude0|0:359.}), (/jlat,ilon/) )

   extract_PV =g2fsh_Wrap(PV(lat|:,lon|:), (/jlat,ilon/) )
   NG=num(ismissing(PV))
   ;printVarSummary(extract_PV)

   extr_PV_S=extract_PV(lat|0:(floattointeger(jlat/2)),lon|:)
   extr_PV_N=extract_PV(lat|(floattointeger(jlat/2)):(jlat-1),lon|:)

   lat_new_S=extr_PV_S&lat
   lat_new_N=extr_PV_N&lat


   ;print("****************selection of maximum and minimum PV value in the hemisphere***********")
   ;**********************************************************************************************
   countour_S=new((/num_eq_lat/),"double")
   countour_N=new((/num_eq_lat/),"double")

   max_PV_S=max(extr_PV_S)
   min_PV_S=min(extr_PV_S)
   max_PV_N=max(extr_PV_N)
   min_PV_N=min(extr_PV_N)

   ;***********************************************************************************************
   ;print("*********************selection of 100 level between max and min************************")
   ;***********************************************************************************************
   if((any(ismissing(max_PV_N))).or.(any(ismissing(min_PV_N)))) then
      print("Your data contains some missing values. Beware.")
      countour_N= extract_PV@_FillValue        
   else
      countour_N=fspan(min_PV_N,max_PV_N,num_eq_lat)
   end if
   if(any(ismissing(max_PV_S)).or.any(ismissing(min_PV_S))) then
      print("Your data contains some missing values. Beware.")
      countour_S= extract_PV@_FillValue        
   else
      countour_S=fspan(min_PV_S,max_PV_S,num_eq_lat)
   end if
   ;**********************************************************************************************
   ;print("*********************geographic area enclosed by each of 100 PV isolines************")
   ;**********************************************************************************************
   AJ_S=new((/num_eq_lat/),double)
   AT_S=new((/num_eq_lat/),double)
   AJ_N=new((/num_eq_lat/),double)
   AT_N=new((/num_eq_lat/),double)
   Num_N=new((/num_eq_lat/),integer)
   Num_S=new((/num_eq_lat/),integer)
   Tot_N=new(1,double)
   Tot_S=new(1,double)
   AT_S=new((/num_eq_lat/),double)
   AT_N=new((/num_eq_lat/),double)

   if (any(ismissing(countour_S))) then
      e_lat_S= extract_PV@_FillValue
   else
      dsizes_a_S = dimsizes(extr_PV_S)
      dsizes_a_N = dimsizes(extr_PV_N)
      extract_PV_S  = ndtooned(extract_PV(lat|:(floattointeger(jlat/2)),lon|:))     ; reduce to 1D array  
      extract_PV_N  = ndtooned(extract_PV(lat|(floattointeger(jlat/2)):(jlat-1),lon|:))     ; reduce to 1D array

      do i=0,num_eq_lat-1
         if (i.eq.0) then 
            indices_S = ind((extract_PV_S.eq.countour_S(i)))                    
            indices_N = ind((extract_PV_N.eq.countour_N(i)))              
         else
            indices_S = ind((extract_PV_S.ge.countour_S(i-1)).and.(extract_PV_S.lt.countour_S(i)))
            indices_N = ind((extract_PV_N.gt.countour_N(i-1)).and.(extract_PV_N.le.countour_N(i)))
         end if
         if(any(ismissing(indices_S))) then
            print("Your data contains some missing values. Beware.")
            AJ_S(i)= 0.          
         else 
            ind_res_S=ind_resolve(indices_S,dsizes_a_S)
            dim1_S=dimsizes(ind_res_S)
            A1_S=new((/dim1_S(0)/),double)
            lat_p_S=new((/dim1_S(0)/),double)
            lat_m_S=new((/dim1_S(0)/),double)
               do k=0,dim1_S(0)-1
                  if ((lat_new_S(ind_res_S(k,0))).le.-89.9) then		                      
                     lat_m_S(k)=-89.375-0.625
                     lat_p_S(k)=-89.375+0.625	    
                  else			 
                     if ((lat_new_S(ind_res_S(k,0))).ge.-0.1) then
                        lat_p_S(k)=-0.625+0.625
                        lat_m_S(k)=-0.625-0.625
                     else
                        lat_p_S(k)=lat_new_S(ind_res_S(k,0))+1.25
                        lat_m_S(k)=lat_new_S(ind_res_S(k,0))-1.25
                     end if
                  end if 

                  A1_S(k)=2.5*degreetoradians*(abs(sin(lat_p_S(k)*degreetoradians)-sin(lat_m_S(k)*degreetoradians)))
                  ;A1_S_pr=gc_qarea((/lat_p_S(k),lat_p_S(k),lat_m_S(k),lat_m_S(k)/),(/0.,360.,0.,360./))
                  ;  print(A1_S(k))
                  ;  print(A1_S_pr)
               end do      
            AJ_S(i)=sum(A1_S)

            delete(ind_res_S)
            delete(dim1_S)
            delete(lat_p_S)
            delete(lat_m_S)
            delete(A1_S)
         end if

         if(any(ismissing(indices_N))) then
            print("Your data contains some missing values. Beware.")
            AJ_N(i)= 0.        
         else                  
            ind_res_N=ind_resolve(indices_N,dsizes_a_N)   		  
            dim1_N=dimsizes(ind_res_N)                              
            A1_N=new((/dim1_N(0)/),double)
            lat_p_N=new((/dim1_N(0)/),double)
            lat_m_N=new((/dim1_N(0)/),double)                
            do k=0,dim1_N(0)-1                                
               if ((lat_new_N(ind_res_N(k,0))).ge.89.9) then
                  lat_m_N(k)=89.375-0.625
                  lat_p_N(k)=89.375+0.625	    
               else			 
                  if ((lat_new_N(ind_res_N(k,0))).le.0.1) then
                     lat_p_N(k)=0.625+0.625
                     lat_m_N(k)=0.625-0.625               
                  else
                     lat_p_N(k)=lat_new_N(ind_res_N(k,0))+1.25
                     lat_m_N(k)=lat_new_N(ind_res_N(k,0))-1.25
                  end if
               end if 

               A1_N(k)=2.5*degreetoradians*(abs(sin(lat_p_N(k)*degreetoradians)-sin(lat_m_N(k)*degreetoradians)))
               ; A1_N_pr=gc_qarea((/lat_p_N(k),lat_p_S(k),lat_m_N(k),lat_m_N(k)/),(/0.,360.,0.,360./))
            end do

            Num_N(i)=dimsizes(A1_N) 
            AJ_N(i)=sum(A1_N) 

            delete(dim1_N)
            delete(lat_p_N)
            delete(lat_m_N)
            delete(A1_N)
            delete(ind_res_N)
         end if

         AT_S(i)=sum(AJ_S(0:i)) 
         delete(indices_S)
         AT_N(i)=sum(AJ_N(0:i)) 
         delete(indices_N)

         Tot_S=sum(AJ_S(:))
         Tot_N=sum(AJ_N(:))

         if(any(ismissing(Tot_S))) then
            if ((fabs(Tot_S-(2*pi))).gt.0.0001) then
               AT_S=AT_S*(2*pi/Tot_S)
            end if
         end if    
         if(any(ismissing(Tot_N))) then
            if ((fabs(Tot_N-(2*pi))).gt.0.0001) then
               AT_N=AT_N*(2*pi/Tot_N)
            end if
         end if
      end do
   end if  
   ;print(AT_S)
;  print(AT_N)

;**************************Trasformo in Latitudine Equivalente***************************

;****************************************************************************************
   e_lat_S=new((/num_eq_lat/),"double")
   e_lat_N=new((/num_eq_lat/),"double")
   totn_S=new(1,"double")
   totn_N=new(1,"double")

   do i=0,num_eq_lat-1
      e_lat_S(i)=57.29577951*(asin(1-(AT_S(i)/(2*pi))))
      e_lat_N((num_eq_lat)-i-1)=57.29577951*(asin(1-(AT_N(i)/(2*pi))))
      ;print(e_lat_S)
      ;print(countour_S)
   end do
   countour_S!0="Equivalent_latitude"
   countour_S&Equivalent_latitude=e_lat_S
   countour_N!0="Equivalent_latitude"
   countour_N&Equivalent_latitude=e_lat_N

   if (he.eq."SH") then
      return (countour_S)
   end if
   if (he.eq."NH") then
     return (countour_N)
   end if

   ;return (e_lat_N)

   totn_S=sum(Num_S(i))
   totn_N=sum(Num_N(i))

;******************************************************************************************
end



