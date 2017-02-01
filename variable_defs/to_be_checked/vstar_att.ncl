;
; Requires: PS:C2Ms,ua:C3M,va:C3M,OMEGA:C3M,ta:C3M,VTH3d:C3M,UV3d:C3M,UW3d:C3M
; ======================================================================
;
; NOT WORKING!!!

;
; Created February 7, 2009
; by J. H. Richter and R. R. Garcia
;
; Calculate the TEM budget in log-pressure coordindates from h0 output files
; One time sample in an input file is assummed (each h0 file contains one monthly mean)
;


info = True
info@derived = True
info@long_name = "wstar vstar"

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

  if (field_number.eq."T2Ms") then
    va_all = read_data(index,"va","T3M")
    ta_all = read_data(index,"ta","T3M")
    ua_all = read_data(index,"ua","T3M")
    vth_all = read_data(index,"VTH3d","T3M")
    uv_all = read_data(index,"UV3d","T3M")
    uw_all = read_data(index,"UW3d","T3M")
    omega_all = read_data(index,"OMEGA","T3M")
    ps_all = read_data(index,"PS","T2Ms")
  end if
  if (field_number.eq."T2Is") then
    va_all = read_data(index,"va","T3I")
    ta_all = read_data(index,"ta","T3I")
    ua_all = read_data(index,"ua","T3I")
    vth_all = read_data(index,"VTH3d","T3I")
    uv_all = read_data(index,"UV3d","T3I")
    uw_all = read_data(index,"UW3d","T3I")
    omega_all = read_data(index,"OMEGA","T3I")
    ps_all = read_data(index,"PS","T2Is")
  end if
    if (field_number.eq."T2Ms") then
      p0 = read_from_file(index,"ua","p0","T3M")
      hyai = read_from_file(index,"ua","hyai","T3M")
      hybi = read_from_file(index,"ua","hybi","T3M")
      hyam = read_from_file(index,"ua","hyam","T3M")
      hybm = read_from_file(index,"ua","hybm","T3M")
    end if
    if (field_number.eq."T2Is") then
      p0 = read_from_file(index,"ua","p0","T3I")
      hyai = read_from_file(index,"ua","hyai","T3I")
      hybi = read_from_file(index,"ua","hybi","T3I")
      hyam = read_from_file(index,"ua","hyam","T3I")
      hybm = read_from_file(index,"ua","hybm","T3I")
    end if

; define dimension 

 dims=dimsizes(ua_all)
 print(dim(0)+"=lon")
 print(dim(1)+"=lat")
 print(dim(2)+"=lev")
 nx=dims(0)
 ny=dims(1)
 nz=dims(2)

;-----------------------------------------------------------------------
; Interpolate zonally averaged  VTH, UV, UW to midpoints
; These are already zonally averaged but are stored as 3D arrays because
;  of the inflexible CAM3 outout structure 
; Values at longitude "0" are interpolated below, but all longitudes
;  actually contain the same values
;-----------------------------------------------------------------------

vthbar = new((/ny,nz/),"float","No_FillValue")
uvbar  = new((/ny,nz/),"float","No_FillValue")
uwbar  = new((/ny,nz/),"float","No_FillValue")

do j = 0, ny-1 

;caf   vthbar(j,*)= interpol(vth(0,j,*), ilev, lev)
;caf   uvbar(j,*) = interpol(uv(0,j,*),  ilev, lev)
;caf   uwbar(j,*) = interpol(uw(0,j,*),  ilev, lev)

   vthbar(j,*)= ftcurv(ilev,vth(0,j,*),  lev)
   uvbar(j,*) = ftcurv(ilev,uv(0,j,*),   lev)
   uwbar(j,*) = ftcurv(ilev,uw(0,j,*),   lev)

; Mask out all regions with unreasonable values
; including after interpolation

   fmask = ind(vth(0,j,*).eq.1.e35) 	; 1.e35 was stored at (lat,lev) where VTH, etc.,
					; could not be computed because too many points
					; were "below ground" at the given pressure
   if (fmask(0).gt.-1) then 
       vthbar(j,fmask) = -999
       uvbar(j,fmask) = -999
       uwbar(j,fmask) = -999
       fmaskm1 = min(fmask)-1
       vthbar(j,fmaskm1) = -999 
       uvbar(j,fmaskm1) = -999
       uwbar(j,fmaskm1) = -999

   end if

end do

;-------------------------------------------------------------------
; Interpolate U, V, OMEGA, T to pressure co-ordinates
; Points "below ground" will have values of -999
;-------------------------------------------------------------------

varout=vinth2p(U,hyam,hybm,lev,ps,1,p0,1,True)
  delete(U)
  U=varout
  delete(varout)
varout=vinth2p(V,hyam,hybm,lev,ps,1,p0,1,True)
  delete(V)
  V=varout
  delete(varout)
varout=vinth2p(OMEGA,hyam,hybm,lev,ps,1,p0,1,True)
  delete(OMEGA)
  OMEGA=varout
  delete(varout)
varout=vinth2p(ta,hyam,hybm,lev,ps,1,p0,1,True)
  delete(ta)
  ta=varout
  delete(varout)

;intop,  hyam, hybm, P0, lev, U, PS, nx, ny, nz, varout
;U = varout

;intop,  hyam, hybm, P0, lev, V, PS, nx, ny, nz, varout
;V = varout

;intop,  hyam, hybm, P0, lev, OMEGA, PS, nx, ny, nz, varout
;OMEGA = varout

;intop,  hyam, hybm, P0, lev, T, PS, nx, ny, nz, varout
;T = varout

; ---------------------------------------------------------------
; Define some quantities needed for the TEM calculations
; ---------------------------------------------------------------

; Pressure (in Pa, same units as reference pressure P0)
  P = 100.*lev

; log-Pressure height (m)
  H = 7.e3		; scale height
  zp = log(P0/P) * H

; Density (kg/m3) on log-pressure surface
  rho = P/(9.81*H)

; ---------------------------------------------------------------
; Calculate potential temperature (masking it the same as T)
; ---------------------------------------------------------------

  th = new((/nx,ny,nz/),"float","No_FillValue")
  kappa = 2./7.
  do  iz = 0, nz-1  
    th(*,*,iz) = T(*,*,iz) * (P0/P(iz))^kappa
  end do
  mask     = ind(T.eq.-999) 
  if (mask(0).ne.-1) then 
    th(mask) = -999
  end if

; ----------------------------------------------------------------------
; Calculate the vertical derivative of theta at every point and then 
; take the zonal average.  This produces a more reasonable calculation 
; of d(th)/thp above the mountain regions than calculating the vertical
; gradient of thbar
; ----------------------------------------------------------------------

thz = add((/nx,ny,nz/),"float","No_FillValue"))

do j = 0, ny-1 	; local d(th)/dz
   do i = 0, nx-1 

     mask22 = ind(th(i,j,*).ne.-999)
;caf     thz(i,j,mask22) = deriv( zp(mask22), th(i,j,mask22) )
     thz(i,j,mask22) = center_finite_diff( zp(mask22), th(i,j,mask22),False,0)

     mask23 = ind(th(i,j,*).eq.-999)
     if(mask23(0).ne.-1) then 
       thz(i,j,mask23) = -999
     end if

   end do
end do

thzbar = add((/ny,nz/),"float","No_FillValue")

do j = 0, ny-1 	; zonal-mean d(th)/dz
   do k = 0, nz-1 

       mask = ind(th(*,j,k).ne.-999)
       npoints = dimsizes(mask)

       if(npoints.ge.nx/3) then  
          thzbar(j,k) = avg(thz(mask,j,k))
       else
         thzbar(j,k) = -999 
       end if

   end do
end do

;-------------------------------------------------------------------
; Calculate Zonally Averaged  U, V, OMEGA, TH
;-------------------------------------------------------------------

ubar = add((/ny,nz/),"float","No_FillValue")
vbar = add((/ny,nz/),"float","No_FillValue")
wbar = add((/ny,nz/),"float","No_FillValue")
thbar = add((/ny,nz/),"float","No_FillValue")

do j = 0, ny-1 
   do k = 0, nz-1 

       mask = ind(u(*,j,k).ne.-999)
       npoints = dimsizes(mask)

       if (npoints.ge.nx/3) then  
           ubar(j,k)=avg(u(mask,j,k))
           vbar(j,k)=avg(v(mask,j,k))
           wbar(j,k)=avg(-omega(mask,j,k)*H/P(k)) ; omega is now converted to w
           thbar(j,k)=avg(th(mask,j,k))

       else
           ubar(j,k)   = -999
           vbar(j,k)   = -999
           wbar(j,k)   = -999
           thbar(j,k)  = -999

       end if

   end do
end do

;---------------------------------------
; Calculate Vertical Derivative of ubar
;---------------------------------------

uzbar   = add((/ny,nz/),"float","No_FillValue")

do j = 0, ny-1 

   mask2 = ind(ubar(j,*).ne.-999)
   mask3 = ind(ubar(j,*).eq.-999)

   uzbar(j,mask2) = center_finit_diff( zp(mask2), ubar(j,mask2),False,0)
   if (mask3(0).ne.-1) then 
     uzbar(j,mask3) = -999
   end if

end do

;-------------------------------------------------------------------
; Apply the same masks to all variables
;-------------------------------------------------------------------

fmask = ind(uvbar.eq.-999)

ubar(fmask)   = -999
vbar(fmask)   = -999
wbar(fmask)   = -999
thbar(fmask)  = -999
thzbar(fmask) = -999
uzbar(fmask)  = -999

mask = ind(ubar.eq.-999)

uvbar(mask)  = -999
uwbar(mask)  = -999
vthbar(mask) = -999

;print, mask 
;print, fmask

;-------------------------------------------------------------------
; Start of the TEM calculation
;-------------------------------------------------------------------

Fphi    = add((/ny,nz/),"float","No_FillValue")      ; Merdional component of EP flux
Fphiphi = add((/ny,nz/),"float","No_FillValue")	     ; d(Fphi)/dphi
Fz      = add((/ny,nz/),"float","No_FillValue")      ; Vertical component of EP flux 
Fzz     = add((/ny,nz/),"float","No_FillValue")      ; d(Fz)/dz
DELF    = add((/ny,nz/),"float","No_FillValue")      ; EP flux divergence
vres    = add((/ny,nz/),"float","No_FillValue")      ; Meridional component of residual velocity
wres    = add((/ny,nz/),"float","No_FillValue")      ; Vertical component of residual velocity

; latitudes in radians
  rlats = lat*!dtor

; a*cos(phi)
  ac = 6.37e6*cos(rlats)           

; rho*a*cos(phi)
  rac = add((/ny,nz/),"float","No_FillValue")
  do k = 0, nz-1 
    rac(*,k) = rho(k)*ac
  end do

; Coriolis parameter (f)
  f = 2.*(2*3.14159/86400.)*sin(rlats)

; ----------------------------------------
; Calculate meridional EP flux component
; ----------------------------------------

Fphi = rac * (uzbar * vthbar/thzbar - uvbar)
mask = ind(ubar.eq.-999)
if(mask(0).ne.-1) then 
  Fphi(mask) = -999
end if

; ----------------------------------------
; Calculate vertical EP flux component
; ----------------------------------------

temp = add(ny,"float","No_FillValue")

do k = 0, nz-1 

   mask2 = ind(ubar(*,k).ne.-999)
   mask3 = ind(ubar(*,k).eq.-999)

   temp(mask2) = center_finite_diff( rlats(mask2), ubar(mask2,k)*cos(rlats(mask2)) ,False,0)
   if (mask3(0).ne.-1) then 
     temp(mask3) = -999
   end if

   Fz(mask2,k) = rac(mask2,k) * ( (f - temp(mask2)/ac )  $
               * vthbar(mask2,k)/thzbar(mask2,k) - uwbar(mask2,k))
   if (mask3(0).ne.-1) then 
     Fz(mask3,k) = -999
   end if

end do 

; -----------------
; Calculate DELF
; -----------------

temp = add(ny,"float","No_FillValue")

do k = 0, nz-1 

   mask2 = ind(ubar(*,k).ne.-999)
   mask3 = ind(ubar(*,k).eq.-999)

   temp(mask2) = center_finite_diff( rlats(mask2), Fphi(mask2,k)*cos(rlats(mask2)),False,0 )
   if (mask3(0).ne.-1) then 
     temp(mask3) = -999
   end if

   Fphiphi(mask2,k) = temp(mask2)/ac
   if (mask3(0).ne.-1) then 
     Fphiphi(mask3,k) = -999
   end if

end do

do j = 0, ny-1 

   mask2 = ind(ubar(j,*).ne.-999)
   mask3 = ind(ubar(j,*).eq.-999)
   Fzz(j,mask2) = center_finite_diff( zp(mask2), Fz(j,mask2),False,0 )
   if (mask3(0).ne.-1) then 
     Fzz(j,mask3) = -999
   end if

end do

DELF = (Fphiphi + Fzz) / rac   ; in m/s/s
mask = ind(ubar.eq.-999)
if(mask(0).ne.-1) then 
  DELF(mask) = -999
end if


; --------------------------------
; Calculate Residual Velocities
; --------------------------------

temp = add(nz,"float","No_FillValue")

do j = 0, ny-1 

   mask2 = ind(ubar(j,*).ne.-999)
   mask3 = ind(ubar(j,*).eq.-999)

   temp(mask2) = center_finite_diff( zp(mask2), rho(mask2)*vthbar(j,mask2)/thzbar(j,mask2),False,0 )
   if (mask3(0).ne.-1) then 
     temp(mask3) = -999
   end if
   vres(j,mask2) = vbar(j,mask2)-(1./rho(mask2))*temp(mask2)
   if (mask3(0).ne.-1) then 
     vres(j,mask3) = -999
   end if

end do

temp = add(ny,"float","No_FillValue")

do k = 0, nz-1 

   mask2 = ind(ubar(*,k).ne.-999)
   mask3 = ind(ubar(*,k).eq.-999)

   temp(mask2) = center_finite_diff( rlats(mask2), cos(rlats(mask2))*vthbar(mask2,k)/thzbar(mask2,k),False,0 )
   if (mask3(0).ne.-1) then 
     temp(mask3) = -999
   end if
   wres(mask2,k) = wbar(mask2,k) + temp(mask2)/ac(mask2)
   if (mask3(0).ne.-1) then 
     wres(mask3,k) = -999
   end if

end do


; ----------------------------------------------------------------
  pro intop,  hyam, hybm, P0, lev, varin, PS, nx, ny, nz, varout
; ----------------------------------------------------------------

; Interpolates variables to pressure grid corresponding to model midpoints, lev
; Points below ground will have values of -999

varout = add((/nx,ny,nz/),"float","No_FillValue")

do ix = 0, nx-1 
  do iy = 0, ny-1 

    press = (hyam(*)*P0 + hybm(*)*PS(ix,iy))/100.       ; actual pressure at (lon=ix, lat=iy)

;caf    varout(ix,iy,*) = interpol(varin(ix,iy,*), press(*), lev(*))
    varout(ix,iy,*) = ftcurv(press(*),varin(ix,iy,*),  lev(*))

    mask = ind(lev.gt.press(nz-1))               	; mask out locations below ground
    if (mask(0).gt.-1) then 
      varout(ix,iy,mask) = -999
    end if

  end do
end do

c=stringtochar(field_number)
c(1)=integertochar(50) ; "2"
c(3)=integertochar(115) ; "s"
field_number=chartostring(c)
print(field_number)
delete(c)
return(data_new)


end

