subroutine global_rain_sm(nx,ny,nt,nyr,monthlypr,prbef,praft,monthlysm,smbef,smaft,sm_clim,topo,lon,mn,filedir,daysperyear)
implicit none

! Vars that need inputting (single month, all years)
integer,parameter::npd=8
integer,intent(in)::nx,ny,nt,nyr,mn,daysperyear !nyr = 1 later, but needed for testing now
real,intent(in),dimension(nyr,nt,ny,nx)::monthlypr
real,intent(in),dimension(nyr,nt,ny,nx)::monthlysm
real,dimension(ny,nx)::sm_clim,topo
character*100::filedir

! rainfall and sm on last day of previous month (bef) and first day of following
! month (aft)
real,intent(in),dimension(nyr,npd,ny,nx)::prbef,praft,smbef,smaft
real,intent(in),dimension(nx)::lon
integer:: ic,ic_old,ndom   


integer::k,i,j,ii,jj,n,d1,d2,yr,l

integer,parameter::box_del2=1       !3x3 box (code searches over a box this size
                        !centred on maximum rainfall

integer,parameter::max_events=80  !LGC: HadGEM2 value - check if/how to amend                                                           
! undef now used for both sm and sm_clim (both real variables)
real,parameter::undef=-999.
real,dimension(ny,nx)::sm
real,dimension(ny,nx)::dsm,ndel,nwetter,ndrier,dsm_clim
real,dimension(ny,nx)::pm_rain
real,dimension(npd*3,ny,nx)::rain !this array holds 3 days of rain (t-1,t,t+1)
real,dimension(npd*3,ny,nx)::soilmoisture ! same for soil moisture
real,dimension(npd,ny,nx)::rain_lt!this array holds 1 diurnal cycles of rain based on local time

integer,dimension(ny,nx)::irain_max,irain_min
real::smwet,smdry,smdry_clim,nsmdry,box_mean_sm,nbox
integer,dimension(ny,nx)::event
logical,dimension(ny,nx)::nearby_event
!
!max(min) rain per pixel between am1 and am2 (pm1 and pm2)
real::max_am_rain=.1,min_pm_rain=3.
!threshold for topographic variability in models
real::topo_max=500.
!
integer::am1,am2,pm1,pm2
!
! event_op is event data: days since start, rain, rain_dry, sm, sm_dry, smcl, smcl_dry, box_mean_sm, ndry
! bevent_op is idry(ndry)
!
integer,parameter::nvar=9,nbvar=8 ! eqn gives an error on f2py: (2*box_del2+1)**2-1 !5x5 box
!
integer,dimension(ny,nx,max_events,nvar)::event_opD
integer,dimension(ny,nx,max_events,nbvar)::bevent_opD
!
! characters for filenames
character*4::cx
character*3::cy
character*2::cmon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CODE STARTS HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialize ndel
!
do j=1,ny
    do i=1,nx
        dsm(j,i)=0.
        ndel(j,i)=0.
        nwetter(j,i)=0.
        ndrier(j,i)=0.
        dsm_clim(j,i)=0.
    enddo
enddo

! note time from python script starts at 0300, not 0000
!
! these numbers define the "morning" and "afternoon" windows
! within the diurnal cycle
! they cannot be changed without changing the code in sample_events.f90
! and in subroutine sm_local_lgc
!
am1=2 ; am2=2 !6 - 9am
pm1=3 ; pm2=5 !9 - 18

event_opD=0 ; event=0  ; bevent_opD=0

ic_old=-999
!
ndom=nt/npd 
!

!!!!!!!!!!!!!!!
! LOOP OVER YEARS AND DAYS
!!!!!!!!!!!!!!!

do yr=1,nyr 
    do k=1,ndom

        ! day counter
!        call icount(yr,mn,k,ic)
        if(daysperyear.eq.360) then
          call icount360(yr,mn,k,ic) 
        elseif(daysperyear.eq.365) then
          call icount365(yr,mn,k,ic)
        else
          print*,'daysperyear not set correctly for global_rain_sm.f90 ',daysperyear ; stop
        endif

        ! day indices (d1 = k-1, d2 = k+1, to make 3 diurnal cycles)
        d1 = npd*(k-2) + 1
        d2 = npd*(k+1)


        ! Calculate morning soil moisture as local time (6-9)
        ! first get 3 diurnal cycles of sm (soilmoisture)
        if(k.eq.1) then
            soilmoisture(1:npd,:,:) = smbef(yr,:,:,:)
            soilmoisture(npd+1:npd*3,:,:) = monthlysm(yr,1:d2,:,:)
        else if(k.eq.ndom) then
            soilmoisture(1:npd*2,:,:) = monthlysm(yr,d1:nt,:,:)
            soilmoisture(npd*2+1:npd*3,:,:) = smaft(yr,:,:,:)
        else
            soilmoisture(1:npd*3,:,:) = monthlysm(yr,d1:d2,:,:)
        endif

        ! then get grid of sm at 06-09 LT
        call sm_local_lgc(ny,nx,npd,soilmoisture,sm,lon)

        ! read in rainfall
        ! need rain from previous evening and following morning
        ! for first and last day of month, use prbef and praft which gives rain from other months
        if(k.eq.1) then
            rain(1:npd,:,:) = prbef(yr,:,:,:)
            rain(npd+1:npd*3,:,:) = monthlypr(yr,1:d2,:,:)
        else if(k.eq.ndom) then
            rain(1:npd*2,:,:) = monthlypr(yr,d1:nt,:,:)
            rain(npd*2+1:npd*3,:,:) = praft(yr,:,:,:)
        else
            rain(1:npd*3,:,:) = monthlypr(yr,d1:d2,:,:)
        endif

        ! create single diurnal cycle of rain using local time
        call rain_local_lgc(ny,nx,npd,rain,rain_lt,lon,am2)


        ! find rainfall maxima
        ! NB returns all potential rain maxima regardless of sm or topography.
        ! Maxima affected by these two are excluded below, but can still affect
        ! other maxima becoming events.
        call find_rain_max(ny,nx,npd,rain_lt,undef,irain_max,irain_min, &
                    am1,am2,pm1,pm2,max_am_rain,min_pm_rain,box_del2,pm_rain, lon)
    !
        nearby_event=.false.

        do j=1,ny
            big_loop: do i=1,nx

                ! cycle if certain conditions not met 
                if(irain_max(j,i).eq.0) cycle
!
!               use sm_clim rather than sm to avoid ocean points
!               also applied later when removing ocena points within box_del2
!
                if(sm_clim(j,i).lt.-100.) cycle
!
                if(topo(j,i).gt.topo_max) cycle    !exclude rain maxima over high topography
        
                box_mean_sm=0. ; nbox=0. ; smdry=0. ; nsmdry=0. ; smdry_clim=0.

                n=irain_max(j,i)
                
                ! smwet = rain over rainfall peak
                smwet = real(sm(j,i))
        
                if(smwet.lt.0.) smwet=0.     !avoids small -ve sm
        
        
        !       calculate box mean sm and reject cases where 50% or less of box contains valid sm
        
                do jj=j-box_del2,j+box_del2
                    do ii=i-box_del2,i+box_del2

                    ! 2 cases for rejection (nearby_event, or undefined sm/sm_mask over rain_minimum)

                        if(nearby_event(jj,ii)) cycle big_loop
                    
                        if(sm_clim(jj,ii).lt.-100.) cycle big_loop
!
                        !exclude rain maxima where one of rain minima has either undefined sm or is a flagged sm pixel
                        if(irain_min(jj,ii).eq.n) then 
                                                                
                            if(sm(jj,ii).lt.-100.) cycle big_loop
        
                            if(topo(jj,ii).gt.topo_max) cycle big_loop
                        endif


                        if(sm(jj,ii).gt.-100.) then
                            box_mean_sm=box_mean_sm+(sm(jj,ii))/100. 
                            nbox=nbox+1.
                        endif

                    enddo
                enddo

                ! reject if 50% or less of box contains valid sm
                if(nbox.le.0.5*((2*box_del2+1)**2)) cycle

                do jj=j-box_del2,j+box_del2
                    do ii=i-box_del2,i+box_del2
                        if(irain_min(jj,ii).eq.n) then

                            smdry = smdry + real(sm(jj,ii)) 
                            if(smdry.lt.0.) smdry=0.     ! avoids small -ve sm
                            nsmdry = nsmdry + 1.

                            if(nint(nsmdry).eq.1.) then
                                event(j,i) = event(j,i) + 1
                                event_opD(j,i,event(j,i),2) = nint(pm_rain(j,i)*100.)
                                event_opD(j,i,event(j,i),3) = nint(pm_rain(jj,ii)*100.)
                            endif

                            bevent_opD(j,i,event(j,i),nint(nsmdry)) = &
                                    (jj-(j-box_del2))*(box_del2*2+1)+ii-(i-box_del2)+1

                            !to avoid running out of memory, do calculation as real
                            !and set event_op(7) only when correct average computed
                            smdry_clim = smdry_clim + 100.*sm_clim(jj,ii)   
                                                                
                        endif
                    enddo
                enddo
    
                if(nsmdry.eq.0.) cycle
    !
                smdry=smdry/nsmdry
                box_mean_sm=box_mean_sm/nbox
                dsm(j,i) = dsm(j,i) + smwet - smdry
                ndel(j,i) = ndel(j,i) + 1.
                if(smwet.gt.smdry) then
                    nwetter(j,i) = nwetter(j,i) + 1.
                elseif(smdry.gt.smwet) then
                    ndrier(j,i) = ndrier(j,i) + 1.
                endif
        !
                event_opD(j,i,event(j,i),1) = ic
                event_opD(j,i,event(j,i),4) = nint(smwet*100.)
                event_opD(j,i,event(j,i),5) = nint(smdry*100.)
                event_opD(j,i,event(j,i),6) = nint(100.*sm_clim(j,i))
                event_opD(j,i,event(j,i),7) = nint(smdry_clim/nsmdry)   
        !
                event_opD(j,i,event(j,i),8) = nint(box_mean_sm*100.)
                event_opD(j,i,event(j,i),9) = nint(nsmdry)
                dsm_clim(j,i) = dsm_clim(j,i) + (smwet-sm_clim(j,i)) - (smdry-event_opD(j,i,event(j,i),7)/100.)
    

    !
                nearby_event(j,i)=.true.
    !
                 
            enddo big_loop
        enddo
   
        ic_old=ic
    enddo
enddo


!!!!!!!!!!!!
! WRITE STUFF OUT
!!!!!!!!!!!!


write(cmon,'(i2.2)')mn 
!
! this monthly file contains details about the number of events per grid cell
! which is used sample_events.f90, as well as additional diagnostics
! for checking the code at an interim point
!
open(4,file=trim(filedir)//'event_output/global_rain_sm_'//cmon//'.gra',&
             status='replace',access='direct', form='unformatted', recl=4*nx*ny)

write(4,rec=1)dsm
write(4,rec=2)ndel
write(4,rec=3)nwetter
write(4,rec=4)ndrier
write(4,rec=5)dsm_clim
close(4)
dsm=0. ; ndel=0. ; nwetter=0. ; ndrier=0. ; dsm_clim=0.
!
! output event ascii files
!
! for every grid point that contains at least one event, create an ascii
! file with the key characteristics of that event. These are read by sample_events.f90
!
do j=1,ny
    do i=1,nx
        if(event(j,i).gt.0) then
            ii=i
            jj = j
            write(cx,'(i4.4)') ii 
            write(cy,'(i3.3)') jj 
            open(20,file=trim(filedir)//&
                    'event_output/mon'//cmon//'/x'//cx//'_y'//cy//'_mon'//cmon//'.txt')

            do k=1,event(j,i)
                l=event_opD(j,i,k,nvar)
                write(20,*) event_opD(j,i,k,:),bevent_opD(j,i,k,1:l)
            enddo

            close(20) 
        endif
    enddo
enddo
!
return
end


!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

subroutine rain_local_lgc(ny,nx,npd,rain,rain_lt,lon,am2)
implicit none

! takes entire rainfall data-set and loops around longitude to get 
! a diurnal cycle of rainfall based on local time
! NB lon must be -180 - 180

integer :: ny,nx,am2,i,k,step_sm,npd
real :: lon(nx), rain(npd*3,ny,nx), rain_lt(npd,ny,nx)

do i=1,nx

    step_sm = ceiling(-lon(i)/(360./npd))+am2

    if(step_sm.lt.1) step_sm=step_sm+npd

    do k=1,npd
        rain_lt(k,:,i) = rain(step_sm-am2+npd+k,:,i)
    enddo
enddo

return

end

!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

subroutine sm_local_lgc(ny,nx,npd,soilmoisture,sm,lon)
implicit none

! takes 3 diurnal cycles of soil moisture and loops around longitude to get 
! soil moisture at 06-09LT 
! NB lon must be -180 - 180
! NB2: index of time to extract is hard-coded here

integer :: ny,nx,i,step_sm,npd
integer :: t ! time to extract
real :: lon(nx), soilmoisture(npd*3,ny,nx), sm(ny,nx)

t = 3 !6-9LT
do i=1,nx

    step_sm = ceiling(-lon(i)/(360./npd))  !!!! belgal: 360!!!

    sm(:,i) = soilmoisture(step_sm+npd+t,:,i)

enddo

return

end

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

subroutine find_rain_max(ny,nx,npd,rain,undef,irain_max,irain_min, &
                   am1,am2,pm1,pm2,max_am_rain,min_pm_rain,box_del2,pm_rain, lon)
implicit none
!
! takes sequence of 3hourly rain images and finds local maxima/minima in rain between pm1 and pm2
! having first removed rain-affected areas between am1 and am2

integer::ny,nx,npd,i,j,k,box_del2,ii,jj,nmax
real,dimension(npd,ny,nx)::rain
integer,dimension(ny,nx)::irain_max,irain_min
real::max_am_rain,min_pm_rain,local_min,rain_step,undef
logical,dimension(ny,nx)::max1
real,dimension(ny,-box_del2+1:nx+box_del2)::pm_rain360  !pm_rain array expanded to provide continuous cover across x=1,nx
integer::am1,am2,pm1,pm2
real,dimension(ny,nx)::am_rain,pm_rain  !am_rain(j,i) can take 2 values per day, 1 for each overpass
real,dimension(nx)::lon

am_rain=0. ; pm_rain=0. ; max1=.false. ; irain_max=0 ; irain_min=0
rain_step=24./npd !conversion factor from rainfall rate in mm/hr over 3hours to mm in 3 hours
! first remove data where there is am rain>threshold
! and find total pm rain for remaining cases

do j=1,ny
  do i=1,nx

    am_loop : &
    do k=am1,am2
      if((rain(k,j,i)).lt.-100.) then
        am_rain(j,i)=undef
        exit am_loop
      endif
            
      am_rain(j,i) = am_rain(j,i) + rain(k,j,i)*rain_step

      
    enddo am_loop
!
    if((am_rain(j,i)).lt.-100..or.am_rain(j,i).gt.max_am_rain) then
      pm_rain(j,i)=undef
    else
      pm_loop : &
      do k=pm1,pm2
        if((rain(k,j,i)).lt.-100.) then
          pm_rain(j,i)=undef
          exit pm_loop
        endif
        pm_rain(j,i) = pm_rain(j,i) + rain(k,j,i)*rain_step
      enddo pm_loop
    endif
  enddo
enddo
!
! at this point all pixels with at least 1 valid sm_mask within the box have had their
! am and pm rain totals computed, regardless of whether the pixel itself has valid sm

nmax=0
!
! set-up pm_rain360
do i=1-box_del2,nx+box_del2
  if(i.ge.1.and.i.le.nx) pm_rain360(:,i)=pm_rain(:,i)    !!!! belgal: 360!!! 
  if(i.lt.1)  pm_rain360(:,i)=pm_rain(:,i+nx)
  if(i.gt.nx) pm_rain360(:,i)=pm_rain(:,i-nx)
enddo
!
! for all non-zero pm rain pixels, find local maximum
! no longer reject maxima at this point in regions of strong topography
!
do j=1+box_del2,ny-box_del2
  do i=1,nx
!
    if(pm_rain(j,i).le.min_pm_rain) cycle
    
! remove if there is a temporal discontinuity (i.e. around 45, 90 etc.)
    if(floor(lon(i)/45.).ne.floor(lon(i+box_del2)/45.)) cycle
    if(floor(lon(i)/45.).ne.floor(lon(i-box_del2)/45.)) cycle

    max1(j,i)=.true.

    
    do jj=j-box_del2,j+box_del2
!if max1 gets set to false within the jj loop then cycle remaining part of loop
      if(.not.max1(j,i)) cycle
      do ii=i-box_del2,i+box_del2
        if(ii.eq.i.and.jj.eq.j) cycle   
        if(pm_rain360(jj,ii).ge.pm_rain(j,i).or.(pm_rain360(jj,ii)).lt.-100.) then   !!!! belgal: 360!!!
            max1(j,i)=.false.
        endif
      enddo
    enddo
  enddo
enddo

! for each local maximum, find all local minima and label each with a counter (nmax)
!
do j=1+box_del2,ny-box_del2
  do i=1,nx
    if(max1(j,i)) then
      nmax=nmax+1
      irain_max(j,i)=nmax
      local_min=minval(pm_rain360(j-box_del2:j+box_del2,i-box_del2:i+box_del2))
      do jj=j-box_del2,j+box_del2
        do ii=i-box_del2,i+box_del2
          if(pm_rain360(jj,ii).eq.local_min) then
            if(ii.ge.1.and.ii.le.nx) then
              irain_min(jj,ii)=nmax
            else
              if(ii.lt.1) irain_min(jj,ii+nx)=nmax
              if(ii.gt.nx) irain_min(jj,ii-nx)=nmax
            endif
          endif
!
        enddo
      enddo
    endif
  enddo
enddo

!
return
end
!
!-------------------------------------------------
!
! Luis' code for 360 day year
!
subroutine icount360(yr,mon,day,ic)

implicit none

INTEGER, INTENT(IN) :: yr,mon,day
INTEGER, INTENT(OUT) :: ic

ic = yr*360 + (mon-1)*30 + day

return
end
!
!-------------------------------------------------
!
! code for 365 day year
!
subroutine icount365(yr,mon,day,ic)
!
! returns day number since Jan 1st of first year
! for 365 day year
!
implicit none
!
INTEGER, INTENT(IN) :: yr,mon,day
INTEGER, INTENT(OUT) :: ic
integer,parameter::nmon=12
integer,dimension(nmon)::day1
data day1/1,32,60,91,121,152,182,213,244,274,305,335/
!
ic = (yr-1)*365 + day1(mon)-1 + day
!
return
end
