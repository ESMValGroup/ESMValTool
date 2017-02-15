subroutine sample_events(nx,ny,nt,nyr,monthlysm,smbef,smaft,lon,lat,mn,filedir,daysperyear,sfiledir)
implicit none


!----------------------------------------
! reads in text files containing details of all events defined by global_rain_sm_lgc.f90
! and outputs ascii files containing large control sample of soil moisture anomaly differences from each event configuration
! each filename defines x,y location and month and contains info on event at that pixel 
! (w=wet pixel, d=dry pixels)
! 1: number of days since 31/12/02, 2: rain_w, 3: rain_d, 4: sm_w, 5: ave sm_d, 6: sm_clim_w, 7:sm_clim_d
! 8: box-mean sm, 9:n_d, 10 onwards :pixel numbers for n dry pixels
!------------------------------------------


! Vars that need inputting (single month, all years)
integer,parameter::npd=8
integer,intent(in)::nx,ny,nt,nyr,mn,daysperyear !nyr = 1 later, but needed for testing now
!
real,intent(in),dimension(nyr,nt,ny,nx)::monthlysm
! sm on last day of previous month (bef) and first day of following
! month (aft)
real,intent(in),dimension(nyr,npd,ny,nx)::smbef,smaft
real,intent(in),dimension(nx)::lon
real,intent(in),dimension(ny)::lat
character*150::filedir,sfiledir

integer,parameter::box_del2=1

! output resolution
real,parameter::dlon_out=5.,dlat_out=dlon_out,lon1_out=-180.+dlon_out/2.,lat1_out=-60.+dlat_out/2.
integer,parameter::nxout=360/dlon_out,dyout=120/dlat_out
real::dlat
integer::i,j,l,iop,ii,jj,s,iout,jout,jjout,jjout_old=0,ndom

real::a,b, find_max_contrast,dsm,dsm_max=10.*2*1000.
!
character*150::event_fname
character*2::cmon,ich,jch
integer::ic,rw,rd,smw,smd,smw_clim,smd_clim,n,box_mean_sm,day
integer,dimension((box_del2*2+1)**2)::k
real,dimension(nyr,nt/npd,1-box_del2:ny+box_del2,1-box_del2:nx+box_del2)::sm
real,dimension(1-box_del2:ny+box_del2,1-box_del2:nx+box_del2)::sm_mean,sm_mean_all
integer::n_mean
real,dimension(nyr,1-box_del2:ny+box_del2,1-box_del2:nx+box_del2)::sm_tot_yr,n_yr
real,parameter::undef=-999.
integer,parameter::big=100000*2*2*2*4


!set y dimension of these arrays to 1 and perform calculations on 1 output row at a time
real,dimension(big,1,nxout)::event,sample
integer,dimension(1,nxout)::nevent,nsample

logical::absolute_sm, all_n=.true.,l_exclude_year=.true.
real,dimension(ny,nx)::nsm!=0.!,smclim=0.,nclim=0.

! new vars
real,dimension(npd*3,ny,nx)::soilmoisture !3 diurnal cycles of sm
real,dimension(1-box_del2:ny+box_del2,1-box_del2:nx+box_del2)::smlocal
integer :: d1, d2, dayev, yrev, mnev, yr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CODE STARTS HERE
!
ndom=nt/npd 
!
absolute_sm=.false.
!
dlat = lat(2)-lat(1)


! This code runs for a single month (unlike original version)
write(cmon,'(i2.2)')mn 


!----------------------------
! READ number of events per grid cell per month (nsm)
! written out by global_rain_sm_lgc.f90
!----------------------------

open(2,file=trim(filedir)//'event_output/global_rain_sm_'//cmon//'.gra',&
             status='old', form='unformatted', access='direct', recl=4*nx*ny)

read(2,rec=2) nsm(:,:)
close(2)


! set up some stuff
sm_mean=0. ; n_mean=0  ; sm_mean_all=0. ; sm_tot_yr=0. ; n_yr=0.
!
j=1

nsample=0 ; nevent=0


!-----------------------
! CALCULATE SOIL MOISTURE AT 06-09 LOCAL TIME from the input SM
!-----------------------

do yr=1, nyr
    do day=1, ndom

        ! day indices (d1 = day-1, d2 = day+1, to make 3 diurnal cycles)
        d1 = npd*(day-2) + 1
        d2 = npd*(day+1)

        ! first get 3 diurnal cycles of sm (soilmoisture)
        if(day.eq.1) then
            soilmoisture(1:npd,:,:) = smbef(yr,:,:,:)
            soilmoisture(npd+1:npd*3,:,:) = monthlysm(yr,1:d2,:,:)
        else if(day.eq.ndom) then
            soilmoisture(1:npd*2,:,:) = monthlysm(yr,d1:nt,:,:)
            soilmoisture(npd*2+1:npd*3,:,:) = smaft(yr,:,:,:)
        else
            soilmoisture(1:npd*3,:,:) = monthlysm(yr,d1:d2,:,:)
        endif

        ! then get grid of sm at 06-09 LT
        call sm_local_lgc(ny,nx,npd,soilmoisture,smlocal,lon,box_del2)
        sm(yr,day,:,:) = smlocal

        sm_mean_all = sm_mean_all + smlocal
        n_mean = n_mean+1
        sm_tot_yr(yr,:,:) = sm_tot_yr(yr,:,:) + smlocal

    enddo
enddo


!-------------------
! CALCULATE MONTHLY AND YEARLY MEANS - sm_mean and sm_tot_yr
if(absolute_sm) then
    sm_mean_all=0.
    print*,'setting sm_mean_all=0 ie no climatology to subtract!!!!!!'
else
    do i=1-box_del2,nx+box_del2
        do j=1-box_del2,ny+box_del2
            if(abs((sm(1,1,j,i))+999).ge.epsilon(sm)) then
                sm_mean_all(j,i) = sm_mean_all(j,i)/n_mean
                sm_tot_yr(:,j,i) = sm_tot_yr(:,j,i)/nyr
            else
                sm_mean_all(j,i) = undef
                sm_tot_yr(:,j,i) = undef
            endif
        enddo
    enddo
endif




! LOOP OVER ALL POINTS LOOKING FOR EVENTS
do j=1, ny

!   output data on single output y row at a time to accomodate ncl
!   jjout is the original counter of output row whilst jout in following code is set to 1
!   except where output files named

    ! jjout is the index for the downscaled data
    ! CT version seems to create slight issue at edges (5 points at start, 3 at the end)
    ! Mine is exactly equivalent for now
    jjout = nint( (lat(j) - lat1_out) / dlat_out ) + 1

    ! Confused (again) by need of loop
    if(jjout.ne.jjout_old) then
        do i=1,nxout
            if(nsample(1,i).gt.0) sample(1:nsample(1,i),1,i)=0.
            if( nevent(1,i).gt.0)  event(1:nevent(1,i),1,i)=0. 
        enddo
        nsample=0 ; nevent=0
    endif
!
    do i=1, nx

!
        if(abs(nsm(j,i)).lt.epsilon(nsm(j,i))) cycle ! nsm=0: no event at this point

        !-------------
        ! OPEN EVENT FILE
        !-------------

        open(1,file=trim(event_fname(i,j,mn,filedir)),status='old',iostat=iop)

        if(iop.ne.0) then
            close(1)
            print*,'problem with code!!! missing event file',&
                    trim(event_fname(i,j,mn,filedir))
            stop
        endif
!
        ! loop over all events at single location
        do l=1,nint(nsm(j,i))
!
            read(1,*)ic,rw,rd,smw,smd,smw_clim,smd_clim,box_mean_sm,n
            backspace(1)
            read(1,*)ic,rw,rd,smw,smd,smw_clim,smd_clim,box_mean_sm,n,k(1:n)

            ! calculate day/year from ic
!
            if(daysperyear.eq.360) then
                call icount2date360(ic,dayev,mnev,yrev)
            elseif(daysperyear.eq.365) then
               call icount2date365(ic,dayev,mnev,yrev)
            else
              print*,'daysperyear not set correctly for sample_events.f90',daysperyear ; stop
            endif
!

            !set up sm_mean specifically excluding data in month of event
            if(l_exclude_year.and.nyr.gt.1) then    
!
                do jj=j-box_del2,j+box_del2
                    do ii=i-box_del2,i+box_del2
                        if(abs((sm(1,1,j,i))+999).ge.epsilon(sm)) then
                            sm_mean(jj,ii) = ( sm_mean_all(jj,ii)*n_mean &
                                                - sm_tot_yr(yrev,jj,ii) ) &  
                                            / ( n_mean - nyr )
                        else
                            sm_mean(jj,ii) = undef
                        endif
                    enddo
                enddo

            else
                sm_mean=sm_mean_all
            endif
            

            ! not sure what the following bit does
            dsm=find_max_contrast(sm_mean(j-1:j+1,i-1:i+1))
            if(dsm.gt.dsm_max) then
                print*,'wet point',i,j,sm_mean(j-1:j+1,i-1:i+1),dsm
                cycle
            endif

            ! LGC: changed jj definitions to reflect that lat is not inverted now
            do s=1,n
                ii=mod(k(s)-1,2*box_del2+1)-box_del2 
                jj=int((k(s)-1)/(2*box_del2+1)) - box_del2
                dsm=find_max_contrast(sm_mean(j+jj-1:j+jj+1,i+ii-1:i+ii+1))
                if(dsm.gt.dsm_max) exit
            enddo

            !already checked if dsm.gt.dsm_max?
            if(dsm.gt.dsm_max) then
                print*,'dry point',i,j,i+ii,j+jj,sm_mean(j+jj,i+ii),dsm
                cycle
            endif


            a=0. ; b=0. ! a and b store dry area sm on event day and monthly mean

            ! LGC: same change to jj as above
            do s=1,n
                ii=mod(k(s)-1,2*box_del2+1)-box_del2 
                jj=int((k(s)-1)/(2*box_del2+1)) - box_del2

                a=a+sm(yrev,dayev,j+jj,i+ii)/n
                if(.not.absolute_sm) b=b+sm_mean(j+jj,i+ii)/real(n)
            enddo

            ! checking that soil moisture values read in here agree with output
            ! from previous code (within rounding errors)
            if(nint(100.*sm(yrev,dayev,j,i)).ne.smw .or. nint(a*100).ne.smd) then
                if(abs(100.*sm(yrev,dayev,j,i)-smw).gt.50 .or. abs(100.*a-smd).gt.50) then
                    print*,'possible problems in code... sm on event day ic=',&
                            ic,' day=',dayev,' mn=',mnev,' yr=',yrev
                    print*,'in global_rain_sm.f90 smw,smd=',&
                            smw,smd,' and here (j,i)',j,i, &
                            nint(100.*sm(yrev,dayev,j,i)), &
                            nint(a*100),' with ndry=',n
                        
                    if(abs(smd/100.-a).gt.1.0) then
                        print*,n,k(1:n),nsm(i,j)
                        print*,'stopping...',smd,nint(smd/100.),a ; stop
                    endif
                endif
            endif

! accomodate creation of multiple output boxes within longitude band in here
! NOT SURE IF NECESSARY
            iout = nint( (lon(i) - lon1_out) / dlon_out ) + 1

            jout=1


! event is output array containing event sm contrast (contrast wrt climatology if .not.absolute_sm) 
            nevent(jout,iout) = nevent(jout,iout) + 1
            event(nevent(jout,iout),jout,iout) = sm(yrev,dayev,j,i) - sm_mean(j,i) - (a-b)

            do yr=1, nyr
                do day=1, ndom

                    if (abs(sm(yr,day,j,i)-undef).lt.epsilon(undef)) cycle

                    ! only take sample if all sm data points valid over +/-2 grid boxes


                    ! don't take sample if within +/- 5days of event 
                    ! (note +/-5 slots in array may not be +/-5 days, could be adjacent year)
                    ! FOR NOW DON'T OVERLAP YEARS (MAY CAUSE SOME DISCREPANCIES
                    !if(.not.l_exclude_year .and. abs(dayev-day).le.5 .and. yrev.eq.yr)  cycle
                    ! now, the same as original code
                    if(.not.l_exclude_year) then
                        if(abs(dayev-day).le.5.and.yrev.eq.yr) cycle
                        if((yrev-yr).eq.1.and.(dayev+30-day).le.5) cycle
                        if((yr-yrev).eq.1.and.(day+30-dayev).le.5) cycle
                    endif



                    ! if excluding year then don't take any samples from that year
                    ! AS I HAVE YREV NO NEED FOR COMPLICATED LOOP
                    if(l_exclude_year.and.yrev.eq.yr) cycle

                    a=0. ; b=0.
                
                    do s=1,n
                        ii=mod(k(s)-1,2*box_del2+1)-box_del2 
                        jj=int((k(s)-1)/(2*box_del2+1)) - box_del2

                        if(abs(sm(yr,day,j+jj,i+ii)-undef).gt.epsilon(undef)) then

                            a=a+real(sm(yr,day,j+jj,i+ii)) - sm_mean(j+jj,i+ii)
                            b=b+1.
                        endif
                    enddo

! 5/1/12 introduce logical all_n - if true=>only use cases where all n pixels are defined
! REMOVE ALL_N LOGICAL (ALWAYS TRUE) AND LOOPS?
                    if(all_n) then
                        if(b.lt.n) cycle
                    else
                        if(b.eq.0.) cycle
                    endif
    

! when stratifying cases according to event gridbox mean sm,
! need to calculate gridbox mean sm for samples as well

                    nsample(jout,iout)=nsample(jout,iout)+1
                    sample(nsample(jout,iout),jout,iout) = &
                                real(sm(yr,day,j,i))-sm_mean(j,i)-a/b
!
                enddo  !end dom loop
            enddo ! end yr loop
!

        enddo ! end 'l' loop across event files
        close(1) ! close events file
    enddo   ! end of x loop


    jjout_old=jjout

    !if jjout on next row ne current jjout then output data
    if(jjout.eq.nint((j)*dlat/dlat_out)+1) cycle  

    jout=1
!
    do iout=1,nxout
!
        if(nevent(jout,iout).eq.0) cycle
!
        write(ich,'(i2)') iout ; if(iout.lt.10)ich(1:1)='0'
        write(jch,'(i2)') jjout ; if(jjout.lt.10)jch(1:1)='0'

        open(3,file=trim(sfiledir)//'5x5G_mon'//cmon//'/i'//ich//'j'//jch//'_fort.3')
        open(4,file=trim(sfiledir)//'5x5G_mon'//cmon//'/i'//ich//'j'//jch//'_fort.4')
        do i=1,nsample(jout,iout)
            write(3,*) sample(i,jout,iout)
        enddo

        do i=1,nevent(jout,iout)
            write(4,*) event(i,jout,iout)
        enddo
        close(3) ; close(4)
    enddo

enddo   !end of y loop
!

end
!
!==============================================================
!
function find_max_contrast(sm)
implicit none
!
! takes in 3x3 array
! and returns maximum absolute difference between central point
! and its 8 neighbours 
!
real,dimension(-1:1,-1:1)::sm
real::find_max_contrast,smax,smin
!
smax=maxval(sm) ; smin=minval(sm)
find_max_contrast=max(smax-sm(0,0),sm(0,0)-smin)
return
end
!
!========================================================
!
! Convert ic to day/yr
! Luis' code for 360 day year
!
subroutine icount2date360(ic,day,mn,yr)

implicit none

INTEGER :: ic, day, mn, yr


day = mod(mod(ic,360),30)
if (day.eq.0) day = 30

mn = mod(ic-day, 360)/30 + 1
yr = (ic-day-(mn-2)*30)/360

end
!
!========================================================
!
! Convert ic to day/yr for 365 day year runs
!
subroutine icount2date365(ic,day,mn,yr)

implicit none

INTEGER :: ic, day, mn, yr, icount, mon, ndays, year
integer,parameter::nmon=12
integer,dimension(nmon+1)::day1
data day1/1,32,60,91,121,152,182,213,244,274,305,335,366/
!
icount=0
do year=1,100
  if(ic.le.icount+365) then
    yr=year
    exit
  endif
  icount=icount+365
enddo
!
do mon=1,nmon
  ndays=day1(mon+1)
  if(ic.lt.icount+ndays) then
    mn=mon
    day=ic-icount-day1(mon)+1
    exit
  endif
enddo
!
end

!========================================================
function event_fname(i,j,imon,filedir)
implicit none

! filename for the events output by global_rain_sm_lgc.f90

integer::i,j,imon
character*150::event_fname
character*150::filedir
character*4::cx
character*3::cy
character*2::cmon

write(cx,'(i4)') i
if(i.lt.1000) cx(1:1)='0' 
if(i.lt.100) cx(2:2)='0'
if(i.lt.10) cx(3:3)='0'

write(cy,'(i3)') j
if(j.lt.100) cy(1:1)='0'
if(j.lt.10) cy(2:2)='0'

write(cmon,'(i2)') imon ; if(imon.lt.10) cmon(1:1)='0'

event_fname=trim(filedir)//'event_output/mon'//cmon//'/x'//cx//'_y'//cy//'_mon'//cmon//'.txt'

return 

end

!================================================================
subroutine sm_local_lgc(ny,nx,npd,soilmoisture,smlocal,lon,box_del2)
implicit none

! takes 3 diurnal cycles of soil moisture and loops around longitude to get 
! soil moisture at 06-09LT 
! NB lon must be -180 - 180
! NB2: index of time to extract is hard-coded here
! NB3: smlocal wraps around both dimensions by box_del2

integer :: ny,nx,i,step_sm,npd,box_del2
integer :: t ! time to extract
real :: lon(nx), soilmoisture(npd*3,ny,nx)
real :: smlocal(1-box_del2:ny+box_del2,1-box_del2:nx+box_del2)

t = 3 !6-9LT
do i=1,nx

    step_sm = ceiling(-lon(i)/(360./npd))  !!! belgal: 360 !!!!

    !if(step_sm.lt.1) step_sm=step_sm+8

    smlocal(1:ny,i) = soilmoisture(step_sm+npd+t,:,i)

enddo

! wrap around longitude/latitude
do i=1,box_del2
    smlocal(:,1-i) = smlocal(:,1-box_del2+nx)
    smlocal(:,nx+i) = smlocal(:,i)
    
    smlocal(1-i,:) = smlocal(1-box_del2+ny,:)
    smlocal(ny+i,:) = smlocal(i,:)
enddo


end


