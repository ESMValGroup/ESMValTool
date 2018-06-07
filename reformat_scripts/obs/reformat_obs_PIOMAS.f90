!#############################################################################
! REFORMAT SCRIPT FOR PIOMAS ARCTIC SEA ICE THICKNESS REANALYSES
!#############################################################################
!
! Tier
!    2 (freely available data set other than obs4MIPs and ana4MIPs)
!
! Source
!    ftp://pscftp.apl.washington.edu/zhang/PIOMAS/data/v2.1/heff/
!    Reference:
!        Zhang, Jinlun and D.A. Rothrock: Modeling global sea ice with a thickness
!        and enthalpy distribution model in generalized curvilinear coordinates,
!        Mon. Wea. Rev. 131(5), 681-697, 2003.
!
! Last access
!    05/2018
!
! Download and processing instructions
!    Download all files from above website.
!    Download 'grid.dat', 'grid.dat.pop' and 'io.dat_360_120.output' from 
!      http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid
!    Adapt 'inpath' & 'outpath' in this script to your system.
!    Process this script (reformat_obs_PIOMAS.f90).
!
! Caveats
!    Fortran netCDF libraries required for compilation.
!
! Modification history
!    20180514-A_laue_ax: written.
!
!############################################################################

program convert_piomas

   use netcdf

   implicit none

   character (len = *), parameter :: inpath = &
      '/export/pa_data01/ESMVal/obs/RAW/Tier2/PIOMAS/'
   character (len = *), parameter :: outpath = &
!      '/export/pa_data02/ESMVal/obs/Tier2/PIOMAS/'
      './'

   integer, parameter :: nx1 = 360
   integer, parameter :: ny1 = 120
   integer, parameter :: nx = nx1
   integer, parameter :: ny = ny1
   integer, parameter :: imt = nx1
   integer, parameter :: jmt = ny1

   integer, parameter :: nyear1 = 1978
   integer, parameter :: nyear2 = 2017

   real, dimension(imt, jmt) :: heff, uice, vice
   real, dimension(imt, jmt) :: clon, clat
   real, dimension(imt, jmt) :: ulat, ulon, HTN, HTE
   real, dimension(imt, jmt) :: HUS, HUW, angle, dxt, dyt
   real, dimension(imt, jmt) :: area

   integer, dimension(imt, jmt) :: kmt

   integer, dimension(imt) :: iVar
   integer, dimension(jmt) :: jVar

   character *80 fopen(5), f1, f2, f3
   character *4 cyear(1900:2100), cyear1(1900:2100)
   character *12 char
   integer SLEN

   integer :: i, j, t
   integer :: imon, iyear

   double precision, dimension(1) :: time, time_ref

   integer status
   integer :: ncID
   integer :: unlimitedDimID
   integer :: LonDimID, LatDimID, TimeDimID
   integer :: LonVarID, LatVarID, TimeVarID, outVarID, outAreaID
   integer :: iDimID, jDimID, iVarID, jVarID
   integer :: oldfillmode

   character (len = 256) :: outfile
   character (len = 9) :: period

   real, parameter :: fillval = 1.e20

   f2 = inpath // 'heff.H'

   ! read lon and lat for scalar fields (like sea ice thickness and concentration)
   open(20, file = inpath // 'grid.dat')
   read(20,'(10f8.2)') ((clon(i, j), i = 1, nx1), j = 1, ny1)
   read(20,'(10f8.2)') ((clat(i, j), i = 1, nx1), j = 1, ny1)
   close(20)

   ! read lon and lat for vector fields (like sea ice and ocean veclocities)
   open(24, file = inpath // 'grid.dat.pop')
   read(24,'(10f8.2)') ((ulat(i, j), i = 1, nx), j = 1, ny)
   read(24,'(10f8.2)') ((ulon(i, j), i = 1, nx), j = 1, ny)
   ! HTN, HTE, HUS, HUW are lengths of the 4 sides of a grid cell in km
   ! HTN, HTE are lengths of the northern and eastern sides of a scaler grid cell in km, HTN*HTE is
   ! the area of a scaler grid cell in km**2 and can be used to calculate sea ice volume and volumes
   ! of other variables
   ! HUS, HUW are lengths of the southern and western sides of a vector grid cell in km
   read(24,'(10f8.2)') ((HTN  (i, j), i = 1, nx), j = 1, ny)
   read(24,'(10f8.2)') ((HTE  (i, j), i = 1, nx), j = 1, ny)
   read(24,'(10f8.2)') ((HUS  (i, j), i = 1, nx), j = 1, ny)
   read(24,'(10f8.2)') ((HUW  (i, j), i = 1, nx), j = 1, ny)
   ! angle is the angle between latitude line and  grid cell x-coordinate line, needed for plotting
   ! vectors in spherical coordinate system
   read(24,'(10f8.2)') ((angle(i, j), i = 1, nx), j = 1, ny)
   close(24)

   ! read model grid mask; ocean levels > 0, land <= 0
   open(20, file = inpath // 'io.dat_360_120.output')
   read(20,'(360i2)') kmt
   close(20)

   area(:, :) = HTN(:, :) * HTE(:, :) * 1.0e6  ! m2

   ! --------------------------------------------------------------------------

   write(outfile, '(a,i4,a,i4,a)') outpath // 'OBS_PIOMAS_reanaly_2-1_T2Ms_sit_', &
      nyear1, '01-', nyear2, '12.nc'

   write(period, '(i4,a,i4)') nyear1, '-', nyear2

   ! open new NetCDF file (write)

   status = NF90_CREATE(trim(outfile), NF90_Write, ncID)
      if (status /= nf90_NoErr) call handle_err(status)

   ! set fill mode to NO_FILL

   status = NF90_SET_FILL(ncID, NF90_NOFILL, oldfillmode)
      if (status /= nf90_NoErr) call handle_err(status)

   ! create dimensions

   status = NF90_DEF_DIM(ncID, "i", nx1, iDimID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_DIM(ncID, "j", ny1, jDimID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_DIM(ncID, "time", NF90_UNLIMITED, TimeDimID)
      if (status /= nf90_NoErr) call handle_err(status)

   ! create dimension variables

   status = NF90_DEF_VAR(ncID, "i", NF90_INT, iDimID, iVarID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_VAR(ncID, "j", NF90_INT, jDimID, jVarID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_VAR(ncID, "time", NF90_DOUBLE, TimeDimID, TimeVarID)
      if (status /= nf90_NoErr) call handle_err(status)

   ! set attributes of dimension variables

   status = NF90_PUT_ATT(ncID, iVarID, "long_name", "cell index along second dimension")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, iVarID, "units", "1")
      if (status /= nf90_NoErr) call handle_err(status)

   status = NF90_PUT_ATT(ncID, jVarID, "long_name", "cell index along first dimension")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, jVarID, "units", "1")
      if (status /= nf90_NoErr) call handle_err(status)

   status = NF90_PUT_ATT(ncID, TimeVarID, "standard_name", "time")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, TimeVarID, "units", "months since 1950-01-01 00:00:00")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, TimeVarID, "axis", "T")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, TimeVarID, "long_name", "time")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, TimeVarID, "calendar", "standard")
      if (status /= nf90_NoErr) call handle_err(status)

   ! define variables

   status = NF90_DEF_VAR(ncID, "lat", NF90_FLOAT, (/iDimID, jDimID/), LatVarID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_VAR(ncID, "lon", NF90_FLOAT, (/iDimID, jDimID/), LonVarID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_VAR(ncID, "sit", NF90_FLOAT, (/iDimID, jDimID, TimeDimID/), outVarID)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_DEF_VAR(ncID, "areacello", NF90_FLOAT, (/iDimID, jDimID/), outAreaID)
      if (status /= nf90_NoErr) call handle_err(status)

   ! set attributes of variables

   status = NF90_PUT_ATT(ncID, LonVarID, "standard_name", "longitude")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, LonVarID, "long_name", "longitude")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, LonVarID, "units", "degrees_east")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, LonVarID, "_CoordinateAxisType", "Lon")
      if (status /= nf90_NoErr) call handle_err(status)

   status = NF90_PUT_ATT(ncID, LatVarID, "standard_name", "latitude")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, LatVarID, "long_name", "latitude")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, LatVarID, "units", "degrees_north")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, LatVarID, "_CoordinateAxisType", "Lat")
      if (status /= nf90_NoErr) call handle_err(status)

   status = NF90_PUT_ATT(ncID, outVarID, "standard_name", "sea_ice_thickness")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outVarID, "units", "m")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outVarID, "long_name", "Sea Ice Thickness")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outVarID, "comment", &
      "the mean thickness of sea ice in the ocean portion of the grid cell " &
      // "(averaging over the entire ocean portion, including the ice-free fraction). " &
      // "Reported as 0.0 in regions free of sea ice.")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outVarID, "cell_methods", "time: mean area: mean where sea")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outVarID, "cell_measures", "area: areacello")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outVarID, "_FillValue", fillval)
      if (status /= nf90_NoErr) call handle_err(status)

   status = NF90_PUT_ATT(ncID, outAreaID, "standard_name", "cell_area")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outAreaID, "units", "m2")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outAreaID, "long_name", "Ocean Grid-Cell Area")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, outAreaID, "_FillValue", fillval)
      if (status /= nf90_NoErr) call handle_err(status)

   ! set global attributes

   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "conventions", "CF/CMOR")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "title", "PIOMAS Arctic Sea Ice Volume Reanalysis")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "reference", "Zhang, Jinlun and D.A. Rothrock: " &
      // "Modeling global sea ice with a thickness and enthalpy distribution model in " &
      // "generalized curvilinear coordinates, Mon. Wea. Rev. 131(5), 681-697, 2003.")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "source", "http://psc.apl.uw.edu/research/" &
      // "projects/arctic-sea-ice-volume-anomaly/data/model_grid")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "tier", 2)
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "field", "T2Ms")
      if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "period", period)
      if (status /= nf90_NoErr) call handle_err(status)
!   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "user", "xxxx")
!      if (status /= nf90_NoErr) call handle_err(status)
!   status = NF90_PUT_ATT(ncID, NF90_GLOBAL, "history", "Created on xxxx")
!      if (status /= nf90_NoErr) call handle_err(status)

   ! end define mode

   status = NF90_ENDDEF(ncID)
   if (status /= nf90_NoErr) call handle_err(status)

   ! write variables i, j, lat, lon, areacello

   do i = 1, imt
      iVar(i) = i - 1
   end do

   do j = 1, jmt
      jVar(j) = j - 1
   end do

   status = NF90_PUT_VAR(ncID, iVarID, iVar)
       if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_VAR(ncID, jVarID, jVar)
       if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_VAR(ncID, lonVarID, clon)
       if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_VAR(ncID, latVarID, clat)
       if (status /= nf90_NoErr) call handle_err(status)
   status = NF90_PUT_VAR(ncID, outAreaID, area)
       if (status /= nf90_NoErr) call handle_err(status)

   ! --------------------------------------------------------------------------

   t = 1  ! time step counter

   do iyear = nyear1, nyear2

      write(unit = cyear(iyear), fmt = '(i4)') iyear
      i = slen(f2)
      open(2, file = f2(1:i) // cyear(iyear), &
         access = 'direct', form = 'unformatted', recl = nx1 * ny1 * 4, &
         status='unknown')

      do imon = 1, 12

         time_ref = (1950 * 12) + 1  ! reference month = 195001
         time     = (iyear * 12) + imon - time_ref ! = months since 195001

         read(2, rec = imon)((heff(i, j), i = 1, nx1), j = 1, ny1)

         where (kmt <= 0)
            heff = fillval
         endwhere

         status = NF90_PUT_VAR(ncID, outvarID, heff(:,:), start = (/1, 1, t/), &
            count = (/nx1, ny1, 1/))
         if (status /= nf90_NoErr) call handle_err(status)

         status = NF90_PUT_VAR(ncID, timeVarID, time, start=(/t/), count = (/1/))
         if (status /= nf90_NoErr) call handle_err(status)

         t = t + 1

      end do

      close(2)
   end do

   status = NF90_CLOSE(ncID)
   if (status /= nf90_NoErr) call handle_err(status)

   stop

end program convert_piomas

! -----------------------------------------------------------------------------

INTEGER FUNCTION slen (string)
   ! ---
   ! --- this function computes the length of a character string less
   ! --- trailing blanks
   ! --- slen > 0, length of string less trailing blanks
   ! ---      = 0, character string is blank
   ! ---
   CHARACTER*(*) string
   CHARACTER*1 cblank
   INTEGER i
   DATA cblank/' '/
   ! ---
   DO 50 i = LEN(string), 1, -1
      IF (string(i:i) .NE. ' ') GO TO 100
50 CONTINUE
   i = 0
100 CONTINUE
   slen = i
   RETURN
END

! -----------------------------------------------------------------------------

subroutine handle_err(status)

   use netcdf

   implicit none

   integer status

   print *, "NetCDF error: ", trim(NF90_STRERROR(status))
   print *, "Abort."
   stop

end subroutine handle_err

