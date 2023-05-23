; #############################################################################
; ESMValTool CMORizer for DARDAR-CLOUD (level 2) data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset (registration required).
;
; Source
;    ftp.icare.univ-lille1.fr/SPACEBORNE/CLOUDSAT/DARDAR-CLOUD.v3.10/
;
; Last access
;    20230517
;
; Download and processing instructions:
;
;      1) Create an account at the ICARE Data and Services Center
;         (https://www.icare.univ-lille.fr/asd-content/register)
;      2) Download the DARDAR-CLOUD .hdf files from
;         ftp.icare.univ-lille1.fr (via ftp)
;      3) Put all files for a specific year under a single directory
;         named like the year (e.g. 2007), no subdirectories with
;         days etc.
;
; Description:
;
; This script bins 3-dim cloud ice water content from individual
; orbits (level 2 data) onto a 2x2 degree grid with
; 40 vertical levels ("alt40"). Individual orbits and days are then averaged
; to monthly means. The DARDAR-CLOUD data are converted from kg/m3 to kg/kg
; using the air density from the US standard atmosphere.
;
; Note: one year of data consists of approximately 2400 GB and 5000 files.
;       Processing one year takes about xxxx hours. It is therefore recommended
;       to process only one year at a time instead of multiple years. The
;       year(s) to process are specified when running "esmvaltool" with
;       the options "--start=year1" and "--end=year2", e.g. (w/o line break)
;
;       esmvaltool data format --config_file [CONFIG_FILE] --start=2015
;                  --end=2015 DARDAR-CLOUD
;
; Modification history
;    20230517-lauer_axel: written.
;
; #############################################################################

loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")

begin
  ; Script name (for logger)
  diag_script = "dardar-cloud.ncl"

  ; Source name
  obsname = "DARDAR-CLOUD"

  ; Tier
  tier = 3

  ; Period (range of complete years: 2007-2016)
  year1 = get_year(start_year, 2007)
  year2 = get_year(end_year, 2007)

  if ((year2 - year1) .gt. 0) then
    log_info("Note: processing of " + tostring(year2 - year1 + 1) + \
             " years requested (start_year = " + sprinti("%04d", year1) + \
             ", end_year = " + sprinti("%04d", year2) + ").")
    log_info("One year of data consists of approximately 2400 GB and 5000 " \
             + "files, processing may take up to several hours.")
    log_info("It is therefore recommended to process only one year at a time.")
    log_info("The year(s) to process are specified when running esmvaltool " \
             + "with the options --start=year1 and --end=year2, e.g.")
    log_info("esmvaltool data format --config_file [CONFIG_FILE] " + \
             "--start=2015 --end=2015 DARDAR-CLOUD")
  end if

  ; output variable (standard name)
  outvar = (/"cli"/)

  ; input variables
  var = (/"iwc"/)

  ; MIP
  mip = "Amon"

  ; Frequency
  freq = "mon"

  ; CMOR table

  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_" + mip/)

  ; Type
  type = "sat"

  ; Version
  version = "v3.10"

  ; Global attributes
  source = "ftp.icare.univ-lille1.fr"
  ref = "Delanoë et al., J. Geophys. Res., " \
        + "doi:10.1029/2009JD012346, 2010"
end

begin
  ; define 2x2 output grid

  nx = 180
  ny = 90
  dx = 360. / nx
  dy = 180. / ny

  ; alt40 from CMIP6 table "coordinates"
  gridhgt = (/240., 720., 1200., 1680., 2160., 2640., 3120., 3600., 4080.,   \
             4560., 5040., 5520., 6000., 6480., 6960., 7440., 7920., 8400., \
             8880., 9360., 9840., 10320., 10800., 11280., 11760., 12240.,   \
             12720., 13200., 13680., 14160., 14640., 15120., 15600., 16080.,\
             16560., 17040., 17520., 18000., 18480., 18960./)

  nz = dimsizes(gridhgt)
  nt = 12

  ; define output longitudes and latitudes

  gridlat = latGlobeFo(ny, "lat", "latitude", "degrees_north")
  gridlon = lonGlobeFo(nx, "lon", "longitude", "degrees_east")

  ; work arrays (NCL (v6.6.2) procedure bin_sum did not work properly with
  ; subsetting the output arrays)
  work = new((/ny, nx/), float)
  workpts = new((/ny, nx/), integer)
  work = 0.
  workpts = 0

  ; "in-cloud" (all points)
  grid = new((/nt, nz, ny, nx/), float)
  gridpts = new((/nt, nz, ny, nx/), integer)
  grid2d = new((/nt, 2, ny, nx/), float)
  gridpts2d = new((/nt, 2, ny, nx/), integer)

  ; "grid-box average"
  grid_avg = new((/nt, nz, ny, nx/), float)
  gridpts_avg = new((/nt, nz, ny, nx/), integer)
  grid2d_avg = new((/nt, 2, ny, nx/), float)
  gridpts2d_avg = new((/nt, 2, ny, nx/), integer)

  ; boundaries of vertical (height) bins

  gridhgt_bnds = new(nz + 1, float)
  gridhgt_bnds(0) = 0.
  do i = 1, nz - 1
    gridhgt_bnds(i) = (gridhgt(i - 1) + gridhgt(i)) / 2.
  end do
  gridhgt_bnds(nz) = gridhgt(nz - 1) + gridhgt(nz - 1) - gridhgt_bnds(nz - 1)

  ; air density of vertical bins (from US standard atmosphere) in kg/m3
  tdp = stdatmus_z2tdp(gridhgt)
  air_density = tdp(1, :)
  delete(tdp)

  ; define coordinates

  grid!0 = "time"
  grid!1 = "alt40"
  grid&alt40 = gridhgt
  grid!2 = "lat"
  grid&lat = gridlat
  grid!3 = "lon"
  grid&lon = gridlon

  grid2d!0 = "time"
  grid2d!2 = "lat"
  grid2d&lat = gridlat
  grid2d!3 = "lon"
  grid2d&lon = gridlon

  copy_VarMeta(grid, grid_avg)
  copy_VarMeta(grid2d, grid2d_avg)

end

begin
  do calyear = year1, year2
    ; reset work arrays
    grid = 0.
    grid_avg = 0.

    grid2d = 0.
    grid2d_avg = 0.

    gridpts = 0
    gridpts_avg = 0

    gridpts2d = 0
    gridpts2d_avg = 0

    ; create time variable
    gridtime = create_timec(calyear, calyear)

    grid&time = gridtime
    grid_avg&time = gridtime

    grid2d&time = gridtime
    grid2d_avg&time = gridtime

    inpath = input_dir_path + tostring(calyear) + "/"

    ; process all days of the year

    do dayofyear = 1, 31 ;;day_of_year(calyear, 12, 31)
      print("day = " + tostring(dayofyear))
      caldate = yyyyddd_to_yyyymmdd(calyear * 10000 + dayofyear)
      calmonth = mod(caldate / 100, 100)
      calday = mod(caldate, 100)

      itime = calmonth - 1

      fn = tostring(calyear) + sprinti("%0.3i", dayofyear)
      fs = systemfunc("ls " + inpath + "DARDAR-CLOUD_" + fn + "*.nc")

      if (ismissing(fs(0))) then
        log_info("Warning: no input file(s) for " + tostring(calyear) + \
                 "-" + sprinti("%03d", dayofyear) + ", skipping day.")
        delete(fs)
        continue
      end if

      norbits = dimsizes(fs)

      ; start loop over all orbits for one day

      do orbit = 0, norbits - 1
        begTime = get_cpu_time()
        print("orbit = " + tostring(orbit + 1))
        ssplit = str_split(fs(orbit), "/")
        strs = str_split(ssplit(dimsizes(ssplit) - 1), "_")
        sorbit = strs(1)
        delete(ssplit)
        delete(strs)

        ; ***************************************************************
        ; Read input
        ; ***************************************************************
        print("processing " + fs(orbit) + "...")

        f = addfile(fs(orbit), "r")
        x2d = f->$var(0)$

        if (isatt(x2d, "_FillValue")) then
          delete(x2d@_FillValue)
        end if

        lat   = f->latitude
        lon   = f->longitude
        hgt = f->height
;        time  = f->time

        ; convert longitudes from -180...180 to 0...360
        lon = where(lon .lt. 0., lon + 360., lon)

;        ; convert height from km to m
;        hgt = hgt * 1000.0

        ; --------------------------------------------------------------------

        dimx = dimsizes(x2d)
        x = ndtooned(x2d)
        delete(x2d)

        lon_ext = ndtooned(conform_dims(dimx, lon, 0))
        lat_ext = ndtooned(conform_dims(dimx, lat, 0))
        hgt_ext = ndtooned(conform_dims(dimx, hgt, 1))

        ; =========================
        ; 3-dim cloud ice content
        ; =========================

        do ilev = 0, nz - 1
          ; Find all elements that are within a given height level
          ; (output grid) that contain valid (x > 0) or missing (x = 0) values;
          ; invalid values (x < 0) are filtered out.
          ii = ind((hgt_ext .ge. gridhgt_bnds(ilev)) .and. \
                   (hgt_ext .lt. gridhgt_bnds(ilev + 1)) .and. (x .ge. 0.))
          ; check if there are input data on the given output level and
          ; output gridbox
          if (.not. ismissing(ii(0))) then
            xselected = x(ii)
            xlon = lon_ext(ii)
            xlat = lat_ext(ii)
            ; *** gridbox average ***
            ; Includes all points with missing cloud liquid water
            ; content (x = 0) when calculating the average. This assumes
            ; that all grid points with "missing" but not invalid values
            ; are cloud-free. Invalid values (x < 0) are filtered out.
            bin_sum(work, workpts, gridlon, gridlat, xlon, xlat, xselected)
            grid_avg(itime, ilev, :, :) = grid_avg(itime, ilev, :, :) + work
            gridpts_avg(itime, ilev, :, :) = gridpts_avg(itime, ilev, :, :) \
              + workpts
            work = 0.
            workpts = 0
            ; *** in-cloud average ***
            ; Filter invalid (x < 0) and missing values (x = 0).
            idx = ind(xselected .gt. 0.)
            if (.not. ismissing(idx(0))) then
              bin_sum(work, workpts, gridlon, gridlat, \
                      xlon(idx), xlat(idx), xselected(idx))
              grid(itime, ilev, :, :) = grid(itime, ilev, :, :) + work
              gridpts(itime, ilev, :, :) = gridpts(itime, ilev, :, :) \
                + workpts
              work = 0.
              workpts = 0
            end if
            delete(idx)
            delete(xselected)
            delete(xlon)
            delete(xlat)
          end if
          delete(ii)
        end do  ; loop over output levels

        delete(x)
        delete(lat_ext)
        delete(lon_ext)
        delete(hgt_ext)

        delete(lat)
        delete(lon)
        delete(hgt)
;        delete(time)

        print("CPU time for orbit: " + tostring(get_cpu_time() - begTime) \
              + " s")
      end do  ; loop over all orbits per day
      delete(fs)
    end do  ; loop over days

    ; =========================
    ; 3-dim cloud water content
    ; =========================

    ; calculate monthly averages

    gridpts = where(gridpts .eq. 0, gridpts@_FillValue, gridpts)
    gridpts_avg = where(gridpts_avg .eq. 0, gridpts_avg@_FillValue, \
                        gridpts_avg)

    grid = grid / gridpts
    grid_avg = grid_avg / gridpts_avg

    ; convert from kg/m3 to kg/kg using air density from
    ; US standard atmosphere

    grid = grid / conform(grid, air_density, 1)
    grid_avg = grid_avg / conform(grid_avg, air_density, 1)

    ; ------------------------------------------------------------------------

    ; save results to files

    ; calculate coordinate bounds
    bounds = guess_coord_bounds(grid, freq)

    ; *** in-cloud ***

    format_coords(grid, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(grid, outvar(0), CMOR_TABLE(0))
;    delete(grid)
    ; Set global attributes
    comment = "Cloud ice water content (in-cloud)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-in-cloud", \
                 mip, outvar(0), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(0), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; *** gridbox average ***

    format_coords(grid_avg, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(grid_avg, outvar(0), CMOR_TABLE(0))
    ; Set global attributes
    comment = "Cloud ice water content (gridbox average)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-gridbox-average", \
                 mip, outvar(0), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(0), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

  end do  ; loop over years

end
