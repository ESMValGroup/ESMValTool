; #############################################################################
; ESMValTool CMORizer for CloudSat Level 2 data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset (registration required).
;
; Source
;    https://www.cloudsat.cira.colostate.edu/
;
; Last access
;    20220810
;
; Download and processing instructions:
;
;    1) Create an account at the CloudSat Data Processing Center
;       (https://www.cloudsat.cira.colostate.edu/)
;
;    2) Download the datasets '2B-CWC-RO' and '2C-PRECIP-COLUMN' from
;       www.cloudsat.cira.colostate.edu (via sftp) following the instructions
;       given here: https://www.cloudsat.cira.colostate.edu/order/sftp-access
;
;    3) Put all files for a specific year under a single directory
;       named like the year (e.g. 2007), no subdirectories with
;       days etc.
;
; Description:
;
; This script bins 3-dim cloud liquid water content from individual CloudSat
; orbits (level 2 data: 2B_CWC_RO "radar only") onto a 2x2 degree grid with
; 40 vertical levels ("alt40"). Individual orbits and days are then averaged
; to monthly means. Precipitating gridboxes are identified using the CloudSat
; level 2 data for column precipitation (2C_PRECIP_COLUMN). The CloudSat data
; are converted from kg/m3 to kg/kg using the air density from the US
; standard atmosphere. In addition, liquid and ice water path are binned onto
; the same spatial grid. The script calculates four types of datasets for the
; 3-dim cloud liquid water content, liquid water path, ice water path:
;
;   1) in-cloud values averaged over all available data points
;   2) in-cloud values averaged over all data points for which no
;      precipitation has been detected
;   3) grid-box averages assuming that missing values are cloud-free, i.e.
;      cloud liquid water content = 0
;   4) grid-box averages calculated over all data points for which no
;      precipitation has been detected
;
; The results for each variable are written to four different output files:
; 1) in-cloud-all, 2) in-cloud-noprecip, 3) gridbox-average,
; 4) gridbox-average-noprecip.
;
; Note: one year of data consists of approximately 110 GB and 10,000 files.
;       Processing one year takes about 5-6 hours. It is therefore recommended
;       to process only one year at a time instead of multiple years. The
;       year(s) to process are specified when running "esmvaltool" with
;       the options "--start=year1" and "--end=year2", e.g. (w/o line break)
;
;       esmvaltool data format --config_file [CONFIG_FILE] --start=2015
;                  --end=2015 CLOUDSAT-L2
;
; Modification history
;    20230904-lauer_axel: added output of clwvi (iwp + lwp)
;    20220809-lauer_axel: adapted CMORizer to new format introduced in
;                         ESMValTool v2.5.0 and added info message
;    20210924-lauer_axel: added processing of lwp and iwp
;    20201119-lauer_axel: allowing more flexible filenames for precipitation
;                         data
;    20201117-lauer_axel: code completely rewritten to reduce cpu time per
;                         orbit from ~ 120 s to ~ 4 s
;    20201109-lauer_axel: written.
;
; #############################################################################

loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")

begin
  ; Script name (for logger)
  diag_script = "cloudsat-l2.ncl"

  ; Source name
  obsname = "CLOUDSAT-L2"

  ; Tier
  tier = 3

  ; Period (range of complete years: 2007-2017)
  year1 = get_year(start_year, 2007)
  year2 = get_year(end_year, 2007)

  if ((year2 - year1) .gt. 0) then
    log_info("Note: processing of " + tostring(year2 - year1 + 1) + \
             " years requested (start_year = " + sprinti("%04d", year1) + \
             ", end_year = " + sprinti("%04d", year2) + ").")
    log_info("One year of data consists of approximately 110 GB and 10,000 " \
             + "files, processing may take up to several hours.")
    log_info("It is therefore recommended to process only one year at a time.")
    log_info("The year(s) to process are specified when running esmvaltool " \
             + "with the options --start=year1 and --end=year2, e.g.")
    log_info("esmvaltool data format --config_file [CONFIG_FILE] " + \
             "--start=2015 --end=2015 CLOUDSAT-L2")
  end if

  ; output variable (standard name)
  outvar = (/"clw", "lwp", "clivi", "clwvi"/)

  ; input variables
  var = (/"LO_RO_liquid_water_content_2B_CWC_RO", \
          "LO_RO_liquid_water_path_2B_CWC_RO",    \
          "IO_RO_ice_water_path_2B_CWC_RO",       \
          "dummy"/)

  var_flag = "Precip_flag_2C_PRECIP_COLUMN"

  ; MIP
  mip = "Amon"

  ; Frequency
  freq = "mon"

  ; CMOR table

  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_" + mip, \
      "/custom/CMOR_lwp.dat", \
      "/cmip5/Tables/CMIP5_" + mip, \
      "/cmip5/Tables/CMIP5_" + mip/)

  ; Type
  type = "sat"

  ; Version
  version = "P1-R05"

  ; Global attributes
  source = "https://www.cloudsat.cira.colostate.edu/"
  ref = "Stephens et al., Bull. Am. Meteorol. Soc., " \
        + "doi:10.1175/BAMS-D-16-0324.1, 2017"
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

  grid2d = new((/nt, 3, ny, nx/), float)
  gridpts2d = new((/nt, 3, ny, nx/), integer)

  ; "grid-box average" (all points)
  grid_avg = new((/nt, nz, ny, nx/), float)
  gridpts_avg = new((/nt, nz, ny, nx/), integer)
  grid2d_avg = new((/nt, 3, ny, nx/), float)
  gridpts2d_avg = new((/nt, 3, ny, nx/), integer)

  ; "in-cloud" (no precipitation)
  grid_noprecip = new((/nt, nz, ny, nx/), float)
  gridpts_noprecip = new((/nt, nz, ny, nx/), integer)
  grid2d_noprecip = new((/nt, 3, ny, nx/), float)
  gridpts2d_noprecip = new((/nt, 3, ny, nx/), integer)

  ; "grid-box average" (no precipitation)
  grid_avg_noprecip = new((/nt, nz, ny, nx/), float)
  gridpts_avg_noprecip = new((/nt, nz, ny, nx/), integer)
  grid2d_avg_noprecip = new((/nt, 3, ny, nx/), float)
  gridpts2d_avg_noprecip = new((/nt, 3, ny, nx/), integer)

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
  copy_VarMeta(grid, grid_noprecip)
  copy_VarMeta(grid, grid_avg_noprecip)

  copy_VarMeta(grid2d, grid2d_avg)
  copy_VarMeta(grid2d, grid2d_noprecip)
  copy_VarMeta(grid2d, grid2d_avg_noprecip)

end

begin
  do calyear = year1, year2
    ; reset work arrays
    grid = 0.
    grid_avg = 0.
    grid_noprecip = 0.
    grid_avg_noprecip = 0.

    grid2d = 0.
    grid2d_avg = 0.
    grid2d_noprecip = 0.
    grid2d_avg_noprecip = 0.

    gridpts = 0
    gridpts_avg = 0
    gridpts_noprecip = 0
    gridpts_avg_noprecip = 0

    gridpts2d = 0
    gridpts2d_avg = 0
    gridpts2d_noprecip = 0
    gridpts2d_avg_noprecip = 0

    ; create time variable
    gridtime = create_timec(calyear, calyear)

    grid&time = gridtime
    grid_avg&time = gridtime
    grid_noprecip&time = gridtime
    grid_avg_noprecip&time = gridtime

    grid2d&time = gridtime
    grid2d_avg&time = gridtime
    grid2d_noprecip&time = gridtime
    grid2d_avg_noprecip&time = gridtime

    inpath = input_dir_path + tostring(calyear) + "/"

    ; process all days of the year

    do dayofyear = 1, day_of_year(calyear, 12, 31)
      print("day = " + tostring(dayofyear))
      caldate = yyyyddd_to_yyyymmdd(calyear * 10000 + dayofyear)
      calmonth = mod(caldate / 100, 100)
      calday = mod(caldate, 100)

      itime = calmonth - 1

      fn = tostring(calyear) + sprinti("%0.3i", dayofyear)
      fs = systemfunc("ls " + inpath + fn + "*2B-CWC-RO_GRANULE*.hdf")

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
        ; combine filename for precipitation data
        ssplit = str_split(fs(orbit), "/")
        strs = str_split(ssplit(dimsizes(ssplit) - 1), "_")
        sorbit = strs(1)
        fs_precip = systemfunc("ls " + inpath + strs(0) + "_" + sorbit + \
                               "_CS_2C-PRECIP-COLUMN_GRANULE_P1_R05_*.hdf")
        delete(ssplit)
        delete(strs)

        ; File extension must be .hdfeos for NCL to correctly read input
        ; files. File extension of cloudsat data downloaded from cloudsat
        ; data processing center is .hdf.
        ; --> workaround: set symbolic links with extension .hdfeos before
        ;     reading input files

        ssplit = str_split(fs(orbit), "/")
        slink1 = ssplit(dimsizes(ssplit) - 1) + "eos"
        delete(ssplit)

        ssplit = str_split(fs_precip, "/")
        slink2 = ssplit(dimsizes(ssplit) - 1) + "eos"
        delete(ssplit)

        ; remove old (possibly broken) links before creating new ones

        if (fileexists(slink1)) then
          system("rm -f " + slink1)
        end if
        if (fileexists(slink2)) then
          system("rm -f " + slink2)
        end if

        ; create symbolic links

        system("ln -s " + fs(orbit) + " " + slink1)
        system("ln -s " + fs_precip + " " + slink2)

        ; ***************************************************************
        ; Read hdf-eos
        ; ***************************************************************
        print("processing " + fs(orbit) + " + " + fs_precip + "...")

        f = addfile(slink1, "r")
        x2d = f->$var(0)$

        if (isatt(x2d, "_FillValue")) then
          delete(x2d@_FillValue)
        end if

        l1d = f->$var(1)$

        if (isatt(l1d, "_FillValue")) then
          delete(l1d@_FillValue)
        end if

        i1d = f->$var(2)$

        if (isatt(i1d, "_FillValue")) then
          delete(i1d@_FillValue)
        end if

        f_precip = addfile(slink2, "r")
        x_flag = f_precip->$var_flag$

        ; remove attribute '_FillValue' to prevent NCL from handling
        ; missing values

        if (isatt(x_flag, "_FillValue")) then
          delete(x_flag@_FillValue)
        end if

        ; check if array dimensions match

        dimx = dimsizes(x2d)
        if (dimx(0) .ne. dimsizes(x_flag)) then
          log_info("Warning: dimsizes do not match. " + fs(orbit) + ", " + \
                   fs_precip + " Skipping orbit " + sorbit + ".")
          continue
        end if

        lat   = f->Latitude_2B_CWC_RO
        lon   = f->Longitude_2B_CWC_RO
        hgt1d = ndtooned(f->Height_2B_CWC_RO)
;        time  = f->Profile_time_2B_CWC_RO

        ; convert longitudes from -180...180 to 0...360
        lon = where(lon .lt. 0., lon + 360., lon)

        ; --------------------------------------------------------------------

        x = ndtooned(x2d)
        delete(x2d)

        x_flag_ext = ndtooned(conform_dims(dimx, x_flag, 0))

        lon_ext = ndtooned(conform_dims(dimx, lon, 0))
        lat_ext = ndtooned(conform_dims(dimx, lat, 0))

        ; =========================
        ; 3-dim cloud water content
        ; =========================

        do ilev = 0, nz - 1
          ; Find all elements that are within a given height level
          ; (output grid) that contain valid (x > 0) or missing (x = 0) values;
          ; invalid values (x < 0) are filtered out.
          ii = ind((hgt1d .ge. gridhgt_bnds(ilev)) .and. \
                   (hgt1d .lt. gridhgt_bnds(ilev + 1)) .and. (x .ge. 0.))
          ; check if there are input data on the given output level and
          ; output gridbox
          if (.not. ismissing(ii(0))) then
            xselected = x(ii)
            xlon = lon_ext(ii)
            xlat = lat_ext(ii)
            x_flag_selected = x_flag_ext(ii)
            ; *** gridbox average, all ***
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
            ; *** in-cloud average, all ***
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
            ; *** gridbox average, no precipitation ***
            ; Filter cases that have been flagged with
            ; precipitation or as unknown. The following CloudSat flags
            ; are defined (see https://www.cloudsat.cira.colostate.edu/
            ; sites/default/files/
            ; products/files/2C-PRECIP-COLUMN_PDICD.P2_R04.20130124.pdf):
            ;
            ; flag 0 = no precipitation detected
            ;      1 = rain possible
            ;      2 = rain probable
            ;      3 = rain certain
            ;      4 = snow possible
            ;      5 = snow certain
            ;      6 = mixed precipitation possible
            ;      7 = mixed precipitation certain
            ;      9 = uncertain
            idx = ind(x_flag_selected .eq. 0)
            if (.not. ismissing(idx(0))) then
              bin_sum(work, workpts, gridlon, gridlat, \
                      xlon(idx), xlat(idx), xselected(idx))
              grid_avg_noprecip(itime, ilev, :, :) = \
                grid_avg_noprecip(itime, ilev, :, :) + work
              gridpts_avg_noprecip(itime, ilev, :, :) = \
                gridpts_avg_noprecip(itime, ilev, :, :) + workpts
              work = 0.
              workpts = 0
            end if
            delete(idx)
            ; *** in-cloud average, no precipitation ***
            ; Filter cases that have been flagged with
            ; precipitation or as unknown.
            idx = ind((xselected .gt. 0.) .and. (x_flag_selected .eq. 0))
            if (.not. ismissing(idx(0))) then
              bin_sum(work, workpts, gridlon, gridlat, \
                      xlon(idx), xlat(idx), xselected(idx))
              grid_noprecip(itime, ilev, :, :) = \
                grid_noprecip(itime, ilev, :, :) + work
              gridpts_noprecip(itime, ilev, :, :) = \
                gridpts_noprecip(itime, ilev, :, :) + workpts
              work = 0.
              workpts = 0
            end if
            delete(idx)

            delete(xselected)
            delete(x_flag_selected)
            delete(xlon)
            delete(xlat)
          end if
          delete(ii)
        end do  ; loop over output levels

        delete(x)
        delete(x_flag_ext)
        delete(lat_ext)
        delete(lon_ext)
        delete(hgt1d)

        ; ============================================================
        ; 2-dim cloud liquid / ice water path / total cloud water path
        ; ============================================================

        do ivar = 1, 3
          if (outvar(ivar) .eq. "lwp") then
            x = l1d
            ilev = 0
          else if (outvar(ivar) .eq. "clivi") then
            x = i1d
            ilev = 1
          else if (outvar(ivar) .eq. "clwvi") then
            ; clwvi is calculated from the *output* fields of lwp and iwp
            ; to make sure it is the sum of the two (masking, etc.)
            ; --> nothing to do here
            continue
          else
            log_info("Warning: output variable unknown: " + outvar(ivar) + \
                     ", skipping variable " + outvar(ivar) + ".")
            continue
          end if
          end if
          end if

          ; Find all elements that contain valid (x > 0) or missing (x = 0)
          ; values; invalid values (x < 0) are filtered out.
          ii = ind(x .ge. 0.)
          ; check if there are input data for given output gridbox
          if (.not. ismissing(ii(0))) then
            xselected = x(ii)
            xlon = lon(ii)
            xlat = lat(ii)
            x_flag_selected = x_flag(ii)
            ; *** gridbox average, all ***
            ; Includes all points with missing lwp
            ; content (x = 0) when calculating the average. This assumes
            ; that all grid points with "missing" but not invalid values
            ; are cloud-free. Invalid values (x < 0) are filtered out.
            bin_sum(work, workpts, gridlon, gridlat, xlon, xlat, xselected)
            grid2d_avg(itime, ilev, :, :) = \
              grid2d_avg(itime, ilev, :, :) + work
            gridpts2d_avg(itime, ilev, :, :) = \
              gridpts2d_avg(itime, ilev, :, :) + workpts
            work = 0.
            workpts = 0
            ; *** in-cloud average, all ***
            ; Filter invalid (x < 0) and missing values (x = 0).
            idx = ind(xselected .gt. 0.)
            if (.not. ismissing(idx(0))) then
              bin_sum(work, workpts, gridlon, gridlat, \
                      xlon(idx), xlat(idx), xselected(idx))
              grid2d(itime, ilev, :, :) = grid2d(itime, ilev, :, :) + work
              gridpts2d(itime, ilev, :, :) = gridpts2d(itime, ilev, :, :) \
                + workpts
              work = 0.
              workpts = 0
            end if
            delete(idx)
            ; *** gridbox average, no precipitation ***
            ; Filter cases that have been flagged with
            ; precipitation or as unknown.
            idx = ind(x_flag_selected .eq. 0)
            if (.not. ismissing(idx(0))) then
              bin_sum(work, workpts, gridlon, gridlat, \
                      xlon(idx), xlat(idx), xselected(idx))
              grid2d_avg_noprecip(itime, ilev, :, :) = \
                grid2d_avg_noprecip(itime, ilev, :, :) + work
              gridpts2d_avg_noprecip(itime, ilev, :, :) = \
                gridpts2d_avg_noprecip(itime, ilev, :, :) + workpts
              work = 0.
              workpts = 0
            end if
            delete(idx)
            ; *** in-cloud average, no precipitation ***
            ; Filter cases that have been flagged with
            ; precipitation or as unknown.
            idx = ind((xselected .gt. 0.) .and. (x_flag_selected .eq. 0))
            if (.not. ismissing(idx(0))) then
              bin_sum(work, workpts, gridlon, gridlat, \
                      xlon(idx), xlat(idx), xselected(idx))
              grid2d_noprecip(itime, ilev, :, :) = \
                grid2d_noprecip(itime, ilev, :, :) + work
              gridpts2d_noprecip(itime, ilev, :, :) = \
                gridpts2d_noprecip(itime, ilev, :, :) + workpts
              work = 0.
              workpts = 0
            end if
            delete(idx)

            delete(xselected)
            delete(x_flag_selected)
            delete(xlon)
            delete(xlat)
          end if
          delete(ii)
        end do  ; loop over implemented 2-dim variables (lwp, iwp)

        delete(x)
        delete(x_flag)
        delete(lat)
        delete(lon)
        delete(l1d)
        delete(i1d)
;        delete(time)

        ; remove temporary symbolic links

        system("rm -f " + slink1)
        system("rm -f " + slink2)

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
    gridpts_noprecip = where(gridpts_noprecip .eq. 0, \
                             gridpts_noprecip@_FillValue, gridpts_noprecip)
    gridpts_avg_noprecip = where(gridpts_avg_noprecip .eq. 0, \
                                 gridpts_avg_noprecip@_FillValue, \
                                 gridpts_avg_noprecip)

    grid = grid / gridpts
    grid_avg = grid_avg / gridpts_avg
    grid_noprecip = grid_noprecip / gridpts_noprecip
    grid_avg_noprecip = grid_avg_noprecip / gridpts_avg_noprecip

    ; convert from mg/m3 to kg/kg using air density from
    ; US standard atmosphere

    grid = grid / conform(grid, air_density, 1) * 1.e-6
    grid_avg = grid_avg / conform(grid_avg, air_density, 1) * 1.e-6
    grid_noprecip = grid_noprecip / \
      conform(grid_noprecip, air_density, 1) * 1.e-6
    grid_avg_noprecip = grid_avg_noprecip / \
      conform(grid_avg_noprecip, air_density, 1) * 1.e-6

    ; ------------------------------------------------------------------------

    ; save results to files

    ; calculate coordinate bounds
    bounds = guess_coord_bounds(grid, freq)

    ; *** in-cloud, all ***

    format_coords(grid, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(grid, outvar(0), CMOR_TABLE(0))
;    delete(grid)
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content (in-cloud, all)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-in-cloud-all", \
                 mip, outvar(0), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(0), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; *** gridbox average, all ***

    format_coords(grid_avg, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(grid_avg, outvar(0), CMOR_TABLE(0))
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content " \
      + "(gridbox average, all)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-gridbox-average-all", \
                 mip, outvar(0), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(0), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; *** in-cloud, no precipitation only ***

    format_coords(grid_noprecip, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(grid_noprecip, outvar(0), CMOR_TABLE(0))
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content " \
      + "(in-cloud, no precipitation)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-in-cloud-noprecip", \
                 mip, outvar(0), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(0), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; *** gridbox average, no precipitation only ***

    format_coords(grid_avg_noprecip, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(grid_avg_noprecip, outvar(0), CMOR_TABLE(0))
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content " \
      + "(gridbox average, no precipitation)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-gridbox-average-noprecip", \
                 mip, outvar(0), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(0), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; ===================================
    ; 2-dim cloud liquid / ice water path
    ; ===================================

    ; calculate monthly averages

    gridpts2d = where(gridpts2d .eq. 0, gridpts2d@_FillValue, gridpts2d)
    gridpts2d_avg = where(gridpts2d_avg .eq. 0, gridpts2d_avg@_FillValue, \
                          gridpts2d_avg)
    gridpts2d_noprecip = where(gridpts2d_noprecip .eq. 0, \
                               gridpts2d_noprecip@_FillValue, \
                               gridpts2d_noprecip)
    gridpts2d_avg_noprecip = where(gridpts2d_avg_noprecip .eq. 0, \
                                   gridpts2d_avg_noprecip@_FillValue, \
                                   gridpts2d_avg_noprecip)

    grid2d = grid2d / gridpts2d
    grid2d_avg = grid2d_avg / gridpts2d_avg
    grid2d_noprecip = grid2d_noprecip / gridpts2d_noprecip
    grid2d_avg_noprecip = grid2d_avg_noprecip / gridpts2d_avg_noprecip

    ; convert from g/m2 to kg/m2

    grid2d = grid2d * 1.e-3
    grid2d_avg = grid2d_avg * 1.e-3
    grid2d_noprecip = grid2d_noprecip * 1.e-3
    grid2d_avg_noprecip = grid2d_avg_noprecip * 1.e-3

    ; ------------------------------------------------------------------------

    ; save results to files

    do ivar = 1, 3
      if (outvar(ivar) .eq. "lwp") then
        ilev = 0
      else if (outvar(ivar) .eq. "clivi") then
        ilev = 1
      else if (outvar(ivar) .eq. "clwvi") then
        ; calculate clwvi (lwp+iwp) as sum of output fields for lwp and iwp
        ilev = 2
        grid2d(:, ilev, :, :) = grid2d(:, 0, :, :) \
          + grid2d(:, 1, :, :)
        grid2d_avg(:, ilev, :, :) = grid2d_avg(:, 0, :, :) \
          + grid2d_avg(:, 1, :, :)
        grid2d_noprecip(:, ilev, :, :) = grid2d_noprecip(:, 0, :, :) \
          + grid2d_noprecip(:, 1, :, :)
        grid2d_avg_noprecip(:, ilev, :, :) = grid2d_avg_noprecip(:, 0, :, :) \
          + grid2d_avg_noprecip(:, 1, :, :)
      end if
      end if
      end if

      out_grid = grid2d(:, ilev, :, :)
      out_grid_avg = grid2d_avg(:, ilev, :, :)
      out_grid_noprecip = grid2d_noprecip(:, ilev, :, :)
      out_grid_avg_noprecip = grid2d_avg_noprecip(:, ilev, :, :)

      ; calculate coordinate bounds
      bounds = guess_coord_bounds(out_grid, freq)

      ; *** in-cloud, all ***

      format_coords(out_grid, calyear + "0101", calyear + "1231", freq)
      ; Set variable attributes
      output = format_variable(out_grid, outvar(ivar), CMOR_TABLE(ivar))
      delete(out_grid)
      ; Set global attributes
      comment = "Liquid-only Radar-only Liquid Water Content (in-cloud, all)"
      gAtt = set_global_atts(obsname, tier, source, ref, comment)
      ; output filename
      datestr = calyear + "01-" + calyear + "12"
      fout = output_dir_path + \
        str_join((/"OBS", obsname, type, version + "-in-cloud-all", \
                   mip, outvar(ivar), datestr/), "_") + ".nc"
      ; write variable
      write_nc(fout, outvar(ivar), output, bounds, gAtt)
      delete(output)
      delete(gAtt)

    ; *** gridbox average, all ***

    format_coords(out_grid_avg, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(out_grid_avg, outvar(ivar), CMOR_TABLE(ivar))
    delete(out_grid_avg)
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content " \
      + "(gridbox average, all)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-gridbox-average-all", \
                 mip, outvar(ivar), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(ivar), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; *** in-cloud, no precipitation only ***

    format_coords(out_grid_noprecip, calyear + "0101", calyear + "1231", freq)
    ; Set variable attributes
    output = format_variable(out_grid_noprecip, outvar(ivar), CMOR_TABLE(ivar))
    delete(out_grid_noprecip)
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content " \
      + "(in-cloud, no precipitation)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-in-cloud-noprecip", \
                 mip, outvar(ivar), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(ivar), output, bounds, gAtt)
    delete(output)
    delete(gAtt)

    ; *** gridbox average, no precipitation only ***

    format_coords(out_grid_avg_noprecip, calyear + "0101", calyear + "1231", \
                  freq)
    ; Set variable attributes
    output = format_variable(out_grid_avg_noprecip, outvar(ivar), \
                             CMOR_TABLE(ivar))
    delete(out_grid_avg_noprecip)
    ; Set global attributes
    comment = "Liquid-only Radar-only Liquid Water Content " \
      + "(gridbox average, no precipitation)"
    gAtt = set_global_atts(obsname, tier, source, ref, comment)
    ; output filename
    datestr = calyear + "01-" + calyear + "12"
    fout = output_dir_path + \
      str_join((/"OBS", obsname, type, version + "-gridbox-average-noprecip", \
                 mip, outvar(ivar), datestr/), "_") + ".nc"
    ; write variable
    write_nc(fout, outvar(ivar), output, bounds, gAtt)
    delete(output)
    delete(gAtt)
    end do  ; loop over implemented 2-dim variables (lwp, iwp, lwp+iwp)
  end do  ; loop over years

end
