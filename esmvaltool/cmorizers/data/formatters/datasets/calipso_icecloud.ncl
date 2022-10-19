; #############################################################################
; ESMValTool CMORizer for CALIPSO-Lidar Level 3 Ice Cloud Data, Standard V1-00
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset (registration required).
;
; Source
;    EarthData via https://eosweb.larc.nasa.gov/project/CALIPSO/
;      CAL_LID_L3_Ice_Cloud-Standard_V1-00
;
; Last access
;    20220804
;
; Download and processing instructions
; (requires EarthData login; see https://urs.earthdata.nasa.gov/)
;    1) Go to https://eosweb.larc.nasa.gov/project/CALIPSO/
;       CAL_LID_L3_Ice_Cloud-Standard_V1-00
;    2) Click on "Get Dataset"
;    3) Select Granules "A" containing both, day and night (Day/Night filter
;       = "Both")
;    4) Download selected granules
;    5) Enter EarthData login and password (https://urs.earthdata.nasa.gov/)
;    6) Follow download instructions in email from EarthData and put all
;       files in the same directory
;
; Modification history
;    20210331-lauer_axel: updated download information
;    20200414-lauer_axel: written
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "calipso_icecloud.ncl"

  ; Source name
  OBSNAME = "CALIPSO-ICECLOUD"

  ; Tier
  TIER = 3

  ; Period (complete years only)
  YEAR1 = 2007
  YEAR2 = 2015

  ; Selected variable (standard name)
  VAR = (/"cli"/)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_" + MIP/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "1-00"

  ; Global attributes
  SOURCE = "https://eosweb.larc.nasa.gov/project/CALIPSO/" + \
    "CAL_LID_L3_Ice_Cloud-Standard_V1-00"
  REF = "doi:10.5067/CALIOP/CALIPSO/L3_ICE_CLOUD-STANDARD-V1-00"
  COMMENT = "Atmospheric Science Data Center (ASDC)"

end

begin

  R_gas = 8.3145       ; gas constant (J/mol/K)
  mw_air = 28.9647e-3  ; molecular weight of dry air (kg/mol)
  R_air = R_gas / mw_air

  ; CMIP6 "plev27" pressure levels (Pa)

  plev27 = (/100000., 97500., 95000., 92500., 90000., 87500., 85000., \
            82500., 80000., 77500., 75000., 70000., 65000., 60000., \
            55000., 50000., 45000., 40000., 35000., 30000., 25000., \
            22500., 20000., 17500., 15000., 12500., 10000./)

  vv = 0
  log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

  ; calculate level bounds

  plev27_bnds = new(dimsizes(plev27) + 1, float)
  do i = 1, dimsizes(plev27) - 1
    plev27_bnds(i) = 0.5 * (plev27(i - 1) + plev27(i))
  end do
  plev27_bnds(0) = plev27(0) + 0.5 * (plev27(0) - plev27(1))
  n = dimsizes(plev27)
  plev27_bnds(n) = plev27(n - 1) - 0.5 * (plev27(n - 2) - plev27(n - 1))

  ; initialize time

  time = create_timec(YEAR1, YEAR2)
  date = cd_calendar(time, 1)

  ; Create timeseries
  do yy = YEAR1, YEAR2

    syear = sprinti("%i", yy)
    do mm = 1, 12

      smonth = sprinti("%0.2i", mm)

      print(syear + "-" + smonth)

      ; Read data

      ; There are 3 granules available:
      ;   day:   CAL_LID_L3_Ice_Cloud-Standard-V1-00.yyyy-mmD.hdf
      ;   night: CAL_LID_L3_Ice_Cloud-Standard-V1-00.yyyy-mmN.hdf
      ;   both:  CAL_LID_L3_Ice_Cloud-Standard-V1-00.yyyy-mmA.hdf
      ; Here, we only use the granule "A" (= day + night).

      fname = input_dir_path + "CAL_LID_L3_Ice_Cloud-Standard-V1-00." + \
        syear + "-" + smonth + "A.hdf"

      ; No input file found
      if (ismissing(fname)) then
        error_msg("f", DIAG_SCRIPT, "main", "input file for " + syear + \
                  "-" + smonth + " not found")
      end if

      ; Extract data
      f = addfile(fname, "r")

      iwc_hist = f->Ice_Water_Content_Histogram
      iwc_bnds = f->Ice_Water_Content_Bin_Boundaries
      cldfree_samp = f->Cloud_Free_Samples
      cld_samp = f->Cloud_Samples
      pres = f->Pressure_Mean * 100.0      ; hPa --> Pa
      temp = f->Temperature_Mean + 273.15  ; deg C --> K

      ; Create output array
      if (.not.isdefined("output")) then
;        alt = round(f->Altitude_Midpoint * 1000.0, 0)  ; km, round to 1 m
        lat = f->Latitude_Midpoint
        lon = f->Longitude_Midpoint
        ndims = new(4, integer)
        ndims(0) = dimsizes(time)
        ndims(1) = dimsizes(plev27)
        ndims(2) = dimsizes(lat)
        ndims(3) = dimsizes(lon)
        output = new(ndims, float)
        output!0 = "time"
        output&time = time
        output!1 = "plev"
        output&plev = plev27
        output!2 = "lat"
        output&lat = lat
        output!3 = "lon"
        output&lon = lon

        ; dimensions = lat, lon, plev
        iwc_avg_plev27 = new((/ndims(2), ndims(3), ndims(1)/), float)
        iwc_avg_plev27!0 = "lat"
        iwc_avg_plev27!1 = "lon"
        iwc_avg_plev27!2 = "plev"
      end if

      ; Calculate grid-box average ice water content from ice water content
      ; histogram. Formula from https://www-calipso.larc.nasa.gov/resources/
      ;                calipso_users_guide/qs/cal_lid_l3_ice_cloud_v1-00.php
      ; (section "In-cloud IWC and Grid-averaged IWC").
      ;
      ; Note: 1) ice water content is derived from extiction coefficient. The
      ;       extinction coefficient can be negative, thus negative ice water
      ;       contents are produced. These are included in the level 3 product
      ;       used here to "accurately represent the retrieved population in
      ;       the level 2 data product" (that was used to derive the level 3
      ;       product used here).
      ;       However, these unphysical ice water contents are excluded here
      ;       when calculating the grid box average ice water content for
      ;       comparison with earth system models.
      ;       2) ice water content is only calculated for clouds with
      ;       1.0e-5 <= iwc <= 1.0 g/m3.

      ; find first bin with iwc >= 1e-5 g/m3 (first bin with valid data)

      tmp = ind(iwc_bnds(:, 0) .ge. 1.0e-5)  ; min bound
      if (all(ismissing(tmp))) then
        error_msg("f", DIAG_SCRIPT, "main", "invalid bounds for histogram")
      end if
      i0 = tmp(0)
      delete(tmp)

      ; find last bin with iwc <= 1.0 g/m3 (last bin with valid data)
      tmp = ind(iwc_bnds(:, 2) .le. 1.0)  ; max bound
      if (all(ismissing(tmp))) then
        error_msg("f", DIAG_SCRIPT, "main", "invalid bounds for histogram")
      end if
      i1 = tmp(dimsizes(tmp) - 1)
      delete(tmp)

      iwc = iwc_hist * conform_dims(dimsizes(iwc_hist), iwc_bnds(:, 1), 3)
      samples = tofloat(cldfree_samp + cld_samp)
      samples@_FillValue = -999.
      samples = where(samples.gt.0., samples, samples@_FillValue)
      iwc_avg = dim_sum(iwc(:, :, :, i0:i1)) / samples

      ; *** convert units from g/m3 to kg/kg

      ; 1. calculate density of air (kg/m3)
      rho_air = pres / (R_air * temp)
      ; 2. convert from g/m3 to kg/kg
      iwc_avg = 1.0e-3 * iwc_avg / rho_air

;      fname = "debug.nc"
;      system("rm " + fname)
;      debugfile = addfile(fname, "c")
;      debugfile->rho = rho_air

      ; *** aggregate CALIPSO levels to plev27 vertical levels

      pres1d = ndtooned(pres)  ; convert to 1D array

      do ilev = 0, dimsizes(plev27) - 1
        ; all levels not within plev27_bnds(ilev) - plev27_bnds(ilev+1) are
        ; masked as missing values, which are then ignored when calculating
        ; averages with function dim_avg
        mapind = ind((pres1d .gt. plev27_bnds(ilev)) .or. \
                     (pres1d .lt. plev27_bnds(ilev + 1)))
        iwc1d = ndtooned(iwc_avg)
        iwc1d@_FillValue = -999.
        if (.not.all(ismissing(mapind))) then
          iwc1d(mapind) = iwc1d@_FillValue
        end if
        delete(mapind)
        iwc_select = onedtond(iwc1d, dimsizes(iwc_avg))
        delete(iwc1d)
        ; the layers of the CALIPSO data quite thin, so we do not bother
        ; calculating weights for each vertical level but simply
        ; calculate arithmetic averages when mapping the CALIPSO data to
        ; the plev27 vertical levels
        ; note: dim_avg ignores missing values
        iwc_avg_plev27(:, :, ilev) = dim_avg(iwc_select)
        delete(iwc_select)
      end do

      delete(pres1d)

      ; *** map gridbox average ice water content to output array

      output(ind(toint(yy * 100 + mm).eq.date), :, :, :) = \
        (/iwc_avg_plev27(plev|:, lat|:, lon|:)/)

      ; *** clean up

      delete(fname)
      delete(f)
      delete(iwc_avg)
      delete(iwc)
      delete(samples)
      delete(pres)
      delete(temp)
      delete(rho_air)
      delete(cldfree_samp)
      delete(cld_samp)
    end do
  end do

  ; Set fill value
  output = where(output.eq.-999, output@_FillValue, output)

  ; Format coordinates
  format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ(vv))

  ; Set variable attributes
  tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
  delete(output)
  output = tmp
  delete(tmp)

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(output, FREQ(vv))

  ; Set global attributes
  gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

  ; Output file
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE, VERSION, \
               MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

  ; Write variable
  write_nc(fout, VAR(vv), output, bounds, gAtt)
  delete(gAtt)
  delete(output)
  delete(bounds)

end
