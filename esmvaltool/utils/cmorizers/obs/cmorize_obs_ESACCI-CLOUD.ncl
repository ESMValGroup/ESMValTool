; #############################################################################
; ESMValTool CMORizer for ESACCI-CLOUD data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/
;
; Last access
;    20190201
;
; Download and processing instructions
;    Download the data from:
;      L3C/AVHRR-PM/
;    To fill the gap 199409-199501, also download:
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199409-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199410-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199411-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199412-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1995/ \
;        199501-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;    Put all files under a single directory (no subdirectories with years).
;
; Modification history
;    20190201-A_righ_ma: adapted to v2.
;    20181116-A_laue_ax: changed to use CLOUD-CCI v3.0 data (AVHRR-PM), gaps
;                        (1994/09 - 1995/01) are filled with AVHRR-AM data
;    20180522-A_righ_ma: changed to use AVHRR-PM data.
;    20160519-A_laue_ax: written (reformat_obs_ESACCI-AEROSOL.ncl).
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-CLOUD.ncl"

  ; Source name
  OBSNAME = "ESACCI-CLOUD"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1982
  YEAR2 = 2016

  ; MIP
  MIP = "Amon"

  ; Selected variable (standard name)
  VARS = (/"clt", "cltStderr", "clivi", "clwvi"/)

  ; Name in the raw data
  NAME = (/"cfc", "cfc_unc", "iwp_allsky", "lwp_allsky"/)

  ; Conversion factor
  CONV = (/100., 1., 0.001, 0.001/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + "/cmor/tables/" + \
    (/"cmip5/Tables/CMIP5_Amon", \
      "custom/CMOR_cltStderr.dat", \
      "cmip5/Tables/CMIP5_Amon", \
      "cmip5/Tables/CMIP5_Amon"/)

end

begin

  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1, YEAR2

      syear = sprinti("%i", yy)
      do mm = 1, 12

        smonth = sprinti("%0.2i", mm)

        ; Read file
        fname = systemfunc("ls " + input_dir_path + syear + smonth + \
                           "-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-*-" + \
                           "fv3.0.nc")

        ; No files found
        if (all(ismissing(fname))) then
          continue
        end if

        ; Extract data
        f = addfile(fname, "r")
        xx = f->$NAME(vID)$

        ; Convert units
        xx = xx * CONV(vID)

        ; lwp is not a CMOR variable, derive as clwi = lwp + iwp
        if (VARS(vID).eq."clwvi") then

          ; Read 2nd variable containing iwp (variable "iwp_allsky")
          iwpname = "iwp_allsky"
          xx2 = f->$iwpname$

          ; Convert units
          xx2 = xx2 * 0.001

          ; add iwp to lwp to calculate cmor variable "clwvi"
          xx = xx + xx2
          delete(xx2)
        end if

        ; Assign to global array
        if (.not.isdefined("output")) then
          dims = dimsizes(xx)
          dims(0) = dimsizes(time)
          output = new(dims, float)
          output!0 = "time"
          output&time = time
          output!1 = "lat"
          output&lat = f->lat
          output!2 = "lon"
          output&lon = f->lon
        end if
        output(ind(toint(yy * 100 + mm).eq.date), :, :) = (/xx/)

        delete(fname)
        delete(f)

      end do
    end do

    ; Set fill value
    output = where(output.eq.-999, output@_FillValue, output)

    ; Format time coordinate
    ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
      error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(output&lat) .eq. -1) then
      output = output(:, ::-1, :)
    end if
    clat = lat_attrib(output&lat)
    delete(output&lat)
    output&lat = clat
    delete(clat)

    ; Format longitude coordinate
    output!2 = "lon"
    if (isMonotonic(output&lon) .eq. 0) then
      error_msg("f", DIAG_SCRIPT, "", "non-monotonic longitude coordinate")
    end if
    if (any(output&lon.lt.0.)) then
      output = lonFlip(output)
    end if
    clon = lon_attrib(output&lon)
    delete(output&lon)
    output&lon = clon
    delete(clon)

    ; Set variable attributes
    tmp = var_attrib(output, VARS(vID), CMOR_TABLE(vID))
    delete(output)
    output = tmp
    delete(tmp)

    ; Set global attributes
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/" + \
      "CLD_PRODUCTS/v3.0/", \
      "Stengel et al., Earth Syst. Sci. Data, " + \
      "doi:10.5194/essd-9-881-2017, 2017", "")

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_sat_AVHRR-fv3.0_" + \
      MIP + "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vID), output, gAtt)
    delete(gAtt)
    delete(output)

  end do

end
