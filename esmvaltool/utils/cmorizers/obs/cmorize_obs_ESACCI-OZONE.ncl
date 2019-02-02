; #############################################################################
; ESMValTool CMORizer for ESACCI-OZONE data
; #############################################################################
;
; Tier
;   Tier 2: other freely-available dataset.
;
; Source
;   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/ozone/data/
;
; Last access
;   20190201
;
; Download and processing instructions
;   Download the data from:
;     total_columns/l3/merged/v0100/
;   Put all files under a single directory (no subdirectories with years).
;
; Modification history
;   20190201-A_righ_ma: adapted to v2.
;   20160224-A_wenz_sa: written based on reformat_obs_ESACCI-AEROSOL.ncl.
;
;#############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Source name
  OBSNAME = "ESACCI-OZONE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1997
  YEAR2 = 2010

  ; Field
  FIELD = "T2Ms"

  ; Selected variable (standard name)
  VARS = (/"toz", "tozStderr"/)

  ; Name in the raw data
  NAME = (/"atmosphere_mole_content_of_ozone", \
           "atmosphere_mole_content_of_ozone_standard_error"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_" + VARS + ".dat"

end

begin

  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1, YEAR2
      do mm = 1, 12

        ldate = yy + sprinti("%0.2i", mm)

        ; File name
        fname = systemfunc("ls " + input_dir_path + "ESACCI-OZONE-L3S-TC-" + \
                           "MERGED-DLR_1M-" + ldate + "??-fv0100.nc")

        ; Check
        if (all(ismissing(fname))) then
          error_msg("f", diag_script, "", "no file found for date " + ldate)
        end if

        ; Extract data
        f = addfile(fname(0), "r")
        xx = f->$NAME(vID)$
        xx@_FillValue = FILL
        xx@missing_value = xx@_FillValue
        xx = where(xx.lt.0., xx@_FillValue, xx)

        ; Assign to global array
        if (.not.isdefined("output")) then
          dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
          output = new(dims, typeof(xx))
          output!0 = "time"
          output&time = time
          output!1 = "lat"
          output&lat = f->latitude
          output!2 = "lon"
          output&lon = f->longitude
        end if
        output(ind(toint(ldate).eq.date), :, :) = (/xx/)
        delete(fname)
        delete(xx)

      end do
    end do

    ; Format time coordinate
    ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
      error_msg("f", diag_script, "", "non-monotonic latitude coordinate")
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
      error_msg("f", diag_script, "", "non-monotonic longitude coordinate")
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
    gAtt = True
    gAtt@comment = "L3 total ozone columns merged using the " + \
      "algorithm GTO-ECV as part of the ESA Ozone CCI"
    gAtt@history = "Created on " + systemfunc("date")
    gAtt@host    = systemfunc("echo $HOST")
    gAtt@user    = systemfunc("echo $USER")
    gAtt@period  = YEAR1 + "-" + YEAR2
    gAtt@field   = FIELD
    gAtt@tier    = TIER
    gAtt@source  = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/ozone/data/"
    gAtt@reference = "Loyola et al., International " + \
      "Journal of Remote Sensing, vol. 30, no. 15, pp. 4295-4318, 2009."
    gAtt@title = OBSNAME + " data reformatted for the ESMValTool v2"
    gAtt@conventions = "CF/CMOR"

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_sat_L3_" + FIELD + \
      "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vID), output, gAtt)
    delete(gAtt)
    delete(output)

  end do

end
