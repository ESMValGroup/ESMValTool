; #############################################################################
; ESMValTool CMORizer for PATMOS-x data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://www.ncdc.noaa.gov/cdr/atmospheric/avhrr-cloud-properties-patmos-x
;
; Last access
;    20190210
;
; Download and processing instructions
;    Click on Download and download all the NOAA data, excluding the
;    preliminary, e.g. with:
;      wget -r --accept '*NOAA*.nc' --reject '*preliminary*' <source>
;    Put all files in input_dir_path (no subdirectories with years).
;    Select only complete years for both ascending and descending orbit.
;
; Caveats
;    The data are processed by calculating the average of the ascending and the
;    descending orbit on each day. Multiple files are available for some days,
;    in this case the most recent version (NOAA-vv) is chosen.

; Modification history
;    20190208-righi_mattia: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "patmos_x.ncl"

  ; Source name
  OBSNAME = "PATMOS-x"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = get_year(start_year, 1982)
  YEAR2 = get_year(end_year, 2016)

    ; Selected variable (standard name)
  VAR = (/"clt"/)

  ; Name in the raw data
  NAME = (/"cloud_fraction"/)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "NOAA"

  ; Global attributes
  SOURCE = "https://www.ncdc.noaa.gov/cdr/atmospheric/avhrr-cloud-" + \
    "properties-patmos-x"
  REF = "Heidinger et al., Bull. Amer. Meteor. Soc., " + \
    "doi:10.1175/BAMS-D-12-00246.1, 2013."
  COMMENT = ""

end

begin

  ; Read coordinates
  files = systemfunc("ls " + input_dir_path + "patmosx_*" + YEAR1 + "*.nc")
  f = addfile(files(0), "r")
  tmp = f->latitude
  lat = tmp * tmp@scale_factor + tmp@add_offset
  nlat = dimsizes(lat)
  delete(tmp)
  tmp = f->longitude
  lon = tmp * tmp@scale_factor + tmp@add_offset
  nlon = dimsizes(lon)
  delete(tmp)
  delete(files)
  delete(f)

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yy = YEAR1, YEAR2

      ; Define output monthly-mean array
      output = new((/12, nlat, nlon/), float)
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      output&time = create_timec(yy, yy)
      output&lat = lat
      output&lon = lon

      do mm = 1, 12

        ; Number of days
        nd = days_in_month(yy, mm)

        ; Define local array
        output_temp = new((/nd, nlat, nlon/), float)

        ; Date string for this month
        yyyymm = yy + sprinti("%0.2i", mm)

        do dd = 1, nd

          ; Date string for this day
          yyyymmdd = yy + sprinti("%0.2i", mm) + sprinti("%0.2i", dd)

          ; Ascending orbit
          files_asc = systemfunc("ls " + input_dir_path + \
                                 "patmosx_v??r??_NOAA-??_asc_d" + \
                                 yyyymm + "??_c*.nc | grep asc_d" + yyyymmdd)
          if (.not.all(ismissing(files_asc))) then
            ; Read most recent file
            f = addfile(files_asc(dimsizes(files_asc) - 1), "r")
            tmp = f->$NAME(vv)$
            xasc = tmp * tmp@scale_factor + tmp@add_offset
            delete(tmp)
          end if
          delete(files_asc)

          ; Descending orbit
          files_des = systemfunc("ls " + input_dir_path + \
                                 "patmosx_v??r??_NOAA-??_des_d" + \
                                 yyyymm + "??_c*.nc | grep des_d" + yyyymmdd)
          if (.not.all(ismissing(files_des))) then
            ; Read most recent file
            f = addfile(files_des(dimsizes(files_des) - 1), "r")
            tmp = f->$NAME(vv)$
            xdes = tmp * tmp@scale_factor + tmp@add_offset
            delete(tmp)
          end if
          delete(files_des)

          ; Skip if no data defined (output_temp will stay missing)
          if (.not.isdefined("xasc") .and. .not.isdefined("xdes")) then
            continue
          end if

          if (.not.isdefined("xasc")) then
            output_temp(dd, :, :) = (/xdes/)
            delete(xdes)
            continue
          end if

          if (.not.isdefined("xdes")) then
            output_temp(dd, :, :) = (/xasc/)
            delete(xasc)
            continue
          end if

          ; Replace missing values in one orbit with valid values from the
          ; other orbit, to avoid propagating missing values while averaging
          xasc = where(ismissing(xasc), xdes, xasc)
          xdes = where(ismissing(xdes), xasc, xdes)

          output_temp(dd - 1, :, :) = 0.5 * (xasc + xdes)
          delete(xasc)
          delete(xdes)

        end do  ; day

        ; Monthly mean
        output(mm - 1, :, :) = dim_avg_n(output_temp, 0)
        delete(output_temp)

      end do  ; month

      if (VAR(vv).eq."clt") then
        output = 100. * output  ; [1] --> [%]
      end if

      ; Format coordinates
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      format_coords(output, yy + "0101", yy + "1231", FREQ(vv))

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
      DATESTR = yy + "01-" + yy + "12"
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

    end do  ; year

  end do  ; variable

end
