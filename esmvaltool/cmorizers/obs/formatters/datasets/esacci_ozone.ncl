; #############################################################################
; ESMValTool CMORizer for ESACCI-OZONE data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/ozone/data/
;
; Last access
;    20190201
;
; Download and processing instructions
;    Download the data from:
;      total_columns/l3/merged/v0100/
;    Put all files under a single directory (no subdirectories with years).
;
; Modification history
;    20190201-righi_mattia: adapted to v2 and replace NaN/inf with FillValue.
;    20160224-wenzel_sabrina: written based on reformat_obs_ESACCI-AEROSOL.ncl.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_ozone.ncl"

  ; Source name
  OBSNAME = "ESACCI-OZONE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = (/1997, 1997, 2007, 2007/)
  YEAR2 = (/2010, 2010, 2008, 2008/)

  ; Selected variable (standard name)
  VAR = (/"toz", "tozStderr", "tro3prof", "tro3profStderr"/)

  ; Name in the raw data
  NAME = (/"atmosphere_mole_content_of_ozone", \
           "atmosphere_mole_content_of_ozone_standard_error", \
           "merged_ozone_vmr", \
           "uncertainty_of_merged_ozone"/)

  ; MIP
  MIP = (/"Amon", "Amon", "Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_" + VAR + ".dat"

  ; File name
  FNAME = (/"ESACCI-OZONE-L3S-TC-MERGED-DLR_1M-_DATE_??-fv0100.nc", \
            "ESACCI-OZONE-L3S-TC-MERGED-DLR_1M-_DATE_??-fv0100.nc", \
            "ESACCI-OZONE-L3-LP-MERGED-MZM-_DATE_-fv0002.nc", \
            "ESACCI-OZONE-L3-LP-MERGED-MZM-_DATE_-fv0002.nc"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3"

  ; Global attributes
  SOURCE = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/ozone/data/"
  REF = "Loyola et al., Int. J. Remote Sens. doi:10.1080/" + \
    "01431160902825016, 2009."
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    time = create_timec(YEAR1(vv), YEAR2(vv))
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1(vv), YEAR2(vv)
      do mm = 1, 12

        ldate = yy + sprinti("%0.2i", mm)

        ; File name
        fname = systemfunc("ls " + input_dir_path + \
                           str_sub_str(FNAME(vv), "_DATE_", ldate))

        ; Check
        if (all(ismissing(fname))) then
          error_msg("f", DIAG_SCRIPT, "", "no file found for date " + ldate)
        end if

        ; Extract data
        f = addfile(fname(0), "r")
        xx = f->$NAME(vv)$
        xx@_FillValue = FILL
        xx@missing_value = xx@_FillValue
        xx = where(xx.lt.0., xx@_FillValue, xx)
        xx = where(xx.gt.1e35, xx@_FillValue, xx)  ; get rid of infinity values
        replace_ieeenan(xx, xx@_FillValue, 0)

        ; Assign to global array
        dimnames = getvardimnames(xx)
        if (.not.isdefined("output")) then
          dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
          output = new(dims, typeof(xx))
          output!0 = "time"
          output&time = time
          do ii = 0, dimsizes(dimnames) - 1
            if (dimnames(ii).eq."air_pressure") then
              output!(ii+1) = "plev"
              output&plev = f->$dimnames(ii)$
            elseif (isStrSubset(dimnames(ii), "latitude")) then
              output!(ii+1) = "lat"
              output&lat = f->$dimnames(ii)$
            elseif (dimnames(ii).eq."longitude")
              output!(ii+1) = "lon"
              output&lon = f->$dimnames(ii)$
            end if
          end do
        end if
        output(ind(toint(ldate).eq.date), :, :) = (/xx/)
        delete(fname)
        delete(xx)

      end do
    end do

    ; Format coordinates
    format_coords(output, YEAR1(vv) + "0101", YEAR2(vv) + "1231", FREQ(vv))

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
    DATESTR = YEAR1(vv) + "01-" + YEAR2(vv) + "12"
    fout = output_dir_path + \
      str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                 MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

    ; Write variable
    write_nc(fout, VAR(vv), output, bounds, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)
    delete(time)
    delete(date)

  end do

end
