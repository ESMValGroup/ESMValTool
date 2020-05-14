; #############################################################################
; ESMValTool CMORizer for ESACCI-AEROSOL data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/aerosol/data/
;
; Last access
;    20190124
;
; Download and processing instructions
;    Download the data from:
;      ATSR2_SU/L3/v4.21/MONTHLY/ (1997-2002)
;      AATSR_SU/L3/v4.21/MONTHLY/ (2003-2011)
;    Other years are not considered since they are not complete.
;    Put all files in input_dir_path (no subdirectories with years).
;
; Modification history
;    20190124-righi_mattia: adapted to v2.
;    20160718-lauer_axel: added AOD550 + AOD870 uncertainties.
;    20160525-righi_mattia: updated to v4.21 and adding more variables.
;    20150126-righi_mattia: adding AOD at other wavelengths.
;    20151124-righi_mattia: switched to monthly raw data (now available).
;    20150424-righi_mattia: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_aerosol.ncl"

  ; Source name
  OBSNAME = "ESACCI-AEROSOL"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1997
  YEAR2 = 2011

  ; Selected variable (standard name)
  VAR = (/"od550aer", "od870aer", "od550lt1aer", "abs550aer", \
          "od550aerStderr", "od870aerStderr"/)

  ; Name in the raw data
  NAME = (/"AOD550_mean", "AOD870_mean", "FM_AOD550_mean", "AAOD550_mean", \
           "AOD550_uncertainty", "AOD870_uncertainty"/)

  ; MIP
  MIP = (/"aero", "aero", "aero", "aero", \
          "aero", "aero"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon", \
           "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_aero", \
      "/cmip5/Tables/CMIP5_aero", \
      "/cmip5/Tables/CMIP5_aero", \
      "/cmip5/Tables/CMIP5_aero", \
      "/custom/CMOR_od550aerStderr.dat", \
      "/custom/CMOR_od870aerStderr.dat"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "SU-v4.21"

  ; Global attributes
  SOURCE = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/aerosol/data/"
  REF = "Popp et al., Remote Sens., doi:10.3390/rs8050421, 2016."
  COMMENT = "Combined dataset ERS2-ATSR2 (1997-2002) and ENVISAT-AATSR " + \
    "(2003-2011), based on the University of Swansea algorithm " + \
    "(monthly mean L3 data)"

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1, YEAR2
      do mm = 1, 12

        ldate = yy + sprinti("%0.2i", mm)

        ; Read file
        fname = systemfunc("ls " + input_dir_path + ldate + "*.nc")

        ; No files found
        if (all(ismissing(fname))) then
          continue
        end if

        ; Extract data
        f = addfile(fname, "r")
        xx = f->$NAME(vv)$

        ; Assign to global array
        if (.not.isdefined("output")) then
          dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
          output = new(dims, float)
          output!0 = "time"
          output&time = time
          output!1 = "lat"
          output&lat = f->latitude
          output!2 = "lon"
          output&lon = f->longitude
        end if
        output(ind(toint(ldate).eq.date), :, :) = (/xx/)
        delete(fname)

      end do
    end do

    ; Set fill value
    output = where(output.eq.-999, output@_FillValue, output)

    ; Format coordinates
    output!0 = "time"
    output!1 = "lat"
    output!2 = "lon"
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

  end do

end
