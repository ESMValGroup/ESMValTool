; #############################################################################
; ESMValTool CMORizer for ESACCI-SOILMOISTURE data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/soil_moisture/data/
;
; Last access
;    20190201
;
; Download and processing instructions
;    Download the data from:
;      daily_files/COMBINED/v04.2/
;      ancillary/v04.2/
;    Put all files under a single directory (no subdirectories with years).
;
; Modification history
;    20190201-righi_mattia: adapted to v2, use new input data version 4.2.
;    20160824-lauer_axel: added processing of volumetric soil moisture
;                         content (sm, smStderr).
;    20160721-lauer_axel: use daily files, added processing of uncertainty.
;    20150523-righi_mattia: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_soilmoisture.ncl"

  ; Source name
  OBSNAME = "ESACCI-SOILMOISTURE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1979
  YEAR2 = 2016

  ; Selected variable (standard name)
  VAR = (/"sm", "smStderr", "dos", "dosStderr"/)

  ; Name in the raw data
  NAME = (/"sm", "sm_uncertainty", "sm", "sm_uncertainty"/)

  ; MIP
  MIP = (/"Lmon", "Lmon", "Lmon", "Lmon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_" + VAR + ".dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3S-SSMV-COMBINED-v4.2"

  ; Global attributes
  SOURCE = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/soil_moisture/data/"
  REF = \
    "Liu et al., Hydrol. Earth Syst. Sci., doi:10.5194/hess-15-425-2011, 2011."
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yy = YEAR1, YEAR2

      ; Set list of files
      files = systemfunc("ls " + input_dir_path + \
                         "ESACCI-SOILMOISTURE-L3S-SSMV-" + \
                         "COMBINED-" + yy + "????000000-fv04.2.nc")
      f = addfiles(files, "r")
      delete(files)

      ; Read data
      xx = f[:]->$NAME(vv)$
      if (isatt(xx, "scale_factor")) then
        tmp = tofloat(xx * xx@scale_factor)
        copy_VarAtts(xx, tmp)
        copy_VarCoords(xx, tmp)
        delete(xx)
        xx = tmp
        delete(tmp)
      end if
      delete(f)

      ; Derive dos using porosity
      if (any(VAR(vv).eq.(/"dos", "dosStderr"/))) then
        g = addfile(input_dir_path + \
                    "/ESACCI-SOILMOISTURE-POROSITY_V01.1.nc", "r")
        zz = g->porosity
        xx = xx * 100. / conform(xx, zz, (/1, 2/))
        delete(zz)
      end if

      ; Add a minor time shift for correct extraction of monthly mean below
      xx&time = xx&time + 0.1

      ; Calculate monthly means
      if (isStrSubset(VAR(vv), "Stderr")) then
        xx2 = xx
        xx2 = xx ^ 2  ; save metadata
        tmp = calculate_monthly_values(xx2, "avg", 0, False)
        delete(xx)
        delete(xx2)
        xx = sqrt(tmp)
        copy_VarAtts(tmp, xx)
        copy_VarCoords(tmp, xx)
        delete(tmp)
      else
        tmp = calculate_monthly_values(xx, "avg", 0, False)
        delete(xx)
        xx = tmp
        delete(tmp)
      end if

      ; Append to time-series
      if (.not.isdefined("output")) then
        output = xx
      else
        output := array_append_record(output, xx, 0)
      end if
      delete(xx)

    end do

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
