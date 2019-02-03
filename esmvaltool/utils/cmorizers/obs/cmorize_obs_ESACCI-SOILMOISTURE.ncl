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
;    20190201-A_righ_ma: adapted to v2, use new input data version 4.2.
;    20160824-A_laue_ax: added processing of volumetric soil moisture
;                        content (sm, smStderr).
;    20160721-A_laue_ax: use daily files, added processing of uncertainty.
;    20150523-A_righ_ma: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-SOILMOISTURE.ncl"

  ; Source name
  OBSNAME = "ESACCI-SOILMOISTURE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1979
  YEAR2 = 2016

  ; Field
  FIELD = "T2Ms"

  ; Selected variable (standard name)
  VARS = (/"sm", "smStderr", "dos", "dosStderr"/)

  ; Name in the raw data
  NAME = (/"sm", "sm_uncertainty", "sm", "sm_uncertainty"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_" + VARS + ".dat"

end

begin

  do vv = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vv))

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
      if (any(VARS(vv).eq.(/"dos", "dosStderr"/))) then
        g = addfile(input_dir_path + \
                    "/ESACCI-SOILMOISTURE-POROSITY_V01.1.nc", "r")
        zz = g->porosity
        xx = xx * 100. / conform(xx, zz, (/1, 2/))
        delete(zz)
      end if

      ; Add a minor time shift for correct extraction of monthly mean below
      xx&time = xx&time + 0.1

      ; Calculate monthly means
      if (isStrSubset(VARS(vv), "Stderr")) then
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
    tmp = var_attrib(output, VARS(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    delete(tmp)

    ; Set global attributes
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/soil_moisture/data/", "", "")

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + \
      "_sat_L3S-SSMV-COMBINED-v4.2_" + FIELD + "_" + VARS(vv) + "_" + \
      YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write output
    write_nc(fout, VARS(vv), output, gAtt)
    delete(output)

  end do

end
