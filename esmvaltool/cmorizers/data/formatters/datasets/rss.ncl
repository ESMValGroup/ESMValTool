; #############################################################################
; REFORMAT SCRIPT FOR RSS prw REANALYISIS DATA
; #############################################################################
;
; Tier
;    Tier 2: other freely available dataset.
;
; Source
;    ftp://ftp.remss.com/vapor/monthly_1deg/
;
; Last access
;    2020-11-19
;
; Download and processing instructions
;    A registration is required for downloading the data, but no licence
;    agreement necessary. Download vapor, here ncdf3 file is used.
;
; Caveats
;
; Modification history
;    20170419-A_gier_bettina: written.
;    20201201-A_weigel_katja: portet to ESMValTool v2.
;    20231214-A_weigel_katja: update to current ESMValTool version.
;
; #############################################################################

loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")
loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/utilities.ncl")

begin

    ; Source name
    OBSNAME = "RSS"

    DIAG_SCRIPT = "rss.ncl"

    ; Tier
    TIER = 2

    ; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
      TIER + "/" + OBSNAME + "/"

    ; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
      TIER + "/" + OBSNAME + "/"

    ; Period
    YEAR1 = 1988
    YEAR2 = 2019

    ; Selected variable (standard name)
    VAR = "prw"

    ; Name in the raw data
    NAME = "precipitable_water"
    NAMEANO = "precipitable_water_anomaly"
    CLIMNAME = "precipitable_water_climatology"

    ; Initialize global variable
    FIELD = "T2Ms"

    ; MIP
    MIP = "Amon"

    ; Frequency
    FREQ = "mon"

    ; CMOR table
    CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP
    ; + ".json"

    ; Type
    TYPE = "reanaly"

    ; Version
    VERSION = "v07r01"

    ; Global attributes
    SOURCE = "ftp://ftp.remss.com/vapor/monthly_1deg/"
    REF = "Remote Sensing Systems. " \
          + "Monthly Mean Total Precipitable Water Data Set " \
          + "on a 1 degree grid made from Remote Sensing " \
          + "Systems Version-7 Microwave Radiometer Data, " \
          + "Date created: 20201110T185551Z " \
          + "[accessed on 2020-11-19]. Santa Rosa, CA, USA. " \
          + "Available at www.remss.com"
    COMMENT = ""

end

begin
    ; verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    log_info("Processing " + VAR + " (" + MIP + ")")

    fname = input_dir_path + "tpw_" + VERSION + "_198801_202010.nc4.nc"

    f = addfile(fname, "r")
    output = (/f->$NAME$/)
    outputano = (/f->$NAMEANO$/)
    clim = (/f->$CLIMNAME$/)
    lat = (/f->latitude/)
    lon = (/f->longitude/)

    print(dimsizes(output))
    ; Cut off the first 10 months of 2020
    output_new = output(:dimsizes(output(:, 0, 0))-11, :, :)
    delete(output)
    output = output_new
    outputano_new = outputano(:dimsizes(outputano(:, 0, 0))-11, :, :)
    delete(outputano)
    outputano = outputano_new

    n = 0
    do iii = 0, dimsizes(outputano(:, 0, 0)) - 1, 12
        do nnn = 0, 11
            outputano_new(iii + nnn, :, :) = outputano(iii + nnn, :, :) + \
                                          clim(nnn, :, :)
        end do
    end do

    outputano_new = where(outputano.eq.(-500), outputano@_FillValue, outputano_new)
    delete(outputano)
    outputano = outputano_new

    ; Format coordinates
    output!0 = "time"
    output!1 = "lat"
    output!2 = "lon"
    output&lat = lat
    output&lon = lon
    outputano!0 = "time"
    outputano!1 = "lat"
    outputano!2 = "lon"
    outputano&lat = lat
    outputano&lon = lon
    ; Format latitude coordinate
    if (isMonotonic(output&lat) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(output&lat) .eq. -1) then
        output = output(:, ::-1, :)
    end if
    if (isMonotonic(outputano&lat) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(outputano&lat) .eq. -1) then
        outputano = outputano(:, ::-1, :)
    end if

    ; Format longitude coordinate
    if (isMonotonic(output&lon) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic longitude coordinate")
    end if
    if (any(output&lon.lt.0.)) then
        output = lonFlip(output)
    end if
    if (isMonotonic(outputano&lon) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic longitude coordinate")
    end if
    if (any(outputano&lon.lt.0.)) then
        outputano = lonFlip(outputano)
    end if

    output&time = create_timec(YEAR1, YEAR2)
    outputano&time = create_timec(YEAR1, YEAR2)

    format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ)
    format_coords(outputano, YEAR1 + "0101", YEAR2 + "1231", FREQ)

    ; Calculate days per month
    date = cd_calendar(output&time, 0)
    dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
    dpmc = conform(output, dpm, 0)
    ; Check time range
    if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
        error_msg("f", diag_script, "", "incorrect number of timesteps")
    end if
    
    date = cd_calendar(outputano&time, 0)
    dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
    dpmc = conform(outputano, dpm, 0)
    ; Check time range
    if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
        error_msg("f", diag_script, "", "incorrect number of timesteps")
    end if

    ; Values of -500 are ice -> set these to missing values
    output = where(output.eq.(-500), output@_FillValue, output)
    outputano = where(outputano.eq.(-500), outputano@_FillValue, outputano)

    ; Set variable attributes
    tmp = format_variable(output, VAR, CMOR_TABLE)
    delete(output)
    output = tmp
    delete(tmp)

    tmp = format_variable(outputano, VAR, CMOR_TABLE)
    delete(outputano)
    outputano = tmp
    delete(tmp)

    ; Set global attributes
    gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)
    gAtt@period = YEAR1 + "-" + YEAR2

    bounds = guess_coord_bounds(output, FREQ)
    boundsano = guess_coord_bounds(outputano, FREQ)

    filter = outputano
    filter = where(outputano.ne.(outputano@_FillValue), 1, filter)

    ; Outfile
    DATESTR = YEAR1 + "01-" + YEAR2 + "12"
    fout = output_dir_path + \
      str_join((/"OBS", OBSNAME, TYPE, VERSION, \
               MIP, VAR, DATESTR/), "_") + ".nc"

    filtout = output_dir_path + \
      str_join((/"OBS", OBSNAME, "filter", VERSION, \
               MIP, VAR, DATESTR/), "_") + ".nc"
    write_nc(fout, VAR, output, bounds, gAtt)
    write_nc(filtout, VAR, filter, boundsano, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)

end
