; #############################################################################
; ESMValTool CMORizer for ISCCP-FH data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://isccp.giss.nasa.gov/pub/flux-fh/tar-nc4_MPF/
;
; Last access
;    20191107
;
; Download and processing instructions
;    1) Go to https://isccp.giss.nasa.gov/projects/flux.html and click on
;       "FH On-line Data"
;    2) username and password are provided on this website
;    3) go to directory tar-nc4_MPF/ (= monthly files)
;    3) download and unpack files
;
; Modification history
;    20191108-lauer_axel: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")
begin

  ; Script name (for logger)
  DIAG_SCRIPT = "isccp_fh.ncl"

  ; Source name
  OBSNAME = "ISCCP-FH"

  ; Tier
  TIER = 2

  ; Period (only complete years)
  YEAR1 = 1984
  YEAR2 = 2016

  ; Selected variables (standard names)
  VAR = (/"rlut", "rlutcs", "rsut", "rsutcs", "rsus", \
          "rsds", "rlus", "rlds", \
          "rsdt", "ps", "ts", "tas", "prw", "alb"/)

  ; Name in the raw data
  NAME = (/"txu5fl", "tru5cr", "sxu5fl", "sru5cr", "sxu1fl", \
           "sxd1fl", "txu1fl", "txd1fl", \
           "sxd5fl", "ps____", "ts____", "ta____", "mnpwfl", "al_srf"/)

  ; MIP
  MIP = new(dimsizes(VAR), string)
  MIP = "Amon"

  ; Frequency
  FREQ = new(dimsizes(VAR), string)
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/custom/CMOR_alb.dat"/)

  ; Conversion factor
  CONV = (/1., 1., 1., 1., 1., 1., 1., 1., 1., 100., 1., 1., 10., 1./)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "v0"

  ; Global attributes
  SOURCE = "https://isccp.giss.nasa.gov/projects/flux.html"
  REF = "Zhang et al., Study of Cloud and Water Processes in Weather and " \
        + "Climate through Satellite Observations (submitted)"
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")
    time = create_timec(YEAR1, YEAR2)
    setfileoption("nc", "MissingToFillValue", False)

    do yr = YEAR1, YEAR2
      do mo = 1, 12

        fname = input_dir_path + "ISCCP-FH.MPF.v.0.0.GLOBAL." + yr + "." \
                + sprinti("%0.2i", mo) + ".nc"
        setfileoption("nc", "MissingToFillValue", False)
        f = addfile(fname, "r")
        input = f->$NAME(vv)$

        if (.not.isdefined("output")) then
          lat = (/f->latitude/)
          lon = (/f->longitude/)
          tmp = dimsizes(input)
          dims = (/dimsizes(time), tmp(0), tmp(1)/)
          output = new(dims, float)
          output!0 = "time"
          output!1 = "lat"
          output!2 = "lon"
          output&time = time
          output&lat = lat
          output&lon = lon
          delete(tmp)
          delete(dims)
        end if

        timeidx = (yr - YEAR1) * 12 + mo - 1

        ; apply scaling
        if (isatt(input, "scale_factor")) then
          output(timeidx, :, :) = tofloat(input) / tofloat(input@scale_factor)
        else
          output(timeidx, :, :) = tofloat(input)
        end if

        ; unit conversion
        output(timeidx, :, :) = output(timeidx, :, :) * CONV(vv)

        ; set fill value
        if (isatt(input, "missing_value")) then
          output(timeidx, :, :) = where(input.eq.input@missing_value, \
                                        output@_FillValue, \
                                        output(timeidx, :, :))
        end if

      end do
    end do

    delete(input)

    ; format time coordinate
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
