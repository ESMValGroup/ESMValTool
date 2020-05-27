; #############################################################################
;  ESMValTool CMORizer for HALOE data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    Grooss, J.-U. and Russell III, J. M., Atmos. Chem. Phys., 5, 2797-2807,
;    doi:10.5194/acp-5-2797-2005, 2005.
;
; Last access
;    20200311
;
; Download and processing instructions
;    Download and untar the supplementary material of the above reference.
;
; Caveats
;    The time period includes incomplete years. According to the given
;    reference (Section 3): "Between October 1991 and August 2002, a total of
;    78600 HALOE profiles of O3, CH4, H2O, HCl, and HF mixing ratio were then
;    combined into bins of 5 degree equivalent latitude for 22 pressure levels
;    between 316 and 0.1 hPa and 12 months".
;    In order to create T3M fields a fake longitude coordinate is added.
;
;  Modification history
;     20200311-righi_mattia: ported from v1.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "haloe.ncl"

  ; Source name
  OBSNAME = "HALOE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = get_year(start_year, 1991)
  YEAR2 = get_year(end_year, 2002)

  ; Selected variable (standard name)
  VAR = (/"hus", "tro3"/)

  ; Name in the raw data
  NAME = (/"H2O", "O3"/)

  ; MIP
  MIP = (/"Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon", "/cmip5/Tables/CMIP5_Amon"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE = "https://www.atmos-chem-phys.net/5/2797/2005/" + \
    "acp-5-2797-2005-supplement.tar"
  REF = " Grooss and Russell, Atmos. Chem. Phys., " + \
    "doi:10.5194/acp-5-2797-2005, 2005."
  COMMENT = ""

end

begin

  ; Open input file
  f = addfile(input_dir_path + "haloe_eqlat_monthly.nc", "r")

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    xx = f->$NAME(vv)$

    ; hus: convert mixing ratio to specific humidity
    if (VAR(vv).eq."hus") then
      xx = xx * 18.015 / 28.97  ; [m3/m3] --> [kg/kg]
      xx = xx / (1. + xx)
    end if

    ; tro3: convert mixing ratio to ppbv
    if (VAR(vv).eq."tro3") then
      xx = xx * 1.e9  ; [m3/m3] --> [1.e-9] (i.e., [nmol/mol])
    end if

    ; Define output array
    tmp = new((/dimsizes(time), dimsizes(xx&press), dimsizes(xx&eqlat)/), \
              float)
    tmp!0 = "time"
    tmp&time = time

    ; Extract indexes for valid time range (see caveats)
    idx1 = ind(cd_calendar(time, -1).eq.199110)
    idx2 = ind(cd_calendar(time, -1).eq.200208)
    tmp(idx1:idx2, :, :) = (/xx/)

    ; Create a fake longitude coordinate
    nlon = 50
    dims = array_append_record(dimsizes(tmp), nlon, 0)
    output = new(dims, float)
    output!0 = "time"
    output!1 = "plev"
    output!2 = "lat"
    output!3 = "lon"
    output&time = time
    output&plev = 100 * f->press  ; [hPa] --> [Pa]
    output&lat = f->eqlat
    output&lon = lonGlobeF(nlon, "lon", "latitude", "degrees_East")
    do ii = 0, dimsizes(output&lon) - 1
      output(:, :, :, ii) = (/tmp/)
    end do
    copy_VarAtts(xx, output)
    delete(tmp)
    delete(xx)

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

  end do

end
