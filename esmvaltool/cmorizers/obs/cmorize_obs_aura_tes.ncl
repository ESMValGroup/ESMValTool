; #############################################################################
; ESMValTool CMORizer for AURA-TES data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://search.earthdata.nasa.gov/search?q=TL3O3M
;
; Last access
;    20181208
;
; Download and processing instructions
;    Select the V004 and V005 projects.
;    Download the script file for both projects.
;    To download the data use: wget -i <scriptfile>.
;    From the downloaded data, pick only the *.he5 files and put them in
;    input_dir_path.
;    Data is freely available, but a registration is required.
;
; Modification history
;    20190108-righi_mattia: adapted to v2.
;    20140129-righi_mattia: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_aura_tes.ncl"

  ; Source name
  OBSNAME = "AURA-TES"

  ; Tier
  TIER = 3

  ; Gridtype
  GRIDTYPE = "Nadir"

  ; Period
  YEAR1 = 2005
  YEAR2 = 2011

  ; Selected variable
  VAR = "tro3"

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE = "https://search.earthdata.nasa.gov/search?q=TL3O3M"
  REF =  "Beer, R., IEEE Trans. Geosci. Rem. Sens., " + \
    "doi:10.1109/TGRS.2005.863716, 2006"
  COMMENT = ""

end

begin

  ; Create time coordinate
  timec = create_timec(YEAR1, YEAR2)
  datec = ut_calendar(timec, 0)

  ; Loop over time
  do tt = 0, dimsizes(timec) - 1

    yy = toint(datec(tt, 0))
    mm = toint(datec(tt, 1))

    log_info("Processing date " + yy + sprinti("%0.2i", mm))

    ; Find files
    fname  = input_dir_path + "TES-Aura_L3-O3-M" + \
      yy + "m" + sprinti("%0.2i", mm) + "_F01_10.he5"

    if (.not.isfilepresent(fname)) then  ; alternative name
      fname = input_dir_path + "TES-Aura_L3-O3-M" + \
        yy + "m" + sprinti("%0.2i", mm) + "_C01_F01_10.he5"
    end if

    if (.not.isfilepresent(fname)) then
      log_info("File missing, skipping this date")
      continue
    end if

    ; Read variable
    f = addfile(fname, "r")
    str = "O3_" + GRIDTYPE + "Grid"
    var = f->$str$

    ; Reorder
    d1 = "XDim_" + GRIDTYPE + "Grid"
    d2 = "YDim_" + GRIDTYPE + "Grid"
    d3 = "nLevels_" + GRIDTYPE + "Grid"
    tmp = var($d1$|:, $d2$|:, $d3$|:)
    delete(var)
    var = tmp
    delete(tmp)

    ; Read vertical coordinate
    str = "Pressure_" + GRIDTYPE + "Grid"
    press = f->$str$
    var!2 = "plev"
    var&plev = todouble(press)

    ; Create array
    if (.not.isdefined("outvar")) then
      outdim = array_append_record(dimsizes(timec), dimsizes(var), 0)
      outvar = new(outdim, float)
      outvar@_FillValue = var@_FillValue
    end if

    outvar(tt, :, :, :) = var
    delete(var)

  end do

  ; Reorder
  outvar!0 = "time"
  outvar!1 = "lon"
  outvar!2 = "lat"
  outvar!3 = "plev"
  output = outvar(time|:, plev|:, lat|:, lon|:)
  output@_FillValue = FILL
  delete(outvar)

  ; Assign coordinates
  output&time = timec
  output&plev = output&plev * 100.  ; [hPa] --> [Pa]
  output&plev@units = "Pa"

  ; Convert units [mol/mol] --> [1e9]
  output = output * 1.e9
  output@units = "1e-9"

  ; Format coordinates
  format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ)

  ; Set variable attributes
  tmp = format_variable(output, VAR, CMOR_TABLE)
  delete(output)
  output = tmp
  delete(tmp)

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(output, FREQ)

  ; Set global attributes
  gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

  ; Output file
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE, VERSION, \
               MIP, VAR, DATESTR/), "_") + ".nc"

  ; Write variable
  write_nc(fout, VAR, output, bounds, gAtt)
  delete(gAtt)
  delete(output)

end
