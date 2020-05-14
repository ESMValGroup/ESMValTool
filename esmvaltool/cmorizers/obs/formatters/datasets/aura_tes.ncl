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
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "aura_tes.ncl"

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

    if (.not.fileexists(fname)) then  ; alternative name
      fname = input_dir_path + "TES-Aura_L3-O3-M" + \
        yy + "m" + sprinti("%0.2i", mm) + "_C01_F01_10.he5"
    end if

    if (.not.fileexists(fname)) then
      log_info("File missing, skipping this date")
      continue
    end if

    ; Convert to NetCDF
    destination = output_dir_path + bname(fname) + ".nc"
    system(str_join((/"ncks", fname, destination/), " "))

    ; Read variable
    f = addfile(destination, "r")
    var = f->O3

    ; Read and assign coordinates
    var!0 = "plev"
    var&plev = todouble(f->Pressure)
    var!1 = "lat"
    var&lat = todouble(f->Latitude)
    var!2 = "lon"
    var&lon = todouble(f->Longitude)

    ; Create array
    if (.not.isdefined("output")) then
      outdim = array_append_record(dimsizes(timec), dimsizes(var), 0)
      output = new(outdim, float)
      output@_FillValue = var@_FillValue
    end if

    output(tt, :, :, :) = var
    delete(var)
    delete(f)
    system("rm -f " + destination)

  end do

  ; Assign standard FillValue
  output@_FillValue = FILL

  ; Assign coordinates
  output!0 = "time"
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
