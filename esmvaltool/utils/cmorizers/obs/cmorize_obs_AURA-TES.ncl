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
;    From the downloaded data, pick only the *.he5 files and save them in the
;    input_dir_path.
;    Data is freely available, but a registration is required.
;
; Modification history
;    20190108-A_righ_ma: adapted for v2.
;    20140129-A_righ_ma: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_AURA-TES.ncl"

  ; Source name
  OBSNAME = "AURA-TES"

  ; Tier
  TIER = 3

  ; Gridtype
  GRIDTYPE = "Nadir"

  ; Period
  YEAR1 = 2005
  YEAR2 = 2011

  ; Field
  FIELD = "T3M"

  ; CMOR table
  CMOR_TABLE = \
    getenv("esmvaltool_root") + "/cmor/tables/cmip5/Tables/CMIP5_Amon"

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
    fname  = input_dir_path + "/" + \
      "TES-Aura_L3-O3-M" + yy + "m" + sprinti("%0.2i", mm) + \
      "_F01_10.he5"

    if (.not.isfilepresent(fname)) then  ; alternative name
      fname = input_dir_path + "/" + \
        "TES-Aura_L3-O3-M" + yy + "m" + sprinti("%0.2i", mm) + \
        "_C01_F01_10.he5"
    end if

    if (.not.isfilepresent(fname)) then
      error_msg("w", DIAG_SCRIPT, "", "file missing, skipping this date")
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

  ; Rename coordinates
  outvar!0 = "time"
  outvar!1 = "lon"
  outvar!2 = "lat"

  ; Convert units [mol/mol] --> [1e9]
  outvar = outvar * 1.e9
  outvar@units = "1e9"

  ; Reorder
  output = outvar(time|:, plev|:, lat|:, lon|:)
  output = output(:, ::-1, :, :)
  output@_FillValue = FILL
  output&time = timec
  delete(outvar)

  ; Set variable attributes
  tmp = var_attrib(output, "tro3", CMOR_TABLE)
  delete(output)
  output = tmp
  delete(tmp)

  ; Format time coordinate
  ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
  delete(output&time)
  output&time = ctime
  delete(ctime)

  ; Format vertical coordinate
  output!1 = "plev"
  output&plev = output&plev * 100.  ; [hPa] --> [Pa]
  output&plev@units = "Pa"
  output = output(:, ::-1, :, :)  ; bottom to top
  cplev = plev_attrib(output&plev)
  delete(output&plev)
  output&plev = cplev
  delete(cplev)

  ; Format latitude coordinate
  clat = lat_attrib(output&lat)
  delete(output&lat)
  output&lat = clat
  delete(clat)

  ; Format longitude coordinate
  output = lonFlip(output)
  clon = lon_attrib(output&lon)
  delete(output&lon)
  output&lon = clon
  delete(clon)

  ; Set global attributes
  gAtt = True
  gAtt@history       = "Created on " + systemfunc("date")
  gAtt@host          = systemfunc("echo $HOST")
  gAtt@user          = systemfunc("echo $USER")
  gAtt@period        = YEAR1 + "-" + YEAR2
  gAtt@field         = FIELD
  gAtt@tier          = TIER
  gAtt@source        = "https://search.earthdata.nasa.gov/"
  gAtt@reference     = "Beer, R., IEEE Trans. Geosci. Rem. Sens., " + \
    "44, 1102-1105, doi:10.1109/TGRS.2005.863716, 2006."
  gAtt@title = OBSNAME + " satellite data CMORized for the ESMValTool v2.0"
  gAtt@conventions = "CF/CMOR"

  ; Outfile
  fout = output_dir_path + "/" + \
    "OBS_" + OBSNAME + "_sat_1_" + FIELD + "_tro3_" + \
    YEAR1 + "01-" + YEAR2 + "12.nc"

  ; Write variable
  write_nc(fout, "tro3", output, gAtt)
  delete(gAtt)
  delete(output)

end
