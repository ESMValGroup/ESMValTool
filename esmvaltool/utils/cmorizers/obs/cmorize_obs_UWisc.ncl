; #############################################################################
; ESMValTool CMORizer for UWisc data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    Data provided by Ralf Bennartz.
;
; Last access
;    20150415
;
; Download and processing instructions
;    Contact Ralf Bennartz (Earth and Environmental Sciences, Vanderbilt
;    University, USA).
;
; Modification history
;    20190208-A_righ_ma: adapted to v2.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")
begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_UWISC.ncl"

  ; Source name
  OBSNAME = "UWisc"

  ; Tier
  TIER = 3

  ; Field
  FIELD = "T2Ms"

  ; Period
  YEAR1 = 1988
  YEAR2 = 2007

  ; MIP
  MIP = "Amon"

  ; Selected variable (standard name)
  VARS = (/"clwvi", "lwpStderr"/)

  ; Name in the raw data
  NAME = (/"LWP", "LWP_ERROR"/)

  ; Conversion factor
  CONVERSION = (/1.e-3, 1.e-3/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    (/"/cmor/tables/cmip5/Tables/CMIP5_Amon", \
      "/cmor/tables/custom/CMOR_lwpStderr.dat"/)

end

begin

  do vv = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vv))

    do yr = YEAR1, YEAR2

      fname = input_dir_path + "UWisc_LWPMEAN_" + yr + "_v2.nc"

      f = addfile(fname, "r")

      output = (/f->$NAME(vv)$/)

      lat = (/f->lat/)
      lon = (/f->lon/)

      ; Format time coordinate
      output!0 = "time"
      output&time = create_timec(yr, yr)
      ctime = time_attrib(output&time, yr, yr, "M")
      delete(output&time)
      output&time = ctime
      delete(ctime)

      ; Format latitude coordinate
      output!1 = "lat"
      output&lat = lat
      if (isMonotonic(output&lat) .eq. -1) then
        output = output(:, ::-1, :)
      end if
      clat = lat_attrib(output&lat)
      delete(output&lat)
      output&lat = clat
      delete(clat)

      ; Format longitude coordinate
      output!2 = "lon"
      output&lon = lon
      if (any(output&lon.lt.0.)) then
        output = lonFlip(output)
      end if
      clon = lon_attrib(output&lon)
      delete(output&lon)
      output&lon = clon
      delete(clon)

      ; Convert units
      output = output * CONVERSION(vv)

      ; Set variable attributes
      tmp = var_attrib(output, VARS(vv), CMOR_TABLE(vv))
      delete(output)
      output = tmp
      delete(tmp)

      ; Set global attributes
      gAtt = set_global_atts( \
        OBSNAME, TIER, \
        "Data provided by Ralf Bennartz (Vanderbilt University, USA)", \
        "O'Dell et al., J. Clim., doi:10.1175/2007JCLI1958.1, 2008", "")
      gAtt@model_id = OBSNAME
      gAtt@project_id = "OBS"

      ; Outfile
      fout = output_dir_path + "OBS_" + OBSNAME + "_sat_v2_" + MIP + \
        "_" + VARS(vv) + "_" + yr + "01-" + yr + "12.nc"

      ; Write variable
      write_nc(fout, VARS(vv), output, gAtt)
      delete(gAtt)
      delete(output)

    end do

  end do

end
