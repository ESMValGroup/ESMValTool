; #############################################################################
; ESMValTool CMORizer for WOA data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://data.nodc.noaa.gov/woa/WOA13/DATAv2/
;
; Last access
;    20190131
;
; Download and processing instructions
;    Download the following files:
;      temperature/netcdf/decav81B0/1.00/woa13_decav81B0_t00_01.nc
;      salinity/netcdf/decav81B0/1.00/woa13_decav81B0_s00_01.nc
;      oxygen/netcdf/all/1.00/woa13_all_o00_01.nc
;      nitrate/netcdf/all/1.00/woa13_all_n00_01.nc
;      phosphate/netcdf/all/1.00/woa13_all_p00_01.nc
;      silicate/netcdf/all/1.00/woa13_all_i00_01.nc
;
; Modification history
;    20190131-A_righ_ma: adapted to v2.
;    20190131-A_demo_le: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_WOA.ncl"

  ; Source name
  OBSNAME = "WOA"

  ; Tier
  TIER = 2

  ; Period
  YEAR = 2000

  ; Selected variables (standard name)
  VARS = (/"thetao", "so", "o2", "no3", "po4", "si"/)

  ; Name in the raw data
  NAMES = (/"t_an", "s_an", "o_an", "n_an", "p_an", "i_an"/)

  ; File names
  FNAMES = (/"woa13_decav81B0_t00_01.nc", \
             "woa13_decav81B0_s00_01.nc", \
             "woa13_all_o00_01.nc", \
             "woa13_all_n00_01.nc", \
             "woa13_all_p00_01.nc", \
             "woa13_all_i00_01.nc"/)

  ; MIP
  MIPS = (/"Omon", "Omon", "Oyr", "Oyr", "Oyr", "Oyr"/)

  ; CMOR table
  CMOR_TABLE = \
    getenv("esmvaltool_root") + "/cmor/tables/cmip5/Tables/CMIP5_" + MIPS


end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vv))

    f = addfile(input_dir_path + "/" + FNAMES(vv), "r")

    output = f->$NAMES(vv)$

    ; Convert units
    if (VARS(vv).eq."thetao") then
      output = output + 273.15  ; [degC] --> [K]
    end if

    if (VARS(vv).eq."o2") then
      output = output * 44.661 / 1000.  ; [ml L-1] --> [mol m-3]
    end if

    if (any(VARS(vv).eq.(/"no3", "po4", "si"/))) then
      output = output / 1000.  ; [ml L-1] --> [mol m-3]
    end if

    ; Set variable attributes
    tmp = var_attrib(output, VARS(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    delete(tmp)

    ; Format time coordinate
    opt = 0
    opt@calendar = "standard"
    ctime = todouble(output&time)  ; this also removes attributes
    ctime = cd_inv_calendar(YEAR, 7, 1, 0, 0, 0, TUNITS, opt)
    ctime@long_name = "time"
    ctime@axis = "T"
    ctime@units = TUNITS
    ctime@standard_name = "time"
    if (isatt(ctime, "_FillValue")) then
      delete(ctime@_FillValue)
    end if
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ; Format vertical coordinate
    output!1 = "lev"
    clev = lev_attrib(output&lev)
    delete(output&lev)
    output&lev = clev
    delete(clev)

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
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "https://data.nodc.noaa.gov/woa/WOA13/DATAv2/", \
      "Locarnini et al., World Ocean Atlas 2013, Vol. 1: Temperature, 2013", \
      "")

    ; Outfile
    fout = output_dir_path + \
      "OBS_" + OBSNAME + "_clim_L3_" + MIPS(vv) + "_" + VARS(vv) + "_" + \
      YEAR + "01-" + YEAR + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vv), output, gAtt)
    delete(gAtt)
    delete(output)

  end do

end
