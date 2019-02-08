; #############################################################################
; REFORMAT SCRIPT FOR HadISST REANALYSIS DATA
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
;
; Last access
;    02/2019
;
; Download and processing instructions
;    Download and unzip "HadISST_ice.nc.gz" and "HadISST_sst.nc.gz".
;
; Modification history
;    20190208-A_hass_bg: adapted to v2.
;    20180530-A_Righ_ma: fixed coordinates and metadata
;    20170217-A_senf_da: modified so that SST fields are also written as 'tos'
;    20150422-A_laue_ax: written
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_HadISST.ncl"

  ; Source name
  OBSNAME = "HadISST"
  
  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1870
  YEAR2 = 2017
  
  ; MIP
  MIP = (/"Amon", "Omon", "OImon"/)
  
  ; Selected variable (standard name)
  VARS = (/"ts", "tos", "sic"/)
  
  ; Name in the raw data
  NAME = (/"sst", "sst", "ice"/)
  
  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_" + MIPS

end
; #############################################################################
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

  ; Loop over variables
  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    initialize = 1
    if (isStrSubset(NAME(vID), "sst")) then
        fname = INDIR + "HadISST_sst.nc"
    end if
    if (isStrSubset(NAME(vID), "ice")) then
        fname = INDIR + "HadISST_ice.nc"
    end if

    f_var = addfile(fname, "r")
    
    lat  = (/f_var->latitude/)
    lon  = (/f_var->longitude/)
    time = (/f_var->time/)

    lat@standard_name = "latitude"
    lat@units = "degrees_north"
    lat@axis = "Y"
    lat@long_name = "latitude"
    lat!0 = "lat"

    lon@standard_name = "longitude"
    lon@units = "degrees_east"
    lon@axis = "X"
    lon@long_name = "longitude"
    lon!0 = "lon"

    time@standard_name = "time"
    time@units = "days since 1870-01-01 00:00:00"
    time@axis = "T"
    time@long_name = "time"
    time@calendar = "gregorian"
    if (isatt(time, "_FillValue")) then
        delete(time@_FillValue)
    end if

    if (isStrSubset(NAME(vID), "sst")) then
        new_var = (/f_var->$NAME(vID)$/)
        new_var@_FillValue = -1000
        new_var@missing_value = -1000
        new_var@units = "K"
        if (isStrSubset(VARS(vID), "ts")) then
            new_var@long_name = "Surface Temperature"
            new_var@standard_name = "surface_temperature"
        end if
        if (isStrSubset(VARS(vID), "tos")) then
            new_var@long_name = "Sea Surface Temperature"
            new_var@standard_name = "sea_surface_temperature"
        end if
    end if
    if (isStrSubset(NAME(vID), "ice")) then
        new_var = (/f_var->$NAME(vID)$/)
        new_var@_FillValue = -1.0e30
        new_var@missing_value = -1.0e30
        new_var@units = "%"
        new_var@long_name = "Sea Ice Area Fraction"
        new_var@standard_name = "sea_ice_area_fraction"
    end if

    new_var!0 = "time"
    new_var!1 = "lat"
    new_var!2 = "lon"
    new_var&time = time
    new_var&lat = lat
    new_var&lon = lon
        
    new_var = new_var(:, ::-1, :)
    new_var = lonFlip(new_var)

    if (isStrSubset(NAME(vID), "sst")) then
        ; convert from degrees C to K
        new_var = new_var + 273.15
    end if
    if (isStrSubset(NAME(vID), "ice")) then
        ; Convert from 1 to %
        ice = ice * 100.
    end if

    ; Set global attributes
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "http://www.metoffice.gov.uk/hadobs/hadisst/data/" + \
        "download.html", \
      "Rayner et al. (2013), J. Geophys. Res.", \
      "HadISST reanalysis data reformatted" + \
        " for the ESMValTool")

    if (isStrSubset(NAME(vID), "sst")) then
        ; name substring
        new_var_substr = "1"
    end if
    if (isStrSubset(NAME(vID), "ice")) then
        ; name substring (file naming following Dominik Kunert, DLR)
        new_var_substr = "20130524"
    end if
        
        
    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_reanaly_" + \
           new_var_substr + "_" + MIP + "_" + \
           VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"
    
    ; Write variable
    write_nc(fout, VARS(vID), new_var, gAtt)
    delete(gAtt)
    delete(new_var)

  end do
      

;    Reference:
;       Rayner, N. A., D. E. Parker, E. B. Horton, C. K. Folland, L. V.
;          Alexander, D. P. Rowell, E. C. Kent, A. Kaplan(2003), Global
;          analyses of sea surface temperature, sea ice, and night marine
;          air temperature since the late nineteenth century, J. Geophys.
;          Res., 108(D14), 4407, doi: 10.1029/2002JD002670.


end
