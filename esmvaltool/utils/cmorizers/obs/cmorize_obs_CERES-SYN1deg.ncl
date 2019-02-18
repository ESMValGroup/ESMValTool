; #############################################################################
; ESMValTool CMORizer for CERES-SYN1deg data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://ceres-tool.larc.nasa.gov/ord-tool/jsp/SYN1degSelection.jsp
;
; Last access
;    20190207
;
; Download and processing instructions
;    Monthly data:
;      Expand "Compute TOA Fluxes" and select:
;        Shortwave Flux, Allsky and Clearsky
;        Longwave Flux, Allsky and Clearsky
;        Shortwave Flux Down, Allsky
;      Expand "Computed Surface Fluxes" and select:
;        Shortwave Flux Up, Allsky and Clearsky
;        Shortwave Flux Down, Allsky and Clearsky
;        Longwave Flux Up, Allsky and Clearsky
;        Longwave Flux Down, Allsky and Clearsky
;      then click on "Monthly", "Regional" and "Get data". All fields are saved
;      in CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Subset_200003-201702.nc
;    3hr data:
;      Select the same fields as above, then click on "Daily 3-Hourly" and
;      "Get data". All fields are saved in
;      CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_YYYYMMDD-YYYYMMDD.nc
;    Put all files in input_dir_path (no subdirectories with years).
;    For orders larger than 2 GB a registration is required.
;
; Modification history
;    20190124-A_righ_ma: written based on v1 different bash scripts by laue_ax
;                        and eval_ma.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_CERES-SYN1deg.ncl"

  ; Source name
  OBSNAME = "CERES-SYN1deg"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 2001
  YEAR2 = 2016

  ; CMOR name
  VARS = (/"rsus", "rsds", \
           "rlus", "rlds", \
           "rldscs", "rsut", \
           "rlut", "rsutcs", \
           "rlutcs", "rsdt", \
           "rsutcs", "rsut", \
           "rlutcs", "rlut", \
           "rsuscs", "rsus", \
           "rsdscs", "rsds", \
           "rluscs", "rlus", \
           "rldscs", "rlds"/)

  ; Name in the raw data
  NAME = (/"sfc_comp_sw_up_all_mon", "sfc_comp_sw_down_all_mon", \
           "sfc_comp_lw_up_all_mon", "sfc_comp_lw_down_all_mon", \
           "sfc_comp_lw_down_clr_mon", "toa_comp_sw_up_all_mon", \
           "toa_comp_lw_up_all_mon", "toa_comp_sw_up_clr_mon", \
           "toa_comp_lw_up_clr_mon", "toa_comp_sw_down_all_mon", \
           "toa_comp_sw_up_clr_3h", "toa_comp_sw_up_all_3h", \
           "toa_comp_lw_up_clr_3h", "toa_comp_lw_up_all_3h", \
           "sfc_comp_sw_up_clr_3h", "sfc_comp_sw_up_all_3h", \
           "sfc_comp_sw_down_clr_3h", "sfc_comp_sw_down_all_3h", \
           "sfc_comp_lw_up_clr_3h", "sfc_comp_lw_up_all_3h", \
           "sfc_comp_lw_down_clr_3h", "sfc_comp_lw_down_all_3h"/)

  ; MIP
  MIPS = (/"Amon", "Amon", \
           "Amon", "Amon", \
           "Amon", "Amon", \
           "Amon", "Amon", \
           "Amon", "Amon", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr"/)

  ; Version
  VERSION = "Ed3A"

  ; CMOR table: Amon version is used also for 3hr, since not all variables are
  ; available in 3hr (the tables are identical anyway)
  CMOR_TABLE = new(dimsizes(MIPS), string)
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_Amon"
  CMOR_TABLE(ind(VARS.eq."rluscs")) =  getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_rluscs.dat"

end

begin

  do vv = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vv) + " (" + MIPS(vv) + ")")

    do yy = YEAR1, YEAR2

      ; Read file
      if (MIPS(vv).eq."Amon") then
        f = addfile(input_dir_path + "CERES_SYN1deg-Month_Terra-Aqua-" + \
                    "MODIS_" + VERSION + "_Subset_200003-201702.nc", "r")
        xx = f->$NAME(vv)$
        delete(f)
        freq = "M"
      end if
      if (MIPS(vv).eq."3hr") then
        files = systemfunc("ls " + input_dir_path + "CERES_SYN1deg-3H_" + \
                           "Terra-Aqua-MODIS_" + VERSION + "_Subset_*" + \
                           yy + "*.nc")
        f = addfiles(files, "r")
        xx = f[:]->$NAME(vv)$
        delete(f)
        delete(files)
        freq = "3H"
      end if

      ; Time selection
      date = cd_calendar(xx&time, 0)
      output = xx(ind(date(:, 0).eq.yy), :, :)
      delete(date)
      delete(xx)

      ; Format time coordinate
      ctime = time_attrib(output&time, yy, yy, freq)
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
        "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/aerosol/data/", \
        "Wielicki et al., Bull. Amer. Meteor. Soc., " + \
        "doi: 10.1175/1520-0477(1996)077<0853:CATERE>2.0.CO;2, 1996", "")

      ; Outfile
      fout = output_dir_path + "OBS_" + OBSNAME + "_sat_" + VERSION + "_" + \
        MIPS(vv) + "_" + VARS(vv) + "_" + yy + "01-" + yy + "12.nc"

      ; Write variable
      write_nc(fout, VARS(vv), output, gAtt)
      delete(gAtt)
      delete(output)

    end do

  end do

end
