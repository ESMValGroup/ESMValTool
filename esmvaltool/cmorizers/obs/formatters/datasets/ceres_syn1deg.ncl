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
;    20190124-righi_mattia: written based on v1 different bash scripts by
;                           Axel Lauer and Martin Evaldsson.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "ceres_syn1deg.ncl"

  ; Source name
  OBSNAME = "CERES-SYN1deg"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 2001
  YEAR2 = 2016

  ; CMOR name
  VAR = (/"rsus", "rsds", \
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
  MIP = (/"Amon", "Amon", \
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

  ; Frequency
  FREQ = (/"mon", "mon", \
           "mon", "mon", \
           "mon", "mon", \
           "mon", "mon", \
           "mon", "mon", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr", \
           "3hr", "3hr"/)

  ; CMOR table: Amon version is used also for 3hr, since not all variables are
  ; available in 3hr (the tables are identical anyway)
  CMOR_TABLE = new(dimsizes(MIP), string)
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_Amon"
  CMOR_TABLE(ind(VAR.eq."rluscs")) =  \
    getenv("cmor_tables") + "/custom/CMOR_rluscs.dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "Ed3A"

  ; Global attributes
  SOURCE = "https://ceres-tool.larc.nasa.gov/ord-tool/jsp/SYN1degSelection.jsp"
  REF = "Wielicki et al., Bull. Amer. Meteor. Soc., " + \
    "doi: 10.1175/1520-0477(1996)077<0853:CATERE>2.0.CO;2, 1996"
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yy = YEAR1, YEAR2

      ; Read file
      if (MIP(vv).eq."Amon") then
        f = addfile(input_dir_path + "CERES_SYN1deg-Month_Terra-Aqua-" + \
                    "MODIS_" + VERSION + "_Subset_200003-201702.nc", "r")
        xx = f->$NAME(vv)$
        delete(f)
      end if
      if (MIP(vv).eq."3hr") then
        files = systemfunc("ls " + input_dir_path + "CERES_SYN1deg-3H_" + \
                           "Terra-Aqua-MODIS_" + VERSION + "_Subset_*" + \
                           yy + "*.nc")
        f = addfiles(files, "r")
        xx = f[:]->$NAME(vv)$
        delete(f)
        delete(files)
      end if

      ; Time selection
      date = cd_calendar(xx&time, 0)
      output = xx(ind(date(:, 0).eq.yy), :, :)
      delete(date)
      delete(xx)

      ; Format coordinates
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      format_coords(output, yy + "0101", yy + "1231", FREQ(vv))

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
      DATESTR = yy + "01-" + yy + "12"
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

    end do

  end do

end
