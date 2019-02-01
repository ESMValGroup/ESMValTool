; #############################################################################
; ESMValTool CMORizer for ESACCI-LANDCOVER data
; #############################################################################
;
; Tier
;   Tier 2: other freely-available dataset.
;
; Source
;   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/land_cover/data/land_cover_maps/
;
; Last access
;   20190110
;
; Download and processing instructions
;   Download the 3 NetCDF files for 2000, 2005 and 2010.
;   Download the CCI-LC Tools from:
;     http://maps.elie.ucl.ac.be/CCI/viewer/download/lc-user-tools-3.14.zip
;   Unpack and run the CCI-LC Tools on each of the NetCDF files as follows:
;      bash lc-user-tools-3.14/bin/aggregate-map.sh \
;        -PgridName=GEOGRAPHIC_LAT_LON -PnumMajorityClasses=1 \
;        -PoutputAccuracy=false -PoutputPFTClasses=true \
;        -PoutputLCCSClasses=false -PnumRows=360 <inputfile.nc>
;   Put the resulting processed data in input_dir_path.
;
; Caveat
;   The CCI-LC Tools must be applied before running this script.
;   The CCI-LC Tools require Java Version 7 or higher.
;   The input data are available for a single year and are copied over to
;   generate a time series over their time range of validity.
;
; Modification history
;   20190110-A_righ_ma: rewritten in NCL for v2.
;   20160714-A_muel_bn: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-LANDCOVER.ncl"

  ; Source name
  OBSNAME = "ESACCI-LANDCOVER"

  ; Tier
  TIER = 2

  ; Years
  YEARS = (/2000, 2005, 2010/)

  ; Version
  VERSION = "L4-LCCS-Map-300m-P5Y-aggregated-0.500000Deg"

  ; Field
  FIELD = "T2Ms"

  ; CMOR table
  CMOR_TABLE = \
    getenv("esmvaltool_root") + "/cmor/tables/cmip5/Tables/CMIP5_Lmon"

  ; Variable names
  VARIABLES = \
    (/"baresoilFrac", "cropFrac", "grassFrac", "shrubFrac", "treeFrac"/)

  ; Corresponding aggregation classes in the raw data
  CLASSES = [/"Bare_Soil", \
              "Managed_Grass", \
              "Natural_Grass", \
              (/"Shrub_Broadleaf_Deciduous",    \
                "Shrub_Broadleaf_Evergreen",    \
                "Shrub_Needleleaf_Evergreen"/), \
              (/"Tree_Broadleaf_Deciduous",  \
                "Tree_Broadleaf_Evergreen",  \
                "Tree_Needleleaf_Deciduous", \
                "Tree_Needleleaf_Evergreen"/)/]

end

begin

  do yy = 0, dimsizes(YEARS) - 1

    fname = \
      input_dir_path + "/ESACCI-LC-" + VERSION + "-" + YEARS(yy) + "-v1.6.1.nc"

    f = addfile(fname, "r")

    ; Create time coordinate
    YEAR1 = YEARS(yy) - 2
    YEAR2 = YEARS(yy) + 2
    time = create_timec(YEAR1, YEAR2)

    do vv = 0, dimsizes(VARIABLES) - 1

      log_info("Processing variable " + VARIABLES(vv))

      ; Set classes to be added up
      class = CLASSES[vv]

      ; Save mask before adding up classes
      do cc = 0, dimsizes(class) - 1
        qq = f->$class(cc)$
        replace_ieeenan(qq, FILL, 0)
        qq@_FillValue = FILL
        tmp = ismissing(qq)
        delete(qq)
        if (cc.eq.0) then
          lmask = tmp
        else
          lmask := lmask .and. tmp
        end if
        delete(tmp)
      end do

      ; Add up classes
      do cc = 0, dimsizes(class) - 1
        log_info("  adding class " + class(cc))
        tmp = f->$class(cc)$
        replace_ieeenan(tmp, FILL, 0)
        tmp@_FillValue = FILL
        tmp = where(ismissing(tmp), 0., tmp)
        if (cc.eq.0) then
          xx = tmp
        else
          xx = xx + tmp
        end if
        delete(tmp)
      end do
      delete(class)

      ; Reapply mask of missing values
      xx = where(lmask, xx@_FillValue, xx)

      ; Define output field
      output = \
        new((/dimsizes(time), dimsizes(xx&lat), dimsizes(xx&lon)/), float)
      output!0 = "time"
      output&time = time
      output!1 = "lat"
      output&lat = xx&lat
      output!2 = "lon"
      output&lon = xx&lon
      output = conform(output, xx, (/1, 2/))
      delete(xx)

      ; Set standard fill value
      output@_FillValue = FILL

      ; Convert units
      output = output * 100
      output@units = "%"

      ; Format time coordinate
      ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
      delete(output&time)
      output&time = ctime
      delete(ctime)

      ; Format latitude coordinate
      output = output(:, ::-1, :)
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

      ; Set variable attributes
      tmp = var_attrib(output, VARIABLES(vv), CMOR_TABLE)
      delete(output)
      output = tmp
      delete(tmp)

      ; Special case for baresoilFrac: add auxiliary coordinate
      if (VARIABLES(vv).eq."baresoilFrac") then
        output@coordinates = "type"
      end if

      ; Set global attributes
      gAtt = True
      gAtt@history       = "Created on " + systemfunc("date")
      gAtt@host          = systemfunc("echo $HOST")
      gAtt@user          = systemfunc("echo $USER")
      gAtt@period        = YEAR1 + "-" + YEAR2
      gAtt@field         = FIELD
      gAtt@tier          = TIER
      gAtt@source        = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/" + \
        "land_cover/data/land_cover_maps/"
      gAtt@reference     = "Defourny, P., et al., 2015. ESA Land Cover " + \
        "Climate Change Initiative (ESA LC_cci) data"
      gAtt@title = OBSNAME + " satellite data CMORized for the ESMValTool v2.0"
      gAtt@conventions = "CF/CMOR"

      ; Outfile
      fout = output_dir_path + "/" + \
        "OBS_" + OBSNAME + "_sat_" + VERSION + "_" + FIELD + "_" + \
        VARIABLES(vv) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

      ; Write variable
      write_nc(fout, VARIABLES(vv), output, gAtt)
      delete(gAtt)
      delete(output)

      ; Special case for baresoilFrac: add auxiliary coordinate
      if (VARIABLES(vv).eq."baresoilFrac") then
        type = tochar("bare_ground")
        type!0 = "strlen"
        type@long_name = "surface type"
        type@standard_name = "area_type"
        w = addfile(fout, "w")
        w->type = type
        delete(w)
      end if

    end do
  end do

end
