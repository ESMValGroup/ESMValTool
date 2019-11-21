; #############################################################################
; # GENERAL ROUTINES FOR REGRIDDING OPERATIONS
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function guestimate_average_grid_area
;    function get_lower_limits
;    function get_upper_limits
;    function is_regional
;    function esmf_conserve_wrapper
;    function rect2rect_interp
;    function plev_lat_interp
;    function get_dataset_minus_ref
;    function esmf_conserve_wrapper_time
;    function regrid_3D_to_rectilinear_grid
;
; #############################################################################


load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/constants.ncl"
load "$diag_scripts/../interface_scripts/data_handling.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; #############################################################################
undef("guestimate_average_grid_area")
function guestimate_average_grid_area(data[*][*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, dimsize_lat, dimsize_lon, average_lat, average_lon
begin

  funcname = "guestimate_average_grid_area"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  dimsize_lat = dimsizes(data&lat)
  average_lat = abs(avg(data&lat(1:dimsize_lat - 1) - \
                        data&lat(0:dimsize_lat - 2)))

  if (isdim(data, "lon")) then
    dimsize_lon = dimsizes(data&lon)
    average_lon = abs(avg(data&lon(1:dimsize_lon - 1) - \
                          data&lon(0:dimsize_lon - 2)))
    average_cell = average_lat * average_lon
  end if

  if (isdim(data, "plev")) then
    dimsize_plev = dimsizes(data&plev)
    average_plev = abs(avg(data&plev(1:dimsize_plev - 1) - \
                           data&plev(0:dimsize_plev - 2)))
    average_cell = average_lat * average_plev
  end if

  leave_msg(scriptname, funcname)
  return(average_cell)

end

; #############################################################################
undef("get_lower_limits")
function get_lower_limits(coordinate[*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, d_size, delta_sizes, lower_limits
begin

  funcname = "get_lower_limits"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  d_size = dimsizes(coordinate)
  delta_sizes = abs((coordinate(1:d_size - 1) - \
                     coordinate(0:d_size - 2)) / 2)
  lower_limits = new(d_size, typeof(coordinate))
  lower_limits(1:d_size - 1) = coordinate(1:d_size - 1) - delta_sizes
  lower_limits(0) = coordinate(0) - delta_sizes(0)

  leave_msg(scriptname, funcname)
  return(lower_limits)

end

; #############################################################################
undef("get_upper_limits")
function get_upper_limits(coordinate[*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, d_size, delta_sizes, upper_limits
begin

  funcname = "get_upper_limits"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  d_size = dimsizes(coordinate)
  delta_sizes = abs((coordinate(1:d_size - 1) - \
                     coordinate(0:d_size - 2)) / 2)
  upper_limits = new(d_size, typeof(coordinate))
  upper_limits(0:d_size - 2) = coordinate(0:d_size - 2) + delta_sizes
  upper_limits(d_size - 1) = coordinate(d_size - 1) + delta_sizes(d_size - 2)

  leave_msg(scriptname, funcname)
  return(upper_limits)

end

; #############################################################################
undef("is_regional")
function is_regional(grid : numeric)
;
; Arguments
;    grid: input grid with lat/lon coordinates
;
; Return value
;    logical indicitating whether it is a global (=.False.) or
;    regional grid (=.True.)
;
; Description
;    Run a test to estimate whether the grid at hand is global or
;    regional.
;
; Caveats
;
; Reference
;
; Modification history
;
local funcname, scriptname, lat, lon, lat_size, lon_size, regional_area
begin

  funcname = "is_regional"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  lat = grid&lat
  lon = grid&lon
  lat_size = dimsizes(lat)
  lon_size = dimsizes(lon)
  regional_area = False

  if (abs(lat(lat_size - 1) - lat(0)) .lt. 80) then
    regional_area = True
  end if

  if (abs(lon(lon_size - 1) - lon(0)) .lt. 170) then
    regional_area = True
  end if

  leave_msg(scriptname, funcname)
  return(regional_area)

end

; #############################################################################
undef("esmf_conserve_wrapper")
function esmf_conserve_wrapper(source[*][*]:numeric,
                               destination[*][*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, \
  data_regrid,     d_destination,   d_source, \
  dstGridCorner,   dst_lower_lat,   dst_lower_lon, \
  dst_upper_lat,   dst_upper_lon,   Opt, \
  src_lower_lat,   src_lower_lon, \
  src_upper_lat,   src_upper_lon,   src_regional, \
  dst_regional
begin

  funcname = "esmf_conserve_wrapper"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  Opt                    = True
  Opt@SrcGridLat         = source&lat
  Opt@SrcGridLon         = source&lon
  Opt@DstGridLat         = destination&lat
  Opt@DstGridLon         = destination&lon
  Opt@DstGridType        = source@GridType
  Opt@SrcGridType        = destination@GridType

  log_debug("esmf_conserve_wrapper: SrcGridType is " + source@GridType)
  log_debug("esmf_conserve_wrapper: DstGridType is " + destination@GridType)

  Opt@ForceOverwrite    = True
  Opt@PrintTimings      = False
  Opt@Debug             = False

  ; "bilinear" ; "patch"  ; "conserve"
  Opt@InterpMethod      = "conserve"

  if (.not. fileexists(regridding_dir)) then
    system("mkdir -p " + regridding_dir)
  end if
  Opt@WgtFileName = regridding_dir + "/" + "/weight_file.nc"
  Opt@SrcFileName = regridding_dir + "/" + "/source_file.nc"
  Opt@DstFileName = regridding_dir + "/" + "/destination_file.nc"

  Opt@printTimings = False
  ; ESMF_regrid may issue a warning
  ; (e.g. if Dst is not covered Src grid)
  ; -> set "exit_on_warning" to False (in main nml) to continue to the
  ; next diag

  ; Setting Opt@SrcRegional=True avoids artefacts
  ; from incomplete Src grids
  ; Activate if only part of the sphere is given
  Opt@SrcRegional = source@Regional
  Opt@DstRegional = destination@Regional

  data_regrid = ESMF_regrid(source, Opt)

  leave_msg(scriptname, funcname)
  return(data_regrid)

end

; #############################################################################
undef("rect2rect_interp")
function rect2rect_interp(source[*][*]:numeric,
                          target[*][*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;    Interpolates rectangular grid source (high res) onto target grid
;    (low res) using local area averages.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, interpolated

begin

  funcname = "rect2rect_interp"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  interpolated = esmf_conserve_wrapper(source, target)

  leave_msg(scriptname, funcname)
  return(interpolated)

end

; #############################################################################
undef("plev_lat_interp")
function plev_lat_interp(source[*][*]:numeric,
                         target[*][*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;    Interpolates plev/lat grid source (high res) onto target grid
;    (low res) using local linear interpolation
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, interpolated, grid_type
begin

  funcname = "plev_lat_interp"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  interpolated = linint2(source&lat, source&plev, source, \
                         False, target&lat, target&plev, 0)

  leave_msg(scriptname, funcname)
  return(interpolated)

end

; #############################################################################
undef("get_dataset_minus_ref")
function get_dataset_minus_ref(dataset[*][*]:numeric,
                               ref[*][*]:numeric)
;
; Arguments
;
; Return value
;
; Description
;    Interpolates rectangular grid source onto target grid and returns their
;    difference.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, dataset_grid_size, ref_grid_size, local_dataset, \
  local_ref, difference
begin

  funcname = "get_dataset_minus_ref"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  is_plev = any(isdim(dataset, "plev"))
  dataset_grid_size = guestimate_average_grid_area(dataset)
  ref_grid_size = guestimate_average_grid_area(ref)
  dimdataset = dimsizes(dataset)
  dimref = dimsizes(ref)
  log_debug("get_dataset_minus_ref: dims are, dataset:" + dimdataset + \
            ", ref" + dimref)

  if (dataset_grid_size .gt. ref_grid_size) then
    if (is_plev) then
      local_ref = plev_lat_interp(ref, dataset)
    else
      local_ref = rect2rect_interp(ref, dataset)
    end if
    local_dataset = dataset
    ; Copy meta data
    difference = local_dataset
  else
    if (is_plev) then
      local_ref = plev_lat_interp(dataset, ref)
    else
      local_dataset = rect2rect_interp(dataset, ref)
    end if
    local_ref = ref
    ; Copy meta data
    difference = local_ref
  end if
  dimdataset = dimsizes(local_dataset)
  dimref = dimsizes(local_ref)
  log_debug("get_dataset_minus_ref: dims are, dataset:" + dimdataset + \
            ", ref" + dimref)
  difference = local_dataset - local_ref

  leave_msg(scriptname, funcname)
  return(difference)

end

; #############################################################################
undef("esmf_conserve_wrapper_time")
function esmf_conserve_wrapper_time(source[*][*][*]:numeric,
                                    destination[*][*][*]:numeric,
                                    weight_file[1]:string,
                                    source_file[1]:string,
                                    destination_file[1]:string)
;
; Arguments
;
; Return value
;
; Description
;
; Caveats
;    Assumes regular grid.
;
; References
;
; Modification history
;
local funcname, scriptname, \
  data_regrid,    d_destination,  dstGridCorner,\
  dst_lower_lat,  dst_lower_lon,  dst_upper_lat,\
  dst_upper_lon,  funcname,       Opt,\
  scriptname,     srcGridCorner,  src_lower_lat,\
  src_lower_lon,  src_upper_lat,  src_upper_lon
begin

  funcname = "esmf_conserve_wrapper_time"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  Opt                    = True
  Opt@SrcGridLat         = source&lat
  Opt@SrcGridLon         = source&lon
  Opt@DstGridLat         = destination&lat
  Opt@DstGridLon         = destination&lon

  Opt@ForceOverwrite    = True
  Opt@PrintTimings      = False
  Opt@Debug             = False

  ; "bilinear" ; "patch"  ; "conserve"
  Opt@InterpMethod      = "conserve"

  if (.not. fileexists(regridding_dir)) then
    system("mkdir -p " + regridding_dir)
  end if

  ; Save regrid files for later use: careful with naming!
  Opt@WgtFileName = regridding_dir + "/" + weight_file + ".nc"
  Opt@SrcFileName = regridding_dir + "/" + source_file + ".nc"
  Opt@DstFileName = regridding_dir + "/" + destination_file + ".nc"

  Opt@printTimings = False
  ; ESMF_regrid may issue a warning
  ; (e.g. if Dst is not covered Src grid)
  ; -> set "exit_on_warning" to False (in main nml) to continue to the
  ; next diag

  ; Setting Opt@SrcRegional=True avoids artefacts
  ; from incomplete Src grids
  Opt@SrcRegional = False  ; Activate if only part of the sphere is given
  Opt@DstRegional = False  ; Activate if only part of the sphere is given

  ; reuse regrid files if present
  if (fileexists(Opt@SrcFileName)) then
    Opt@SkipSrcGen = True
  end if
  if (fileexists(Opt@DstFileName)) then
    Opt@SkipDstGen = True
  end if

  ntimes = dimsizes(source(:, 0, 0))  ; get number of time steps

  Opt@SrcMask2D = where(.not.ismissing(source(0, :, :)), 1, 0)
  test_regrid = ESMF_regrid(source(0, :, :), Opt)  ; only to get dim sizes
  newdims = dimsizes(test_regrid)

  dims_3d = (/ntimes, newdims(0), newdims(1)/)
  source_regrid = new(dims_3d, float, test_regrid@_FillValue)

  do t = 0, ntimes - 1
    SrcMask2D = where(.not.ismissing(source(t, :, :)), 1, 0)
    if (any(SrcMask2D - Opt@SrcMask2D .gt. 0)) then
      Opt@SrcMask2D  = SrcMask2D
      source_regrid(t, :, :) = ESMF_regrid(source(t, :, :), Opt)
    else
      source_regrid(t, :, :) = \
        ESMF_regrid_with_weights(source(t, :, :), Opt@WgtFileName, Opt)
    end if
  end do
  copy_VarCoords_2(source, source_regrid)

  leave_msg(scriptname, funcname)
  return(source_regrid)

end

; #############################################################################
undef("regrid_3D_to_rectilinear_grid")
function regrid_3D_to_rectilinear_grid(data_in:numeric,
                                       lon:numeric,
                                       lat:numeric,
                                       grid_resolution[1]:string,
                                       filename[1]:string,
                                       regular[1]:logical)
;
; Arguments
;    data_in: 3D field array with imput data
;    lon: array with longitudes
;    lat: array with latitudes
;    grid_resolution: grid resolution of destination grid
;    filename: file name of dataset file
;    regular: defines grid type
;             True: rectilinear
;             False: curvilinear
;
; Return value
;    An 3D array with new dimensions
;
; Description
;
; Caveats
;    It seems to not work properly with irregular grids.
;
; References
;
; Modification history
;    20151026-righi_mattia: added warning for unavailable lat/lon vertices
;                           in input.
;    20151023-righi_mattia: moved to regridding.ncl and renamed
;                           regrid_3D_data_to_global_rectilinear_grid -->
;                           regrid_3D_to_rectilinear_grid.
;    20150703-wenzel_sabrina: moved to anav13jclim_func.ncl and adapted to
;                             ESMValTool structure.
;    201505??-anav_alessandro: written.
;
local funcname, scriptname, INTERP_METHOD, srcGridName, dstGridName, dimx, \
  ntime, nlat, nlon, temp, Opt, sfile, var
begin

  funcname = "regrid_3D_to_rectilinear_grid"
  scriptname = "diag_scripts/shared/regridding.ncl"
  enter_msg(scriptname, funcname)

  ; Workdir
  workpath = config_user_info@work_dir

  ; Define interpolation method
  INTERP_METHOD = (/"conserve"/)

  ; Temporary files
  srcGridName = workpath + "/CMIP5_SCRIP.nc"
  dstGridName = workpath + "/World_SCRIP.nc"
  wgtFile = workpath + "/CMIP5_World.nc"

  ; Retrieving dimensions of input data
  dimx = dimsizes(data_in)
  ntime = dimx(0)
  nlat = dimx(1)
  nlon = dimx(2)

  ; Convert datasets to SCRIP file.
  temp = dim_sum_n(data_in, 0)  ; This is needed to create the mask

  ; Variables like LAI have a different distribution of NaNs
  ; between winter and summer months

  Opt = True

  Opt@SrcMask2D = where(.not.ismissing(temp), 1, 0)
  Opt@SrcRegional = True  ; Setting True to avoid artefacts from incomplete
                          ; Src grids. Activate if only part of the sphere
                          ; is given
  Opt@PrintTimings = False
  Opt@ForceOverwrite = True
  Opt@Debug = False
  delete(temp)

  if (regular.eq.True) then
    rectilinear_to_SCRIP(srcGridName, lat, lon, Opt)
  else
    sfile = addfile(filename, "r")
    if (all(isfilevar(sfile, (/"lat_vertices", "lon_vertices"/)))) then
      Opt@GridCornerLat = sfile->lat_vertices
      Opt@GridCornerLon = sfile->lon_vertices
    else
      error_msg("w", "regridding.ncl", funcname, \
                "lat_vertices/lon_vertices not available for this " + \
                "dataset, using ESMF estimated values")
    end if
    curvilinear_to_SCRIP(srcGridName, lat, lon, Opt)
  end if
  delete(Opt)

  ; convert to SCRIP file
  Opt = True
  Opt@ForceOverwrite = True
  Opt@PrintTimings = False
  Opt@Debug = False

  if (grid_resolution.eq."0.5deg") then
    dNLON = 720
    dNLAT = 360
    Opt@DstLLCorner = (/-89.75, 0.75 /)
    Opt@DstURCorner = (/ 89.75, 359.75 /)
    latlon_to_SCRIP(dstGridName, "0.5deg", Opt)
  else if (grid_resolution.eq."1deg") then
    dNLON = 360
    dNLAT = 180
    Opt@DstLLCorner = (/-89.75, 0.75 /)
    Opt@DstURCorner = (/ 89.75, 359.75 /)
    latlon_to_SCRIP(dstGridName, "1x1", Opt)
  else
    dNLON = 180
    dNLAT = 90
    Opt@LLCorner = (/-89.5, 0.00/)
    Opt@URCorner = (/89.5, 359.5/)
    latlon_to_SCRIP(dstGridName, "2x2", Opt)
  end if
  end if
  delete(Opt)

  ; Generate weights
  Opt = True
  Opt@ForceOverwrite = True
  Opt@SrcESMF = False
  Opt@DstESMF = False
  Opt@remap_method = INTERP_METHOD
  Opt@InterpMethod = INTERP_METHOD
  Opt@PrintTimings = False

  ESMF_regrid_gen_weights(srcGridName, dstGridName, wgtFile, Opt)
  delete(Opt)

  ; Interpolating datasets from original grid final grid
  Opt = True
  Opt@PrintTimings = False
  data_out = new((/ntime, dNLAT, dNLON/), typeof(data_in))
  copy_VarCoords_2(data_in, data_out)

  ; Regrid over time data
  do t = 0, ntime - 1
    var = data_in(t, :, :)
    data_out(t, :, :) = ESMF_regrid_with_weights(var, wgtFile, Opt)
  end do

  system("rm -f " + srcGridName + " " + dstGridName + " " + \
         wgtFile + " " + "PET0.RegridWeightGen.Log")

  leave_msg(scriptname, funcname)
  return(data_out)

end
