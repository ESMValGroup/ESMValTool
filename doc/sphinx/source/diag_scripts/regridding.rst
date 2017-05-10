:mod:`regridding`
=================
.. function:: find_destination_grid(indexes[*]:integer, var:string, field:string, opt[1]:string)

   :param integer indexes: a 1-D array of model indexes to be considered.
   :param string var: variable name.
   :param string field: field type.
   :param string opt: type of grid to be selected: "coarsest": returns the lowest resolution grid. "finest": returns the highest resolution grid.

   Return value
      A 2-D or 3-D dummy variable representing the grid with the attached
      plev/lat/lon or lat/lon coordinates.
  
   Description:
      Given an array of models, returns the coordinates of the coarsest or
      finest grid, to be used as a destination grid in regridding routines.
      For the vertical coordinate, the extent is considered as first priority
      (to avoid loss of data).
      All models are expect to have the same rank and dimension sizes.
  
   Caveats
      The returned plev, lat, lon coordinates do not necessarily come from the
      same input model, i.e. it is possible to get latitude from one model
      longitude from another model.
  
   Modification history
      20150113-A_gott_kl: check that all models have the same dimensions.
      20150113-A_gott_kl: revised to allow data without lev or lon coordinate.
      20140903-A_righ_ma: revised plev coordinate selection.
      20140311-A_righ_ma: added plev coordinate.
      20140212-A_righ_ma: written.
  
.. function:: guestimate_average_grid_area(data[*][*]:numeric)


   Return value
  
   Description
  
   Caveats
  
   References
  
   Modification history
  
.. function:: get_lower_limits(coordinate[*]:numeric)


   Return value
  
   Description
  
   Caveats
  
   References
  
   Modification history
  
.. function:: get_upper_limits(coordinate[*]:numeric)


   Return value
  
   Description
  
   Caveats
  
   References
  
   Modification history
  
.. function:: is_regional(grid : numeric)

   :param  numeric grid: input grid with lat/lon coordinates

   Return value
      logical indicitating whether it is a global (=.False.) or
      regional grid (=.True.)
  
   Description
      Run a test to estimate whether the grid at hand is global or
      regional.
  
   Caveats
  
   Reference
  
   Modification history
  
.. function:: esmf_conserve_wrapper(source[*][*]:numeric, destination[*][*]:numeric)


   Return value
  
   Description
  
   Caveats
  
   References
  
   Modification history
  
.. function:: rect2rect_interp(source[*][*]:numeric, target[*][*]:numeric)


   Return value
  
   Description
      Interpolates rectangular grid source (high res) onto target grid
      (low res) using local area averages.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: plev_lat_interp(source[*][*]:numeric, target[*][*]:numeric)


   Return value
  
   Description
      Interpolates plev/lat grid source (high res) onto target grid
      (low res) using local linear interpolation
  
   Caveats
  
   References
  
   Modification history
  
.. function:: get_model_minus_ref(model[*][*]:numeric, ref[*][*]:numeric)


   Return value
  
   Description
      Interpolates rectangular grid source onto target grid and returns their
      difference.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: esmf_conserve_wrapper_time(source[*][*][*]:numeric, destination[*][*][*]:numeric, weight_file[1]:string, source_file[1]:string, destination_file[1]:string)


   Return value
  
   Description
  
   Caveats
      Assumes regular grid.
  
   References
  
   Modification history
  
.. function:: regrid_3D_to_rectilinear_grid(data_in:numeric, lon:numeric, lat:numeric, grid_resolution[1]:string, filename[1]:string, regular[1]:logical)

   :param numeric data_in: 3D field array with imput data
   :param numeric lon: array with longitudes
   :param numeric lat: array with latitudes
   :param string grid_resolution: grid resolution of destination grid
   :param string filename: file name of model file
   :param logical regular: defines grid type True: rectilinear False: curvilinear

   Return value
      An 3D array with new dimensions
  
   Description
  
   Caveats
      It seems to not work properly with irregular grids.
  
   References
  
   Modification history
      20151026_A_righ_ma: added warning for unavailable lat/lon vertices
                          in input.
      20151023_A_righ_ma: moved to regridding.ncl and renamed
                          regrid_3D_data_to_global_rectilinear_grid -->
                          regrid_3D_to_rectilinear_grid.
      20150703_A_wenz_sa: moved to anav13jclim_func.ncl and adapted to
                          ESMValTool structure.
      201505??_A_anav_al: written.
  
