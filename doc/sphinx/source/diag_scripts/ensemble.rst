:mod:`ensemble`
===============
.. function:: get_start_year(model_index[1]:numeric)

   :param numeric model_index: index of array "models"; this function is supposed to be called inside a model loop

   Return value
      An integer containing the first year
  
   Description
      Chooses first year of model input with respect to
      diag_script_info@range_option
  
   Caveats
  
   References
  
   Modification history
      20140128-A_senf_da: written.
  
.. function:: get_end_year(model_index:numeric)

   :param numeric model_index: index of array "models"; this function is supposed to be called inside a model loop

   Return value
      An integer containing the last year
  
   Description
      Chooses last year of model input with respect to
      diag_script_info@range_option
  
   Caveats
  
   References
  
   Modification history
      20140128-A_senf_da: written.
  
.. function:: multi_model_stats(idx[*]:integer, var:string, field:string, grd[1]:string, statistics[*]:string, opt[1])

   :param integer idx: a 1-D array of model indexes to be considered.
   :param string var: variable name.
   :param string field: field type.
   :param string grd: type of grid to be selected: "coarsest": returns the lowest resolution grid. "finest": returns the highest resolution grid.
   :param string statistics: a 1-D array of strings for selecting the statistical fields to return. Currently implemented: "model_mean" -> multi model mean "stddev" -> standard deviation over "models" dimension "model_grd" -> models regridded to a common grid
   :param integer opt: not used yet, set to 0

   Return value
      A field with target grid spatial dimensions (lat, lon) or
      (plev, lat, lon), time coordinate, models coordinate.
      The "models" coordinate contains the strings from models@name and
      statistics.
      Covers the overlapping time period of all models only.
  
   Description
      Determines statistical fields of all input models, after regridding them
      to a common grid. Statistics is over the "models" coordinate only.
      See function find_destination_grid (regridding.ncl) for details of
      target grid selection.
  
   Caveats
      Works for rectilinear grids only, but selection of regridding routines
      might be extended (e.g. select automatically, depending on input).
      Might be extended to return multi year monthly mean values.
      Wrap to add result to the structure returned by function read_data?
      Probably not very efficient: try to re-use intermediate data etc.
  
   References
  
   Modification history:
      20140515-A_gott_kl: written.
  
