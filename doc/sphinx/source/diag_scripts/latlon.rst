:mod:`latlon`
=============
.. function:: roi(latlon_dat[4]:float, latlon_roi[4]:float)

   :param float latlon_dat: (/latmin,latmax,lonmin,lonmax/) of the (data) region to check
   :param float latlon_roi: (/latmin,latmax,lonmin,lonmax/) of the region of interest. Ranges: -90 < lat < 90, 0 < lon < 360

   Return value
      A string array containing all true statements from the following list:
      contained / center / overlap / encloses / outside (e.g. "contained"
      implies "center" and "overlap")
  
   Description
      Checks if the area described by latlon_dat is contained in / has its
      center in / overlaps with / encloses / is outside the area described by
      latlon_roi, for use with Emmons.ncl
  
   Caveats
      Not fully tested for areas containing Greenwich or the date line
  
   References
  
   Modification history
      20140129-A_gott_kl: written.
  
.. function:: extract_area(index[1]:integer, data_pointer[1]:logical, requested_param[1]:string, parent_var[1]:string)

   :param integer index: see interface_scripts/data_handling.ncl.
   :param logical data_pointer: see interface_scripts/data_handling.ncl.
   :param string requested_param: the parameter to fetch for parent_var.
   :param string parent_var: name of the variable for which the requested params belong.

   Return value
      The field corresponding to the variable var.
  
   Description
      Extracts additional data from the climo-file.
      The parameters for irregular grids (e.g. area,lat,lon) are written to
      climo-files by the reformat routines.
      The functionality is similar to "extract_data", but reduced (e.g.
      fetching of var sub-sets is not implemented here).
      This new function was written to avoid modifying "extract_data", which is
      designed to extract the target variable only.
  
   Caveats
  
   References
  
   Modification history:
      2013????-A_gott_kl: written.
  
.. function:: gridcell_area(deltax[1]: numeric, lat_lo[1]: numeric, lat_hi[1]: numeric)

   :param  numeric deltax: longitude resolution [deg].
   :param  numeric lat_lo: lower limit of the box [deg].
   :param  numeric lat_hi: upper limit of the box [deg].

   Return value
      The area of the element in units of [m^2].
  
   Description
      Calculates the area of a grid cell on the sphere.
  
   Modification history:
      20121211-A_righ_ma: written.
  
.. function:: map_area(lat[*]:numeric, lon[*]:numeric)

   :param numeric lat: the latitude coordinate of the map [deg].
   :param numeric lon: the longitude coordinate of the map [deg].

   Return value
      A 2D lat-lon area with the area of each gridbox in units of [m^2].
  
   Description
      Calculates the area of each grid cell on a global map.
  
   Caveats
      Assumes a constant resolution in longitude.
  
   Modification history
      20140819-A_righ_ma: modified to support non-global input.
      20121211-A_righ_ma: written.
  
.. function:: area_operations(field:numeric, latmin[1]:numeric, latmax[1]:numeric, lonmin[1]:numeric, lonmax[1]:numeric, opt[1]:string, l_wgt[1]:logical)

   :param numeric field: a numeric array of rank at least 2; second-to-last and last. dimension must be lat and lon, respectively.
   :param numeric latmin: minimum latitude boundary of the region to be selected.
   :param numeric latmax: maximum latitude boundary of the region to be selected.
   :param numeric lonmin: minimum longitude boundary of the region to be selected.
   :param numeric lonmax: maximum longitude boundary of the region to be selected.
   :param string opt: type of operation: "extract": extracts selected region. "average": averages over the selected region. "sum": integrate over the selected region.
   :param logical l_wgt: if True, calculates area-weighted average/sum (has no effect for opt = "extract").

   Return value
      An array of the same rank as field, of rank-1 or of rank-2, depending on
      opt and on the region boundaries.
  
   Description
      Extracts a selected region or point on a global map.
      Performs the (weighted) average over a selected region on a global map.
      Performs the (weighted) sum over a selected region on a global map.
  
   Caveats
      This function assumes that the input field is a global map.
      Mind the order for lonmin and lonmax (e.g., 60,120 is different from.
      120,60: the order is meant eastwards).
      To consider the global domain, use -90,90,0,360 as arguments.
      If lonmin is outside the field boundaries, it is assigned to 0.
      If lonmax is outside the field boundaries, it is assigned to max(lon).
      For latmin=latmax and lonmin=lonmax the single gridbox is extracted, no
      average/sum is possible in this case.
  
   Modification history
      20140116-A_righ_ma: written.
  
.. function:: select_region(region:string)

   :param string region: a string specifying the region to be selected.

   Return value
      An array with the region boundary as (latmin, latmax, lonmin, lonmax)
      with the name of the region as a string attribute @name.
  
   Description
      Translates a region specification into lat/lon boundaries and a region
      name as an attribute.
  
   Modification history
      20141205-A_gott_kl: adjusted names to Righi et al. (2015).
      20140410-A_fran_fr: extended to midlat, equatorial and polar regions.
      20140129-A_fran_fr: written.
  
