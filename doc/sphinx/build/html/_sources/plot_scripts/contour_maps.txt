:mod:`contour_maps`
===================
.. function:: contour_map(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype
      source[*,*]
      source!0 = lat
      source!1 = lon
  
   Return value
      A graphic variable.
  
   Description
      Wrapper for gsn_csm_contour_map.
      Opens default wks, if not provided as argument of type "graphic".
      Defines default resources, which are overridden by argument res.
      Creates plot, according to wks & res.
  
   Caveats
      Selection of defaults for res almost arbitrary
      Please check results of all scripts that use this routine if modifying
      the defaults!
      Input via netCDF not yet implemented
  
   Modification history
      20131104-A_gott_kl: written.
  
.. function:: contour_map_polar(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype
      source[*,*]
      source!0 = lat
      source!1 = lon
  
   Return value
      A graphic variable.
  
   Description
      Wrapper for gsn_csm_contour_map_polar.
      Opens default wks, if not provided as argument of type "graphic".
      Defines default resources, which are overridden by argument res_in.
      Creates plot according to wks & res, unless res_in@gsnDraw is set
        to 'False' in your diag_script (to enable panelling, for example).
  
   Caveats
      Selection of defaults for res almost arbitrary
      Please check results of all scripts that use this routine if modifying
      the defaults!
      Input via netCDF not yet implemented
  
   Modification history
      20140623-A_senf_da: now takes res as attributes of source.
      20131218-A_senf_da: written.
  
.. function:: contour_map_ce(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype
      source[*,*]
      source!0 = lat
      source!1 = lon
  
   Return value
      A graphic variable.
  
   Description
      Wrapper for gsn_csm_contour_map_ce.
      Opens default wks, if not provided as argument of type "graphic".
      Defines default resources, which are overridden by argument res.
      Creates plot, according to wks & res.
  
   Caveats
      Selection of defaults for res almost arbitrary
      Please check results of all scripts that use this routine if modifying
      the defaults!
      Input via netCDF not yet implemented
  
   Modification history
      20140228-A_righ_ma: written.
  
.. function::  add_markers_to_map(wks_in[1], plot[1]: graphic, lat[*]: numeric,  lon[*]: numeric, data[*]: numeric)

   :param integer wks_in: input workstation.
   :param  graphic plot: a graphic object representing a contour plot.
   :param  numeric lat: an array of latitude positions for the marker.
   :param  numeric lon: an array of longitude positions for the marker.
   :param  numeric data: an array of values to be overlaid as colored markers.

   Description
      Overlays markers to an existing map plot, using the same color coding of
      the associated labelbar.
  
   Caveats
  
   Modification history:
      20140214-A_righ_ma: written.
  
