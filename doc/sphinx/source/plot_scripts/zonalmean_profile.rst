:mod:`zonalmean_profile`
========================
.. function:: zonalmean_profile(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype
      source[*,*]
      source!0 = plev
      source!1 = lat
  
   Return value
      A graphic variable
  
   Description
      Draws a pressure-latitude plot
  
   Caveats
  
   Modification history
      20131210-A_fran_fr: written.
  
