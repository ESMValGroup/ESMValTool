:mod:`portrait_plot`
====================
.. function:: portrait_plot(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype:
      source = (ndiag, nmod)
      source!0 = diagnostics
      source!1 = models
  
      source = (ndiag, nmod, nref)
      source(:,:,0) = reference model
      source(:,:,1) = alternative model (optional)
      source!0 = diagnostics
      source!1 = models
  
   Return value:
      A graphic variable.
  
   Caveats
  
   Modification history:
      20140605-A_righ_ma: modified with flexible plot shapes.
      20140204-A_fran_fr: extended.
      20140114-A_righ_ma: written.
  
