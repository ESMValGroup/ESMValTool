:mod:`taylor_plot`
==================
.. function:: taylor_plot(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype:
      source = (ndiag, nmod, statistic)
      source(:,:,0) = normalized standard deviation
      source(:,:,1) = correlation
      source!0 = diagnostics
      source!1 = models
      source!2 = statistic
  
   Return value:
      A graphic variable.
  
   Caveats
      The taylor plot is drawn in two different ways, depending on ndiag:
        ndiag = 1: models are drawn with different colors and markers, a 
                   separate legend file is created.
        ndiag > 1: variables are drawn with different colors, models are
                   marked with a numbers. Two legends are drawn within the
                   plot frame: one for the variables (markers) and one for the
                   models (numbers).
  
   Modification history:
      20150505-A_righ_ma: written based on the original NCL code.
  
