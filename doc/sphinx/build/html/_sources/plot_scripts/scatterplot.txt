:mod:`scatterplot`
==================
.. function:: scatterplot(wks_in[1], source, varname[1]: string, reflines: logical, stats: logical)

   :param integer wks_in: workstation ("graphic" object or default will be used).
   :param integer source: data to be plotted (see source prototype above) or a NetCDF filename with data.
   :param  string varname: variable name in the file. logmode: if true, log scale will be used on both axes.
   :param  logical reflines: show/hide reference lines (1:1 and +/- factor 2).
   :param  logical stats: show/hide summary of statistical values on the plot.

   Source prototype
      source = (2, npoints)
      source(0, :) = x-axis values
      source(1, :) = y-axis values
      source!0 = models
  
   Return value
      A graphic object.
  
   Description
      Ceates a scatter plot and optionally adds lines for the factor-2 range
      and statistical information.
  
   Caveats
      Linear axes are used as default, log axis can be optionally provided by
      the calling diag_script.
  
   Modification history
      20140228-A_righ_ma: written.
  
