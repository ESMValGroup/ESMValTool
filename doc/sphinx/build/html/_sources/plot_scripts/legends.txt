:mod:`legends`
==============
.. function:: legend_lines(wks_in[1],  source,  varname[1] : string)

   :param integer wks_in: workstation ("graphic" object or default will be used)
   :param integer source: * data to be plotted, either passed directly or via netCDF file * legend styles and strings are passed as attributes data@
   :param  string varname: variable name, needed for netCDF files with multiple variables

       The following attributes of the input data are evaluated:
       source@diag_script: name(s) of the calling script(s) 
       source@colors: vector of colors  
       source@thicks: vector of line thicknesses
       source@dashes: vector of dash styles
       source@annots: vector of label strings
       source@nrow: number of rows (optional)
       source@ncol: number of columns (optional)
  
   Returns:
        wks : Workstation with legend
  
   Description:
      * creates an extra plot with a legend, specified by labels and line styles. 
        It will be saved in the outfile directory and returned as a workstation 
  
   Modification history:
      * 20140326 written (Klaus-Dirk.Gottschaldt@dlr.de), 
                 based on code by F. Frank
  
