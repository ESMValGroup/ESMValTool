:mod:`aux_plotting`
===================
.. function::  create_legend_lines(labels:string,  styles,  outfile:string, opt)

   :param string labels:  labels of the lines
   :param integer styles: style of the lines in the plot @colors @dashes  -> optional (line dashes @thicks  -> optional (line/marker thickness) @markers -> optional (marker index) @sizes   -> optional (marker size)
   :param string outfile: outfile directory
   :param integer opt: "lines" or "markers" legend

   Description
      Draws an extra plot with a legend, specified by labels and styles.
  
   Caveats
  
   Modification history:
      20150508-A_righ_ma: added lines/markers option.
      20150120-A_gott_kl: remove pre-existing file type suffix
      20140305-A_righ_ma: modified to plot always as epsi format.
      20140219-A_fran_fr: written.
  
.. function:: output_type()
  
   Return value
      A string with the output file type
  
   Description
      Provides a default, if file type is not explicitely specified
  
   Caveats
  
   Modification history
      20131028-A_gott_kl: written.
  
.. function::  copy_VarAtt_sel(var1, var2, sel: string)

   :param integer var1: variable of any type and dimension with source attributes
   :param integer var2: variable of any type and dimension, shall receive selected attributes from var1
   :param  string sel: string (or list of strings) that specify the BEGINNING letters of the attributes to copy

   Return value
      var2 gets additional attributes 
  
   Description
      Intended to copy selected plot ressources for use in a different 
      plot routine that may not allow all ressources from var1. 
      Written for function legend_lines (legends.ncl).
  
   Caveats:
  
   Modification history:
      * 20141227-A_gott_kl written
  
.. function:: panelling(wks, plots : graphic, nvert[1] : integer, nhori[1] : integer, pres_in[1] : logical)

   :param integer wks: graphics object associated with a file.
   :param  graphic plots: graphics variable containing one or more plots.
   :param  integer nvert: Maximum allowed number of plots per page (vertical).
   :param  integer nhori: Maximum allowed number of plots per page (horizontal).
   :param  logical pres_in: if it is a logical variable, attributes are used to override panelling defaults

   Return value
      A string with the output file name
  
   Description
      Writes a graphics file
  
   Caveats
      Might not be fit to be used with non-NCL routines.
  
   Modification history
      20131112-A_gott_kl: written.
  
.. function:: get_plot_dir()
  
   Return value
      A string root path for plots
  
   Description
      Provides a default, if plot_dir is not explicitely specified
  
   Caveats
  
   Modification history
      * 20131104-A_gott_kl: written.
  
.. function:: get_outfile_name(diag_script[1] : string, add_specs[1] : string)

   :param  string diag_script: name of diagnostic script(s)
   :param  string add_specs: string containing specific elements to be added to the file name if not of type string, only variable name is used.

   Return value
      Complete outfile name incl. path, additional specifications and file type
      output_dir + diag_script_base + add_specs + file_type
  
   Description
      Fetches file_type, plot_dir, diag_script_base, output_dir via other
      scripts
      Fetches string with additional elements specified within diag script
  
   Caveats
      diag_script may need to be extended by other contributing scripts
  
   Modification history
      20131204-A_senf_da: generalized naming.
      20131104-A_gott_kl: written.
  
.. function:: get_wks(wks_in, diag_script[1]: string, add_specs[1]: string)

   :param integer wks_in: dummy or graphic object
   :param  string diag_script: name of diagnostic script
   :param  string add_specs: own specificactions to be added to file name, i.e. variable name, etc. - needed for function get_outfile_name

   Return value
      wks: graphic object
  
   Description
      Provides a default wks, if wks_in is not of type "graphic".
      Attribute wks@fullname is used to transfer the output file name, since
      wks@name cuts off the path to the file name.
  
   Caveats
  
   Modification history
      * 20131113-A_gott_kl: written.
  
.. function::  add_markers(wks[1] : graphic, plot[1] : graphic, res_in[1] : logical, xpos_in : numeric, ypos_in : numeric)

   :param  graphic wks: valid workstation, e.g. created by get_wks
   :param  graphic plot: plot identifier, e.g. created by gsn_*
   :param  logical res_in: plot ressources that may override local function defaults
   :param  numeric xpos_in: horizontal marker position(s)
   :param  numeric ypos_in: vertical marker position(s)

   Return value
      Attaches polyline IDs as attributes to plot.
  
   Description:
      Adds markers to an existing plot. If a horizontal (vertical) coordinate
      has only one element, then this position is used for all markers.
  
   Caveats:
  
   Modification history:
      * 20140224-A_gott_kl: written for use with profile_plev.ncl in Emmons.ncl
  
.. function::  horizontal_whiskers(wks[1] : graphic, plot[1] : graphic, res_in[1] : logical, xmin_in : numeric, xmax_in : numeric, ypos_in: numeric)

   :param  graphic wks: valid workstation, e.g. created by get_wks.
   :param  graphic plot: plot identifier, e.g. created by gsn_*.
   :param  logical res_in: plot ressources that may override local function defaults.
   :param  numeric xmin_in: vector of whiskers' left ends (same size as xmax & y).
   :param  numeric xmax_in: vector of whiskers' right ends (same size as xmin & y).
   :param  numeric ypos_in: vector of whiskers' vertical positions (must have same size as xmax & xmin).

   Return value
      Attaches polyline IDs as attributes to plot.
  
   Description
      Creates vectors suitable as input for gsn_add_polyline:
        x = (/xmin1,xmax1,_FillValue,xmin2,xmax2,_FillValue, .../)
        y = (/ypos1,ypos1,_FillValue,ypos2,ypos2,_FillValue, .../)
      The separation by _FillValue results in individual whiskers.
      No whisker is created where xmin, xmax or ypos is missing.
  
   Caveats
  
   References
      www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_add_polyline.shtml
  
   Modification history
      20140224-A_gott_kl: written.
  
