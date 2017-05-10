:mod:`style`
============
.. function:: unique_labels_min(prio: string)

   :param  string prio: string vector with attribute names (models@*), ordered by priority for annotation (starting with highest).

   Return value
      A vector (string) with one element for each models@name -> each label
      contains the least possible attribute strings.
  
   Description
      Builds the vector by looping over models@name.
      Adds "_attribute" to non-unique labels, until prio is exhausted.
  
   Caveats
       Uses models@*, which is available here anyway.
  
   References
  
   Modification history
      20130422-A_gott_kl: written.
  
.. function:: unique_labels_all(prio:string)

   :param string prio: string vector with attribute names (models@*), ordered by priority for annotation (starting with highest)

   Description
      Builds the vector by looping over models@name.
      Adds "_attribute" until prio is exhausted or until labels are unique.
  
   Return value
      A vector (string) with one element for each models@name -> all labels
      contain the same (least possible) number of attribute strings.
  
   Caveats
      Uses models@*, which is available here anyway.
  
   References
  
   Modification history
      20130422-A_gott_kl: written.
  
.. function:: project_style(info, flag)

   :param integer info: info array, as defined in ./variable_defs.
   :param integer flag: string determining the type of array requested: "annots": annotation strings. "colors": colors (named colors, RGB or RGBA codes) "dashes": line dash patterns. "thicks": line thicknesses. "markers": marker indexes. ;,            "avgstd": average/standard deviation flags (0 = takes part in the calculation of mean and standard deviation, 1 = does not take part; usually 0 is for models and 1 for observations and reanalyses).

   Return value
      An array of the same size of models@name, with the stlye information for
      the given flag. The type depends on the flag.
  
   Description
      Retruns style informations (annotations, colors, line dash patterns, line
      thicknesses, marker indexes and avgstd flat) based on a given styleset.
      The styleset is determined based on the following priority list:
        1st: style information for the given flag explicitely set as
             diag_script_info@$flag$
        2nd: styleset explicitely set as diag_script_info@styleset
        3rd: styleset not defined, set to DEFAULT
  
   Caveats
  
   References
  
   Modification history
      20150512-A_righ_ma: modified to read style info from external style
                          files, instead of using hard-coded values in the
                          code. Functionalities of the project_styleset and
                          project_style_<styleset> functions porteed here.
      20130419-A_gott_kl: written.
  
.. function:: project_style_GO(flag:string)

   :param string flag: = string determining the type of array requested Return value: array of dimsizes(models@name) * Definition of plot attributes; Returns arrays of dimsizes(models@name) * flag = "colors": returns an array of colors (either RGB triples or named colors) * flag = "dashes": returns an array of dash styles (integer numbers) * flag = "thicks": returns an array of line thicknesses (numeric) * flag = "annots": returns an array of annotation strings * flag = "avgstd": returns an array of flags 0 -> (model) takes part in calculation of mean & stddev 1 -> (obs/reanalysis) takes not part in calculation of mean & stddev Description: * Definition of plot attributes: type depending on flag Modification history: * 20130419 written (Klaus-Dirk.Gottschaldt@dlr.de) ; local result, modelstyles, flag begin verbosity  = stringtointeger(getenv("ESMValTool_verbosity")) info_output("<<<<<<<< Entering style_GO.ncl", verbosity, 8)

.. function:: place_debuginfo(wks[1]:graphic, debugstring[1]:string, res[1]:logical, plot[1]:graphic)

   :param graphic wks: current workstation.
   :param string debugstring: string to attach.
   :param logical res: resource settings for display box.
   :param graphic plot: graphic object to draw text onto

   Return value
  
   Description
      Places the text string debugstring onto wks.
  
   Caveats
  
   References
  
   Modification history
  
.. function::  place_description(wks[1]:graphic, description[1]:string, y_ndc_coord[1]:float)

   :param graphic wks: current workstation
   :param string description: string to attach
   :param float y_ndc_coord: vertical placement in ndc space (-1 for default)

   Return value
  
   Description
       Places the text strings in array debugboxes onto wks
  
   Caveats
  
   References
  
   Modification history
  
.. function:: gsnColorRange(lower:numeric, upper:numeric, step:numeric, center:numeric, color_end:integer, center_color:integer)

   :param numeric lower: cnMinLevelValF.
   :param numeric upper: cnMaxLevelValF.
   :param numeric step: cnLevelSpacingF.
   :param numeric center: The numerical value the colormap is centered on. For anomalies or trends, it's common to use 0.0, so blue means cold or cooling and red means warm or warming.
   :param integer color_end: The number of colors in colormap (ex. 97 for BlRe, 253 for BlueRed).
   :param integer center_color: = Color value on the left of the "center" value (see above).

   Description
      Sets the gsnSpreadXXX resources necessary to correctly span a two-color
      colortable.
  
   Return value
      A logical variable with the start and end colors as attributes.
  
   Caveats
  
   References
      http://www.ncl.ucar.edu/Applications/Scripts/contoursym_4.ncl
  
   Modification history
      20130422-A_gott_kl: written.
  
.. function:: format_units(str[1]: string)

   :param  string str: a string.

   Return value
      A string.
  
   Description
      Reformats units string to properly display superscripts
      (e.g. m^2 --> m~S1~2)
  
   Caveats
      Currently convering only very few cases, to be extended.
  
   References
  
   Modification history
      20140320-A_righ_ma: written.
  
.. function:: set_log_ticks(fmin[1]:numeric, fmax[1]:numeric, opt[1]:string)

   :param numeric fmin: the minimum axis value (as specified by trXMinF or trYMinF)
   :param numeric fmax: the maximum axis value (as specified by trXMaxF or trYMaxF)
   :param string opt: "major" or "minor", to return major or minor ticks.

   Return value
      An array of tickmark values, to be used in trXBValues, trYLValues and
      trXBLabels, trYLLabels
  
   Description
      Since NCL only allows for exponential notation in plots with logarithmic
      axes, this function creates an explicit sets of tickmarks with float
      notation.
  
   Caveats
  
   References
  
   Modification history
      20141003-A_righ_ma: written.
  
.. function:: sort_alphabetically(orig_names[*], idx_exclude, dest_exclude)

   :param integer orig_names: the array of model names prior to sorting
   :param integer idx_exclude: the index(es) to be excluded from sorting, -1 to include everything
   :param integer dest_exclude: the position where to put the excluded values after sorting ("begin" or "end")

   Return value
      An integer array of the sime size of orig_names, with the permutation
      index to be used to sort the array in alphabetical order.
  
   Description
      Given an array of model names, this function returns the permutation
      indexes which can be used to sort the array in alphabetical order.
      Certain elements of the array can be excluded from the sorting and
      placed either at the beginning or at the end of the sorted array (e.g.,
      for sorting model alphabetically but leaving observations at the end,
      or multi-model mean at the beginning).
      The function itself does NOT perform any sorting, it just returns the
      permutation indexes. These have to be applied to both the data AND
      the model coordinate to get consistent results.
      For example:
  
          data(models|:, lat|:, lon|:)
          pid = sort_alphabetically(data&models, -1, "")
          sorted_data = data(pid, :, :)
          sorted_data&models = data&models(pid)
  
   Caveats
      Overwriting the original data can lead to incorrect results:
          data = data(pid, :, :)          ; THIS IS WRONG!
          data&models = data&models(pid)  ; THIS IS WRONG!
  
   References
  
   Modification history
      20151028-A_righ:ma: written.
  
