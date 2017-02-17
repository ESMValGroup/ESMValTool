:mod:`style`
============
.. function:: project_styleset(info, flag)

   :param integer info: info array, as defined in ./variable_defs.
   :param integer flag: string to look for in info as first priority.

   Return value
      String of style set, determining
      "./diag_scripts/GENERAL_style/style_" + styleset + ".ncl"
  
   Description
      Determines "./diag_scripts/GENERAL_style/style_" + styleset + ".ncl".
      Considers 4 priority levels.
  
   Caveats
  
   References
  
   Modification history:
      20130419-A_gott_kl: written.
  
.. function:: select_style(modelstyles:string, column:integer, specifiers:string)

   :param string modelstyles: string array, relating models/cases/... with plot attributes (defined in function project_style_*)
   :param integer column: integer determining which column to scan in the above array
   :param string specifiers: array of strings to loop over & look for matches in the first column of modelstyles

   Return value
      Vector (integer or string) with one element for each models@name.
  
   Description
      Builds the vector by looping over models@name.
      1st priority: take value from modelstyles, if first column matches.
      2nd priority: take value from row "unknown", if available in modelstyles.
  
   Caveats
  
   References
  
   Modification history
      20130422-A_gott_kl: written.
      20140320-A_gott_kl: modified to take list of model/case/... specifiers 
                       as input argument rather than always using models@name
  
.. function:: unique_labels_min(prio: string)

   :param  string prio: string vector with attribute names (models@*), ordered by priority for annotation (starting with highest).

   Return value
      A vector (string) with one element for each models@name -> each label
      contains the least possible attribute strings.
  
   Description
      Builds the vector by looping over models@name.
      Adds "_attribute" to non-unique labels, until prio is exhausted.
  
   Caveats
  
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
  
   References
  
   Modification history
      20130422-A_gott_kl: written.
  
.. function:: gsnColorRange(lower: numeric, upper: numeric, step: numeric, center: numeric, color_end: integer, center_color: integer)

   :param  numeric lower: cnMinLevelValF.
   :param  numeric upper: cnMaxLevelValF.
   :param  numeric step: cnLevelSpacingF.
   :param  numeric center: The numerical value the colormap is centered on. For anomalies or trends, it's common to use 0.0, so blue means cold or cooling and red means warm or warming.
   :param  integer color_end: The number of colors in colormap (ex. 97 for BlRe, 253 for BlueRed).
   :param  integer center_color: = Color value on the left of the "center" value (see above).

   Description
      Sets the gsnSpreadXXX resources necessary to correctly span a two-color
      colortable.
  
   Return value
      A logical variable with the start and end colors as attributes.
  
   Caveats
  
   Reference
      http://www.ncl.ucar.edu/Applications/Scripts/contoursym_4.ncl
  
   Modification history:
      20130422-A_gott_kl: written.
  
.. function:: project_style_CMIP5(flag:string)

   :param string flag: string determining the type of array requested.

   Return value
      An array of dimsizes(models@name).
      Definition of plot attributes; Returns arrays of dimsizes(models@name)
        flag: "colors": returns an array of colors
                         (either RGB triples or named colors)
        flag: "dashes": returns an array of dash styles (integer numbers)
        flag: "thicks": returns an array of line thicknesses (numeric)
        flag: "annots": returns an array of annotation strings
        flag: "avgstd": returns an array of flags
                          0 -> (model) takes part in calculation of
                               mean & stddev
                          1 -> (obs/reanalysis) takes not part in calculation
                               of mean & stddev
  
   Description
      Definition of plot attributes for CMIP5.
  
   Caveats
  
   Reference
  
   Modification history:
      * 20130419-A_gott_kl: written.
  
.. function:: project_style_DEFAULT(flag:string)

   :param string flag: string determining the type of array requested.

   Return value
      An array of dimsizes(models@name).
      Definition of plot attributes; Returns arrays of dimsizes(models@name)
        flag: "colors": returns an array of colors
                         (either RGB triples or named colors)
        flag: "dashes": returns an array of dash styles (integer numbers)
        flag: "thicks": returns an array of line thicknesses (numeric)
        flag: "annots": returns an array of annotation strings
        flag: "avgstd": returns an array of flags
                          0 -> (model) takes part in calculation of
                               mean & stddev
                          1 -> (obs/reanalysis) takes not part in calculation
                               of mean & stddev
  
   Description
      Definition of plot attributes for DEFAULT.
  
   Caveats
  
   Reference
  
   Modification history:
      * 20130430-A_gott_kl: written.
  
.. function:: project_style_EMAC(flag:string)

   :param string flag: string determining the type of array requested.

   Return value
      An array of dimsizes(models@name).
      Definition of plot attributes; Returns arrays of dimsizes(models@name)
        flag: "colors": returns an array of colors
                         (either RGB triples or named colors)
        flag: "dashes": returns an array of dash styles (integer numbers)
        flag: "thicks": returns an array of line thicknesses (numeric)
        flag: "annots": returns an array of annotation strings
        flag: "avgstd": returns an array of flags
                          0 -> (model) takes part in calculation of
                               mean & stddev
                          1 -> (obs/reanalysis) takes not part in calculation
                               of mean & stddev
   Description:
      * Definition of plot attributes: type depending on flag
      * consider models@ (/"name", "ensemble", "experiment" or "case_name"/) 
        to assign colors, dashes, thicks
      --> predefined different colors/thicks/dashes for selected EMAC runs  
          (case_name is used to distinguish those selected EMAC simulations)
      --> identical other colors/thicks/dashes for not listed EMAC simulations
          (entry "EMAC" in modelstyles)
      --> identical, but again other colors/thicks/dashes for non EMAC data 
          (entry "unknown" in modelstyles)
      --> modelstyles also contains definitions for some obs/rea/specials
   Modification history:
      * 20140321 complete overhaul towards a combination of CMIP5 and
                 EMAC_default styles (Klaus-Dirk.Gottschaldt@dlr.de)
      * 20140120 changed to plot style for EMAC models (Franziska.Frank@dlr.de)
      * 20130430 written (Klaus-Dirk.Gottschaldt@dlr.de)
  
.. function:: project_style_EMAC_default(flag:string)

   :param string flag: string determining the type of array requested.

   Return value
      An array of dimsizes(models@name), type depending on flag.
        flag: "colors": returns an array of colors
                         (either RGB triples or named colors)
        flag: "dashes": returns an array of dash styles (integer numbers)
        flag: "thicks": returns an array of line thicknesses (numeric)
        flag: "annots": returns an array of annotation strings
        flag: "avgstd": returns an array of flags
                          0 -> (model) takes part in calculation of
                               mean & stddev
                          1 -> (obs/reanalysis) takes not part in calculation
                               of mean & stddev
  
   Description
       * Definition of plot attributes, depending on flag
       * EMAC data are distinguished by default colors, but all get solid line,
       * others get grey, but are distinguished by default dashes 
  
   Caveats
  
   Reference
  
   Modification history
      * 20140320-A_gott_kl: renamed EMAC -> EMAC_default 
      * 20130214-A_gott_kl: written for Emmons_diagnostics.
  
.. function:: project_style_righi15gmd(flag:string)

   :param string flag: string determining the type of array requested.

   Return value
      An array of dimsizes(models@name).
      Definition of plot attributes; Returns arrays of dimsizes(models@name)
        flag: "colors": returns an array of colors
                         (either RGB triples or named colors)
        flag: "dashes": returns an array of dash styles (integer numbers)
        flag: "thicks": returns an array of line thicknesses (numeric)
        flag: "annots": returns an array of annotation strings
        flag: "avgstd": returns an array of flags
                          0 -> (model) takes part in calculation of
                               mean & stddev
                          1 -> (obs/reanalysis) takes not part in calculation
                               of mean & stddev
  
   Description
      Definition of plot attributes for righi15gmd.
  
   Caveats
  
   Reference
  
   Modification history:
      * 20140708-A_gott_kl: Adjusting colors to new trunk version, adding OBS
      * 20131213-A_fran_fr: written.
  
.. function:: project_style(info, flag)

   :param integer info: info array, as defined in ./variable_defs.
   :param integer flag: string determining the type of array requested.

   Return value
      An array of dimsizes(models@name).
      Definition of plot attributes; Returns arrays of dimsizes(models@name)
        flag: "colors": returns an array of colors
                         (either RGB triples or named colors)
        flag: "dashes": returns an array of dash styles (integer numbers)
        flag: "thicks": returns an array of line thicknesses (numeric)
        flag: "annots": returns an array of annotation strings
        flag: "avgstd": returns an array of flags
                          0 -> (model) takes part in calculation of
                               mean & stddev
                          1 -> (obs/reanalysis) takes not part in calculation
                               of mean & stddev
  
   Description
      Definition of plot attributes for righi15gmd.
  
   Caveats
  
   Reference
  
   Modification history:
      * 20130419-A_gott_kl: written.
  
.. function:: place_debuginfo(wks [1] : graphic, debugstring [1] : string, res [1] : logical, plot [1] : graphic)

   :param  graphic wks:       -- current workstation
   :param  string debugstring:   -- string to attach
   :param  logical res:       -- resource settings for display box
   :param  graphic plot:       -- graphic object to draw text onto

   Description:
       Places the text string debugstring onto wks
  
.. function::  place_description(wks [1] : graphic, description [1] : string, y_ndc_coord [1] : float)

   :param  graphic wks:          -- current workstation
   :param  string description:  -- string to attach
   :param  float y_ndc_coord:  -- vertical placement in ndc space (-1 for default)

   Description:
       Places the text strings in array debugboxes onto wks
  
.. function:: format_units(str[1]: string)

   :param  string str: a string.

   Return value
      A string.
  
   Description
      Reformats units string to properly display superscripts
      (e.g. m^2 --> m~S1~2)
  
   Caveats
      Convering only very few cases, to be extended.
  
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
  
