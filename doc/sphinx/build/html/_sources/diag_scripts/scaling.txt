:mod:`scaling`
==============
.. function:: convert_units(var:numeric, units_to[1]:string)

   :param numeric var: a numeric field on any dimensionality, it must have a units attribute.
   :param string units_to: a string, specifying the units to be converted to.

   Return value
      A numeric field of the same dimensionality of var.
  
   Description
      Converts units of var from the one specified by the units attribute to
      the one given by the units_to argument. An error message is issued if
      the requested unit conversion is not available.
  
   Caveats
      Covers only a few cases, to be extended.
  
   References
  
   Modification history:
      20150216-A_righ_ma: written.
  
.. function:: scalfac(u:numeric, digits[1]:numeric)

   :param numeric u: numeric field of any dimension containing the data to be scaled.
   :param numeric digits: number of desired relevant digits left of the decimal point.

   Return value
      Scaling factor of type float, which would need to be multiplied to the
      original data to bring them in the desired range (needed e.g. for
      annotation).
  
   Description
      Calculates a factor for scaling the input data to a range with the
      desired number of relevant digits left of the decimal point.
      The sign is counted as one of those digits if required for negative
      values.
  
   Caveats
  
   References
  
   Modification history:
      20140220-A_gott_kl: written.
  
.. function:: scale_units(u:numeric, digits[1]:numeric)

   :param numeric u: numeric field of any dimension containing the data to be scaled.
   :param numeric digits: number of desired relevant digits left of the decimal point (needed e.g. for annotation).

   Return value
      Scaling factor of type float, which would need to be multiplied to the
      original data to bring them in the desired range.
      The string of the new unit is added as attribute.
      If the "units" string contains function codes, the code used will be
      added in s@FuncCode.
  
   Description
      Calculates a factor for scaling the input data to a range with the
      desired number of relevant digits left of the decimal point.
      The result is further scaled to match a "nice" unit (if available)
      depending on the original unit attribute of u.
      Format codes: http://www.ncl.ucar.edu/Applications/text.shtml
  
   Caveats
  
   References
  
   Modification history
      20140220-A_gott_kl: written.
  
.. function:: rescale(u:numeric, s[1]:numeric)

   :param numeric u: numeric field of any dimension containing the data to be scaled.
   :param numeric s: scaling factor, e.g. as output by function scale_units.

   Return value
      Field of the same dimensions and type of u, rescaled to the units
      of s.
  
   Description
      Rescales the input field according to scaling factor s.
      Metadata of u are kept, but the units attribute will be changed according
      to s@units.
  
   Caveats
  
   References
  
   Modification history
      20140221-A_gott_kl: written.
  
