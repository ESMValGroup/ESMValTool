:mod:`set_operators`
====================
.. function:: UNIQ(a)

   :param integer a: array to analyse, could be any type.

   Return value
      Index array of first occurence of each new element.
  
   Description
      Useful to get an array of unique elements.
  
   Caveats
  
   Reference
  
   Modification history
      20130419-A_gott_kl: written.
  
.. function:: union(set1[*], set2[*])

   :param integer set1: a one-dimensional array.
   :param integer set2: a one-dimensional array.

   Return value
      An array.
  
   Description
      Returns  the union of two sets
  
   Caveats
  
   References
  
   Modification history
  
.. function:: set_inclusive_OR(set1[*]:float, set2[*]:float)

   :param float set1: a one-dimensional array.
   :param float set2: a one-dimensional array.

   Return value
      An array.
  
   Description:
     Returns sorted inclusive 'or' of two sets.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: inlist(item[1]:string, alloweditems[*]:string)

   :param string item:  string to check for in list.
   :param string alloweditems: list with strings to compare with.

   Return value
      A logical: True if "item" is present in "alloweditems", false otherwise.
  
   Description
      Checks for an intem in a list and returns True (False) if it is present
      (missing)
  
   Caveats
  
   References
  
   Modification history
  
.. function:: intersection(array1[*], array2[*])

   :param integer array1: a one-dimensional array.
   :param integer array2: a one-dimensional array.

   Return value
      array: Intersection array.
  
   Description
      Returns the intersection of array1 and array2 or 'False' if no
      intersection is found.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: is_array_subset(subset_array[*], full_array[*])

   :param integer subset_array:  an array of dimension N.
   :param integer full_array: an array of dimension >= N.

   Return value
      A logical: True if "subset_array" a true subset of "full_array", False
      otherwise.
  
   Description
      Checks if an array is a subset of another array
  
   Caveats
  
   Reference
  
   Modification history
  
.. function:: relative_complement(array1[*], array2[*])

   :param integer array1: an array
   :param integer array2: another array

   Return value
      logical: True if there is a complemnet, in this case the complement
                    itself is attached as '@array'
               False  if there isn't a complement
  
   Description
      Substracts all elements in array1 from array2 (array2 - array1)
      See https://en.wikipedia.org/wiki/Complement_(set_theory)#Relative_complement
      for further details and expected behaviour
  
   Caveats
  
   Reference
  
   Modification history
  
.. function:: set_symmetric_difference(array1[*], array2[*])

   :param integer array1: an array
   :param integer array2: another array

   Return value
      array: The symmetric difference of array1 and array2
  
   Description
      Returns all elements only in array1 or only in array2
      See "https://en.wikipedia.org/wiki/Symmetric_difference"
      for furhter details and expected behaviour
  
   Caveats
  
   Reference
  
   Modification history
  
