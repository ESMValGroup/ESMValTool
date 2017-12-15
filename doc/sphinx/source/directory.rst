.. _directory:

Directory structure of the ESMValTool
*************************************

An overview of the directory structure used in the ESMValTool is given in :numref:`tab_direc_struc`.
This section summarizes the underlying principles of the structure.

* Common namelist settings (e.g., models, year ranges, diagnostics) are
  usually stored in one place.
* Less common settings may be hidden deeper in the directory structure (see
  also Section :numref:`config_files`).
* Diagnostic scripts that can be used by namelist entries are also stored in
  one place and it generally possible to combine them in a modular way (e.g.,
  using the output of one routine as input for another).
* Reuse of code is strongly encouraged, i.e., one place for each
  functionality (modularity on the technical level).
* The goal is to centralize functionality in individual functions/procedures
  whenever there is the possibility of reusability. New developers are encouraged
  to consider building on or extending existing routines before introducing new
  ones.
* Routines are sorted into folders according to their functionality. Whenever
  possible, the hierarchy level of routines is reflected by their position in the
  directory structure.

