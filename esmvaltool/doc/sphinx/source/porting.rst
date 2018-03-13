.. _porting:

Porting namelists and diagnostics to ESMValTool v2.0
****************************************************

This guide summarizes the main steps to be taken in order to port an ESMValTool namelist and the corresponding diagnostics from v1.0 to v2.0. ESMValTool v2.0 is being developed in the public git branch REFACTORING_backend. It is strongly recommended to create a branch from REFACTORING_backend for each of the namelist to be ported and to name it REFACTORING_<namelist>. 

REFACTORING_backend contains both v1.0 and v2.0, the latter within the esmvaltool/ directory. It is therefore possible, and recommended, to run both versions of the ESMValTool within the same branch: this will facilitate testing and comparison of the two version as long as the porting process proceeds.
