.. _porting:

Porting namelists and diagnostics to ESMValTool v2.0
****************************************************

This guide summarizes the main steps to be taken in order to port an ESMValTool namelist and the corresponding diagnostics from v1.0 to v2.0. ESMValTool v2.0 is being developed in the public git branch REFACTORING_backend. It is strongly recommended to create a branch from REFACTORING_backend for each of the namelist to be ported and to name it REFACTORING_<namelist>. 

REFACTORING_backend contains both v1.0 and v2.0, the latter within the esmvaltool/ directory. It is therefore possible, and recommended, to run both versions of the ESMValTool within the same branch: this will facilitate testing and comparison of the two version as long as the porting process proceeds.

Converting xml to yml
=====================

In ESMValTool v2.0, the main namelist (hereafter *"the namelist"*) is written in yaml format (`Yet Another Markup Language format <http://www.yaml.org/>`_). It may be useful to activate syntax highlighting for the editor in use. This improves the readability of the namelist file and facilitates the editing, especially concerning the indentations which are essential in the yml format (like in python). Instructions can be easily found online, for example for `emacs <https://www.emacswiki.org/emacs/YamlMode>`_ and `vim <http://www.vim.org/scripts/script.php?script_id=739>`_.

For each of the ESMValTool v1.0 namelists, a very first draft in yml format has already been created and is available in ./nml/. This can be used as a starting point, but keeping in mind that it has been created at a very early stage of v2.0 development and it will certainly need further changes.

The namelist file in yml shall first be moved to the esmvaltool/ directory where the v2.0 is being developed::

        mv ./nml/<namelist>.yml ./esmvaltool/namelists/<namelist>.yml


This will help to keep track of which namelists have been already ported to the new version.

The yaml namelist can now be edited and tested, starting with a few models and one diagnostics and proceed gradually. A corresponding version of the v1.0 namelist (xml) shall also be developed in parallel in order to compare the output of the two versions and to make sure that the porting has been successful (see "Testing" below).
