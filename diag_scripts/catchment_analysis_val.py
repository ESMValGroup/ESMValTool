#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
;;#############################################################################
;; catchment_analysis_val.py
;; Author: Philipp Sommer (MPI,Germany)
;;         Stefan Hagemann (MPI, Germany)
;;         Alexander Loew (LMU, Germany)
;; EMBRACE project
;;#############################################################################
;; Description
;;    The diagnostics calculates water balance components for different
;;    catchments and compares the results against observations.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;     !!!Be aware to have no spaces in the filenames and paths!!!
;;     Otherwise it might produce errors.
;;
;; Modification history
;;    20151029-A_laue_ax: added output of acknowledgements + processed files
;;                        to log-file
;;    20150717-A_somm_ph: written
;;
;; ############################################################################
"""
import os
import pdb
import sys
import logging
import glob
import inspect

sys.path.append("diag_scripts/aux/catchment_analysis")
from catchment_analysis_tool_val import analysecatchments, logger
from esmval_lib import ESMValProject


def main(project_info):
    """
    main interface routine to ESMValTool

    Parameters
    ----------
    project_info : dict
        dictionary that contains all relevant informations
        it is provided by the ESMValTool launcher
    """

    # extract relevant information from project_info using a generic python
    # interface library
    # logging configuration
    verbosity = project_info['GLOBAL']['verbosity']
    if verbosity == 1:
        logger.setLevel(logging.WARNING)
    elif verbosity == 2:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    # create instance of a wrapper that allows easy access to data
    E = ESMValProject(project_info)

    # get information  for runoff/ET script
    # input file directory, returns a dict
    ifile_dict = E.get_raw_inputfile()

    # A-laue_ax+
    diag_script = E.get_diag_script_name()
    res = E.write_references(diag_script,              # diag script name
                             ["A_somm_ph", "A_hage_st"],      # authors
                             ["A_loew_al"],            # contributors
                             ["D_hagemann13jadvmodelearthsyst"],  # diag_references
                             ["E_duemenil00mpi", "E_weedon14waterresourres"], # obs_references
                             ["P_embrace"],            # proj_references
                             project_info,
                             verbosity,
                             False)
    # A-laue_ax-

    # set up of input files dictionary, rewrite ifile_dict to match for
    # catchment_analysis_tool_val
    # ifiles={<<<MODEL>>>:
    #             {<<<VAR>>>:
    #                 {'file':<file>[, 'unit':<unit>, 'cdo':<cdo>,
    #                                'vname':<vname>, 'reffile':<file>,
    #                                'refvname':<vname>]}}}
    # dictionary containing the model names as keys and a second dictionary as
    # value.
    # This (second) dictionaries contain the variable names (eihter 'runoff' or
    # 'ET') as key and a third dictionary as value.
    # This (third) dictionaries contain the key
    #   - 'file'     for the datafile for the climatological mean
    #   and optionally
    #   - 'unit'     for the unit of the data. If 'unit' is not set, it will be
    #                taken from the nc-file multiplied by 'm^2'
    #   - 'cdo'      for additional cdo-commands (e.g. multiplication for
    #                changing
    #                the units)
    #   - 'vname'    for the name of the variable in the datafile (if not set,
    #                the variable name as used for the key will be used).
    #   - 'reffile'  for the reference file (if not set, defaultreffile will be
    #                used and <<<VAR>>> will be replaced by the variable name
    #                as given in the key
    #   - 'refvname' for the name of the variable in the reference file (if not
    #                set, the variable name as used for the key will be used).
    ifiles = {}
    for model, vlst in ifile_dict.items():
        ifiles[str(model)] = {}
        for var in map(str, vlst.keys()):
            if var == 'evspsbl':
                myvar = 'ET'
                # units are are kg m-2 s-1 and therefore must be multiplied
                # by the amount of seconds in one year
                ifiles[str(model)][myvar] = {'unit': 'mm/a', 'vname': var,
                                             'cdo': '-mulc,86400 -muldpy '}
            elif var == 'mrro':
                myvar = 'runoff'
                # units are are kg m-2 s-1 and therefore must be multiplied
                # by the amount of seconds in one year
                ifiles[str(model)][myvar] = {'unit': 'mm/a', 'vname': var,
                                             'cdo': '-mulc,86400 -muldpy '}
            elif var == 'pr':
                myvar = 'precip'
                # units are are kg m-2 s-1 and therefore must be multiplied
                # by the amount of seconds in one year
                ifiles[str(model)][myvar] = {'unit': 'mm/a', 'vname': var,
                                             'cdo': '-mulc,86400 -muldpy '}
            else:
                # only ET, mrro and precip are supported trough reference
                # values. Therefore raise error
                raise ValueError(
                    "Only mrro (runoff), pr (precipitation) and evspsbl "
                    "(evapotranspration) are supported!")
            # try to get reformated data and if this does not work (because no
            # CMIP5 project, use raw input data
            try:
                ifiles[str(model)][myvar]['file'] = map(str, glob.glob(
                    E.get_clim_model_filenames(var)[model]))
            except ValueError:
                indir = str(ifile_dict[model][var]['directory'])
                indir = msd["dir"]
                if indir[-1] != os.sep:
                    indir += os.sep
                infile = str(ifile_dict[model][var]['files'])

                ifiles[str(model)][myvar]['file'] = indir + infile

            # A-laue_ax+
            for file in ifiles[str(model)][myvar]['file']:
                E.add_to_filelist(file)
            # A-laue_ax-

    POUT = str(os.path.join(E.get_plot_dir(), E.get_diag_script_name()))
    POUT = POUT + os.sep
    if not os.path.exists(POUT):
        os.makedirs(POUT)

    # catchment_dir = data_directory + 'cat/'
    catchment_dir = os.path.dirname(os.path.abspath(inspect.getfile(
        inspect.currentframe()))) + '/'
    # set path to nc-file containing the catchment definition
    pcatchment = os.path.join(catchment_dir, "aux/catchment_analysis", "big_catchments.nc")

    # A-laue_ax+
    E.add_to_filelist(pcatchment)
    # A-laue_ax-

    # ---- File input declaration ----
    # input file (needs to contain grid informations). For the
    # analysis the data will be regridded to pcatchments. Data is expected to
    # be in mm/s (if not, see 'cdo' in ifiles dictionary below)

    # ---- switches ----
    # CALCULATION: switch to calculate climatological means (True) or use
    # already existing files
    CALCULATION = True
    # KEEPTIMMEAN: switch to keep files created during calculation. if true,
    # files like 'tmp_<<<MODEL>>>_<<<var>>>_<<<CATCHMENT>>>.<<<fmt>>>' will be
    # stored in the output directory (POUT) containing the timmean data for
    # the catchments (note that the file format <<<fmt>>> is defined by the
    # input file)
    KEEPTIMMEAN = False
    # PLOT: switch to set plotting. If true, one plot for each model named
    # POUT+<<<MODEL>>>_bias-plot.pdf containg all variables will be generated
    PLOT = True
    # SEPPLOTS: Integer to control diagram style. If 2: relative and absolute
    # variables will be plotted in separate files named
    # POUT+<<<MODEL>>>_rel-bias-plot.pdf and
    # POUT+<<<MODEL>>>_abs-bias-plot.pdf.
    # If 3: they will be plotted into one # figure but separated axes in file
    # POUT+<<<MODEL>>>_sep-bias-plot.pdf, if 5: # they will be plotted into one
    # single axes within file POUT+<<<MODEL>>>_bias-plot.pdf. Multiplication
    # (e.g. 6, 15 or 30) is also possible to make more than one option at the
    # same time
    SEPPLOTS = 3
    # ETREFCALCULATION: If True, create reference file
    # ref_ET_catchments.txt from timmean of precfile and runoff data from
    # ref_runoff_catchments.txt.  If None, use default reference values, save
    # them to defaultreffile and delete it afterwards. If False, use existing
    # file defined by defaultreffile or reffile defined in ifiles (see above).
    ETREFCALCULATION = None

    # formatoptions for plot (ax2 modifies diagrams showing relative values)
    # (for keywords see function xyplotformatter in
    # catchment_analysis_tool_val)
    # be aware that the additional setting of a keyword for the plot of the
    # absolute values may also influence the plot of relative value. To prevent
    # from that, define the option for 'ax2' manually.
    # ---
    # set minimal and maximal bounds of y-axis (Syntax: [ymin, ymax]). 'minmax'
    # will cause y-axis limits with minimum and maximum value and makes the
    # plot symmetric around 0 in case of SEPPLOTS%5 == 0. To let the plotting
    # routine (i.e. pyplot) choose the limits, set ylim to None
    # ---
    # Set yticks: integer or list (Default (i.e. without manipulation by the
    # plotting routine of catchment_analysis_tool_val): None). Defines the
    # y-ticks. If integer i, every i-th tick of the automatically

    # other formatoptions as provided by function
    # catchment_analysis_tool_val.xyplotformatter may be included in the
    # dictionary below

    fmt = {
        myvar: {
            'ylim': 'minmax', 'yticks': None,
            'ax2': {'ylim': 'minmax', 'yticks': None}
               }
        for myvar in ifiles.values()[0].keys()
          }

    # start computation
    analysecatchments(
        project_info,
        POUT=POUT,
        pcatchment=pcatchment,
        ifiles=ifiles,
        CALCULATION=CALCULATION,
        KEEPTIMMEAN=KEEPTIMMEAN,
        PLOT=PLOT,
        SEPPLOTS=SEPPLOTS,
        ETREFCALCULATION=ETREFCALCULATION,
        fmt=fmt,
        # precfile=precfile,
        # defaultreffile=defaultreffile
        )
