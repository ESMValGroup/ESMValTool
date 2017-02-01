#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
;;#############################################################################
;; GENERAL ROUTINES FOR PYTHON CATCHMENT ANALYSIS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    def analysecatchments
;;    def computemean
;;    def writefile
;;    def readfile
;;    def check
;;    def plotfile
;;    def make_barplot
;;    def make_coefficient_plots
;;    def xyplotformatter
;;    def round_to_05
;;    def make_reference_calculations
;;    def runoff_ref_calculation
;;    def ET_ref_calculation
;;    def precip_ref_calculation
;;    def setdefaultfromdict
;;    class defaults
;;#############################################################################
"""
import datetime as dt
import glob
import logging
import numpy as np
import os
import pdb
import re
import subprocess as spr
import sys

from copy import deepcopy
from difflib import get_close_matches
from itertools import chain, cycle

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter, ScalarFormatter, MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import netCDF4 as nc
from cdo import Cdo
cdo = Cdo()

__author__ = "Philipp Sommer (philipp.sommer@mpimet.mpg.de)"

logger = logging.getLogger(__name__)
formatter = logging.Formatter(
    "[%(name)s] - %(levelname)s: %(message)s")

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


def analysecatchments(project_info,
                      POUT,
                      pcatchment,
                      ifiles,
                      CALCULATION,
                      KEEPTIMMEAN,
                      PLOT,
                      SEPPLOTS,
                      precfile=None,
                      MRROREFCALCULATION=None,
                      ETREFCALCULATION=None,
                      PRECIPREFCALCULATIN=None,
                      catchments=None,
                      defaultreffile=None,
                      fmt={}):
    """
    Analyse catchments. Currently supported variables: 'runoff' and 'ET'
    INPUT:
        - project_info  Configuration dictionary from input XML-file
        - POUT          output path for results
        - pcatchment    path to nc-file containing the catchment definition
        - ifiles        dictionary. Syntax:
                        ifiles={<<<MODEL>>>:
                                        {<<<VAR>>>:
                                                {'file': <file>[,
                                                 'unit': <unit>,
                                                 'cdo': <cdo>,
                                                 'vname': <vname>,
                                                 'reffile': <file>,
                                                 'refvname': <vname>]}}}
                        dictionary containing the model names as keys and a
                        second dictionary as value.
                        This (second) dictionaries contain the variable names
                        (eihter 'runoff' or 'ET') as key and a third dictionary
                        as value.
                        This (third) dictionaries contain the key
                        - 'file'     for the datafile for the climatological
                                     mean and optionally
                        - 'unit'     for the unit of the data. If 'unit' is not
                                     set, it will be taken from the nc-file
                                     multiplied by 'm^2'
                        - 'cdo'      for additional cdo-commands (e.g.
                                     multiplication for changing the units)
                        - 'vname'    for the name of the variable in the
                                     datafile (if not set, the variable name
                                     as used for the key will be used).
                        - 'reffile'  for the reference file (if not set,
                                     defaultreffile will be used and
                                     <<<VAR>>> will be replaced by the variable
                                     name as given in the key
                        - 'refvname' for the name of the variable in the
                                     reference file (if not set, the
                                     variable name as used for the key will be
                                     used).

   ---- switches ----
        - CALCULATION       True/False. calculate climatological means (True)
                            or use already existing files
        - KEEPTIMMEAN       True/False. If True, files like
                         'tmp_<<<MODEL>>>_<<<var>>>_<<<CATCHMENT>>>.<<<fmt>>>'
                            will be stored in POUT containing the timmean data
                            for the catchments (note that the file format
                            <<<fmt>>> is
                            defined by the input file)
        - PLOT              True/False. If True, plot according
                            SEPPLOTS-options all variables into one file
        - SEPPLOTS          Integer. If 2: relative and absolute variables will
                            be plotted in separate files named
                            POUT+<<<MODEL>>>_rel-bias-plot.pdf and
                            POUT+<<<MODEL>>>_abs-bias-plot.pdf. If 3: they will
                            be plotted into one figure but separated axes in
                            file POUT+<<<MODEL>>>_sep-bias-plot.pdf, if 5: they
                            will be plotted into one single axis within file
                            POUT+<<<MODEL>>>_bias-plot.pdf. Multiplication
                            (e.g. 6, 15 or 30) is also possible to make more
                            than one option at the same time
        - precfile          file with precipitation data stored in variable
                            'Rainf_GPCC' with units in mm/a for the calculation
                            of ET reference (only needed if ETREFCALCULATION is
                            True).
        - MRROREFCALCULATION  True/False/None. If True, create reference file
                            ref_runoff_catchments.txt for runoff using data
                            from ref_runoff_catchments_old.txt and
                            catchment-areas.txt.
                            If None, use default reference values, save them to
                            defaultreffile and delete it afterwards. If False,
                            use existing file defined by defaultreffile or
                            reffile defined in ifiles (see above).
        - ETREFCALCULATION  True/False/None. If True, create reference file
                            ref_ET_catchments.txt from timmean of precfile and
                            runoff data from ref_runoff_catchments.txt.  If
                            None, use default reference values, save them to
                            defaultreffile and delete it afterwards. If False,
                            use existing file defined by defaultreffile or
                            reffile defined in ifiles (see above).
        - catchments        dictionary. Catchments with name as used in REFFILE
                            as key and the catchment number as used in
                            pcatchment as value. If None, use default
                            catchments and numbers.
                            Syntax: catchments = {catchment:catchmentnumber}
        - defaultreffile    txt-file with reference numbers (Observations) for
                            each catchment. <<<VAR>>> will be replaced by the
                            specific variable short name (may be individually
                            specified for each variable in ifiles). Reference
                            files can either be a .txt-file (be aware of the
                            ending!) containing the reference value for each
                            catchment, or be global a grid cell based file
                            (e.g. in nc-format). If None, the defaultreffile
                            will be assumed to be like
                            ref_<<<VAR>>>_catchments.txt stored in POUT.

    --- Plotting formatoptions ---
        - fmt           dictionary defining formatoptions for the plot.
                        Possible keywords are given in xyplotformatter
    OUTPUT:
    txt-Files containing catchment data:
        - POUT + <<<MODEL>>>_<<<VAR>>>_catchments.txt: absolute value for each
            catchment of the model <<<MODEL>>> and variable <<<VAR>>>
        - POUT + <<<MODEL>>>_bias_<<<VAR>>>_catchments.txt: differences
            (absolute and relative) to reference each catchment of the model
            <<<MODEL>>> and variable <<<VAR>>> timmean-files (if KEEPTIMMEAN
            is True)
        - POUT/tmp_<<<MODEL>>>_<<<VAR>>>_<<<CATCHMENT>>>.xxx: timmean-files
            containing the timmean-values.
    plotfiles (see switches PLOT and SEPPLOTS above)
        - POUT+<<<MODEL>>>_bias-plot.pdf
        - POUT+bias-plot.pdf
    """

    # set up default catchments
    if catchments is None:
        catchments = defaults.catchments

    # modify file names if different grids are set
    # perfom reference calculations
    if defaultreffile is None:
        defaultreffile = POUT + defaults.defaultreffile
    PREF = os.path.dirname(defaultreffile) + os.path.sep
    make_reference_calculations(defaultreffile,
                                PREF,
                                MRROREFCALCULATION,
                                ETREFCALCULATION,
                                PRECIPREFCALCULATIN)

    # first check whether everything is okay
    check(catchments, defaultreffile, ifiles)

    # create data an save it to dictionary
    if CALCULATION:
        data = computemean(ifiles, pcatchment, catchments, POUT)
        diff2ref = {model: {var: {} for var in sorted(data[model].keys())}
                    for model in sorted(data.keys())}
        fullrefdata = {model: {var: {} for var in sorted(data[model].keys())}
                       for model in sorted(data.keys())}
        # perform subtraction and save data for each model and variable
        for model in sorted(data.keys()):
            for var in sorted(data[model].keys()):
                # read reference data
                REFFILE = ifiles[model][var].get('reffile', defaultreffile)
                # if reffile is txt-file, use the method readfile to calculate
                # the mean
                if REFFILE[-4:] == '.txt' or type(REFFILE) is dict:
                    refdata = readfile(REFFILE.replace('<<<VAR>>>', var))
                    fullrefdata[model][var] = refdata.copy()
                    # check wheter units match
                    for catchment in sorted(data[model][var].keys()):
                        if (data[model][var][catchment]['unit'] !=
                                refdata[catchment]['unit']):
                                logger.warning(
                                    "Attention! Unit of reference file is "
                                    "%s in catchment %s whereas unit for model"
                                    " %s is %s for variable %s. Be sure that "
                                    "no errors are produced!",
                                    refdata[catchment]['unit'],
                                    catchment,
                                    model,
                                    data[model][var][catchment]['unit'],
                                    var)

                        diff2ref[model][var][catchment] = {
                            'data': data[model][var][catchment]['data'] -
                            refdata[catchment]['data'],
                            'unit': data[model][var][catchment]['unit']}
                        diff2ref[model][var][catchment]['rel'] = \
                            diff2ref[model][var][catchment]['data'] / \
                            abs(refdata[catchment]['data']) * 100.

                # if reffile has full data (e.g. in nc-format) we use the
                # timmeanfile and cdos
                else:
                    # compute reference data for completeness
                    refdata = computemean(
                        {'reference': {
                            var: {'file': REFFILE,
                                  'unit': ifiles[model][var]['unit'],
                                  'vname': ifiles[model][var].get(
                                      'refvname', var)}}},
                        pcatchment, catchments, POUT,
                        KEEPTIMMEAN=False)
                    writefile(refdata['reference'][var], POUT +
                              'ref_' + var + '_catchments.txt',
                              'Observations:')

                    # compute difference to reference data for each grid cell
                    # and perform fldmean
                    diff2ref[model][var] = {catchment: {
                        'data': float(
                            cdo.output(input=(
                                "-fldmean -sub -selname,%s "
                                "%s -selname,%s %s") % (
                                    ifiles[model][var].get('vname', var),
                                    data[model][var][catchment]['timmeanfile'],
                                    ifiles[model][var].get('refvname', var),
                                    REFFILE))[0]),
                        'unit': data[model][var][catchment]['unit']
                        } for catchment in sorted(data[model][var].keys())}

                    # compute relative values
                    for catchment in sorted(diff2ref[model][var].keys()):
                        diff2ref[model][var][catchment]['rel'] = \
                            diff2ref[model][var][catchment]['data'] / \
                            abs(refdata[catchment]['data'] * 100.)
                if not KEEPTIMMEAN:
                    spr.call(
                        ['rm'] +
                        [data[model][var][catchment]['timmeanfile'] for
                            catchment in sorted(data[model][var].keys())])
                # # write data to file for each catchment
                writefile(data[model][var], POUT + model + '_' + var +
                          '_catchments.txt', ' ' + model + ':')
                # write bias data to file for each catchment
                writefile(diff2ref[model][var], POUT + model + '_bias_' + var
                          + '_catchments.txt', ' ' + model + '-bias:')
    else:
        data = {}
        diff2ref = {}
        fullrefdata = {}
        for model, model_dict in ifiles.items():
            data[model] = {}
            diff2ref[model] = {}
            fullrefdata = {model: {}}
            for var, var_dict in model_dict.items():
                data[model][var] = readfile(
                    POUT + model + '_' + var + '_catchments.txt')
                diff2ref[model][var] = readfile(
                    POUT + model + '_bias_' + var + '_catchments.txt')
                REFFILE = ifiles[model][var].get('reffile', defaultreffile)
                fullrefdata[model][var] = readfile(REFFILE.replace('<<<VAR>>>',
                                                                   var))

    if MRROREFCALCULATION is None:
        os.remove(defaultreffile.replace('<<<VAR>>>', 'runoff'))
    if ETREFCALCULATION is None:
        os.remove(defaultreffile.replace('<<<VAR>>>', 'ET'))
    if PRECIPREFCALCULATIN is None:
        os.remove(defaultreffile.replace('<<<VAR>>>', 'precip'))

    # ------------ plotting ------------------------
    if PLOT:
        suffix = project_info['GLOBAL']['output_file_type'].lower()
        if SEPPLOTS % 2 == 0:
            plotfile(diff2ref,
                     os.path.join(POUT, '<<<MODEL>>>_rel-bias-plot.' + suffix),
                     dataorrel='rel',
                     fmt=fmt,
                     output_format=suffix)

            plotfile(diff2ref,
                     os.path.join(POUT, '<<<MODEL>>>_abs-bias-plot.' + suffix),
                     dataorrel='data',
                     fmt=fmt,
                     output_format=suffix)

        if SEPPLOTS % 3 == 0:
            plotfile(diff2ref,
                     os.path.join(POUT, '<<<MODEL>>>_sep-bias-plot.' + suffix),
                     fmt=fmt,
                     output_format=suffix)

        if SEPPLOTS % 5 == 0:
            plotfile(diff2ref,
                     os.path.join(POUT, '<<<MODEL>>>_bias-plot.' + suffix),
                     sep=False,
                     fmt=fmt,
                     output_format=suffix)

        make_coefficient_plots(data,
                               diff2ref,
                               fullrefdata,
                               os.path.join(POUT, 'ro-et_coefficient_biases.' + suffix),
                               os.path.join(POUT, 'ro_coefficient-rel-pr_biases.' + suffix),
                               output_format=suffix)


def computemean(ifiles, pcatchment, catchments, POUT="", KEEPTIMMEAN=True):
    """compute climatological mean of input files
    Input:
        ifiles              dictionary. Syntax:
                        ifiles={<<<MODEL>>>:
                                    {<<<VAR>>>:
                                            {'file':<file>[, 'unit':<unit>,
                                                            'cdo':<cdo>,
                                                            'vname':<vname>,
                                                            'reffile':<file>,
                                                            'refvname':
                                                                <vname>]}}}
        pcatchment  string. Path to nc-file containing the catchment
                      definition
        catchments  {<<<CATCHMENT>>>:<<<CATCHMENTNUMBER>>>}. Catchments with
                      name as used in REFFILE as key and the catchment number
                      as used in pcatchment as value
        POUT        string. Output path for TIMMEAN files
        KEEPTIMMEAN True/False. Keep the timmeanfiles
  """
    logger.info('Computing means...')
    # set up dictionary containing the names of the maskfiles
    maskfiles = {catchment: POUT + 'tmp_' + catchment + '.nc' for catchment in
                 sorted(catchments.keys())}

    # create grid file for remapping
    gridfile = POUT + 'tmp-grid_' + str(np.random.randint(0, 100000)) + '.txt'
    with open(gridfile, 'w') as f:
        f.writelines([name + '\n' for name in cdo.griddes(input=pcatchment)])

    # create mask files
    for catchment in sorted(maskfiles.keys()):
        cdo.eqc(catchments[catchment], input=pcatchment,
                output=maskfiles[catchment])
    # compute climatological means and save it (together with the unit) to
    # output
    output = {model: {
        var: {catchment: {} for catchment in sorted(maskfiles.keys())}
        for var in sorted(ifiles[model].keys())}
        for model in sorted(ifiles.keys())}
    TMPFILES = []
    for model in sorted(ifiles.keys()):
        for var in sorted(ifiles[model].keys()):
            modelfiles = deepcopy(ifiles[model][var]['file'])
            if isinstance(modelfiles, (str, unicode)):
                modelfiles = [modelfiles]
                ending = modelfiles[0].split('.')[-1]
                thisremappedfiles = ["%stmp-remap_%i.%s" % (
                    POUT, np.random.randint(0, 100000),
                    ending)]
                TMPFILES += thisremappedfiles
            else:
                ending = modelfiles[0].split('.')[-1]
                thisremappedfiles = ["%stmp-remap_%i.%s" % (
                    POUT, np.random.randint(0, 100000), ending)
                    for FILE in modelfiles]
                TMPFILES += thisremappedfiles
            for FILE in modelfiles:
                remappedfile = thisremappedfiles[modelfiles.index(FILE)]
                cdo.remapcon(gridfile,
                             input=' -selname,%s %s' % (
                                 ifiles[model][var].get('vname', var),
                                 FILE),
                             output=remappedfile)
            for catchment in sorted(maskfiles.keys()):
                # define temporary (if not KEEPTIMMEAN) file for
                # timmean calculation
                timmeanfile = POUT + 'tmp_' + \
                    '_'.join(name for name in [model, var, catchment]) + \
                    '.' + ending
                # calculate timmean
                if len(thisremappedfiles) == 1:
                    cdo.timmean(
                        input='-ifthen %s %s' % (maskfiles[catchment],
                                                 thisremappedfiles[-1]),
                        output=timmeanfile)
                else:
                    TMPFILES.append(
                        "%stmp-remap_%i.%s" % (
                            POUT, np.random.randint(0, 100000), ending))
                    cdo.mergetime(
                        input='-timmean %s' % (
                            ' -timmean '.join(thisremappedfiles)),
                        output=TMPFILES[-1])
                    cdo.timmean(
                        input='-ifthen %s %s' % (maskfiles[catchment],
                                                 TMPFILES[-1]),
                        output=timmeanfile)
                output[model][var][catchment]['timmeanfile'] = timmeanfile
                output[model][var][catchment]['data'] = float(
                    cdo.output(input=ifiles[model][var].get('cdo', '') +
                               '-fldmean ' + timmeanfile)[0])
                if ifiles[model][var].get('unit', False):
                    output[model][var][catchment]['unit'] = \
                        ifiles[model][var]['unit']
                elif output[model][var][catchment]['unit'][-3:] == '.nc':
                    output[model][var][catchment]['unit'] = \
                        str(
                            nc.Dataset(
                                ifiles[model][var]['file']
                                ).variables[var].units).replace(' ', '')
                else:
                    logger.warning(
                        'Attention: No unit specified for model %s, '
                        'variable %s', model, var)
                    output[model][var][catchment]['unit'] = ''
                if not KEEPTIMMEAN:
                    os.remove(timmeanfile)  # delete timmeanfile
    # delete maskfiles and gridfile
    spr.call(['rm'] +
             [maskfiles[catchment] for catchment in sorted(maskfiles.keys())] +
             [gridfile] +
             TMPFILES)
    return output


def writefile(data, output, title):
    """write data to file as to be read by function readfile"""
    import csv
    with open(output, 'w') as csvfile:
        w = csv.writer(csvfile)
        w.writerow([title])
        w.writerow([' '])
        for key in sorted(data.keys()):
            val = data[key]['data']
            unit = data[key]['unit']
            if not data[key].get('rel', False):
                w.writerow(
                    [key.ljust(30) + ':' +
                     str(str(val) + ' ' + unit).rjust(17)])
            else:
                rel = data[key]['rel']
                w.writerow(
                    [key.ljust(30) + ':' +
                     str(str(val) + ' ' + unit).rjust(17) +
                     str(' (' + '%6.2f' % (rel) + ' %)')])
    return


def readfile(ifile):
    """read data from file
    input data structure must be of the form
    <title> (e.g. 'Observations:')

    <CATCHMENTNAME1>                     : <value> <unit> ...
    <CATCHMENTNAME2>                     : <value> <unit> ...
    ...
    (... standing for everything behind the <unit> will be ignored)
    """
    import csv
    data = {}
    with open(ifile, 'r') as csvfile:
        w = csv.reader(csvfile, delimiter=':', quoting=csv.QUOTE_NONE,
                       escapechar='"')
        # skip first two rows
        HEADER = w.next()
        HEADER = w.next()
        # read data
        for row in w:
            data[row[0].replace(' ', '')] = {
                'data': float(row[1].split()[0]),
                'unit': row[1].split()[1]}
            try:
                data[row[0].replace(' ', '')]['rel'] = float(
                    row[1].replace('(', '').split()[2])
            except:
                pass
    return data


def check(catchments, defaultreffile, ifiles):
    """check whether all catchments, variables and files exist"""
    logger.info('Check input...')
    for model in sorted(ifiles.keys()):
        for var in sorted(ifiles[model].keys()):
            # check whether catchments exist in reference and catchment-size
            # file
            REFFILE = ifiles[model][var].get('reffile', defaultreffile)
            if REFFILE[-4:] == '.txt':
                refdata = readfile(REFFILE.replace('<<<VAR>>>', var))
                for catchment in sorted(catchments.keys()):
                    if not refdata.get(catchment, False):
                        raise ValueError(
                            catchment + ' does not exist in reference files!'
                            )
            # check whether input files exist
            modelfiles = deepcopy(ifiles[model][var]['file'])
            if isinstance(modelfiles, (str, unicode)):
                modelfiles = [modelfiles]
            for FILE in modelfiles:
                if not os.path.isfile(FILE.split(' ')[-1]):
                    raise ValueError(
                        FILE.split(' ')[-1] + ' does not exist')
    return


# -----------------------------------------------------------------------------
#
# ----------------------------- plotting routines -----------------------------
#
# -----------------------------------------------------------------------------


def plotfile(data,
             output,
             sep=True,
             dataorrel='all',
             fmt={},
             output_format="pdf"):
    logger.info('Plot files...')
    """plot the model data from analysecatchments
    INPUT:
      - data    {<<<MODEL>>>:
                  {<<<VAR>>>:
                       {<<<CATCHMENT>>>:
                            {'data':<data>,
                             'rel':<relative data in %>
                             'unit':<unit>}}}}
     - output   <string>. name of the pdf. '<<<MODEL>>>' in output will be
                replaced by <<<MODEL>>> from data
     - sep      True/False. If False, put all models together in one axis with
                two y-axes if False. If True, or make two axes into one figure
                if 'axes'. Does not have an effect if dataorrel!='all'
     - dataorrel      <string> (Default: 'all'; Possible values: 'data',
                      'rel', 'all'). Plots absolute and relative
                      variables or only the one specified by <string>.
     - fmt      dictionary defining formatoptions for the plot. Possible
                keywords are given in xyplotformatter
    """
    def set_ylim(data, target):
        if np.min(data) < -200.0:
            target['ylim'] = [-200.0,
                              ['rounded', np.max(data)]]
        elif np.min(data) >= 0:
            target['ylim'] = [0.0,
                              ['rounded', np.max(data)]]
        else:
            target['ylim'] = ['auto', 'auto']
        if np.max(data) > 200.0:
            target['ylim'][1] = 200.0
            if target['ylim'][0] == 'auto':
                target['ylim'][0] = ['rounded', np.min(data)]
        elif np.max(data) <= 0.0:
            target['ylim'][1] = 0.0
            if target['ylim'][0] == 'auto':
                target['ylim'][0] = ['rounded', np.min(data)]
        else:
            pass

    fmt = deepcopy(fmt)
    for model in sorted(data.keys()):
        plotfmt = deepcopy(fmt)
        ofile = output.replace('<<<MODEL>>>', model)
        if dataorrel == 'all':
            plotdata = {var: [
                [data[model][var][catchment]['data']
                 for catchment in sorted(data[model][var].keys())]]
                for var in sorted(data[model].keys())}
            reldata = {var: [
                [data[model][var][catchment]['rel']
                 for catchment in sorted(data[model][var].keys())]]
                for var in sorted(data[model].keys())}
        else:
            plotdata = {var: [
                [data[model][var][catchment][dataorrel]
                 for catchment in sorted(data[model][var].keys())]]
                for var in sorted(data[model].keys())}
            reldata = None
        formatoptions = {var: {
            'ylabel': "%s [%s]" % (
                var, data[model][var][
                    sorted(data[model][var].keys())[0]]['unit']),
            'xlabel': 'Catchment',
            'xlim': [-1, len(plotdata[var][0])],
            # Catchments as xticks
            'xticks': range(len(sorted(data[model][var].keys()))),
            'xticklabels': sorted(data[model][var].keys()),
            'ticksize': 'small',
            'xrotation': 90,
            'tight': True,
            'title': 'Bias of ' + var + ' for ' + model,
            'ax2': {'ylabel': 'relative bias [%]',
                    'ticksize': 'small'
                    }
            } for var in sorted(data[model].keys())}
        # restrict relative axis to maximum 200%
        if reldata is not None:
            for var in formatoptions:
                set_ylim(reldata[var], formatoptions[var]['ax2'])

        # -- delete default values which are set in this routine but not in
        # xyplotformatter
        defaultfmt = [('ylim', 'minmax')]
        keyswithdefault = [defaultpair[0] for defaultpair in defaultfmt]
        for key, val in plotfmt.items():
            if isinstance(val, dict):
                for key2, val2 in val.items():
                    if isinstance(val2, dict):
                        for key3, val3 in val2.items():
                            if (key3, val3) in defaultfmt:
                                del plotfmt[key][key2][key3]
                    if (key2, val2) in defaultfmt or val2 == {}:
                        del plotfmt[key][key2]
            if (key, val) in defaultfmt or val == {}:
                del plotfmt[key]
        # -- add variables to plotfmt if not existing
        # 1) get original formatoptions without variables (in case keywords
        # are in the dictionary)
        fmtorig = {key: val for key, val in plotfmt.items() if key not in
                   sorted(data[model].keys())}
        # 2) create new formatoptions with the original formatoptions
        fmtnew = {var: fmtorig for var in formatoptions if var not in
                  plotfmt.keys()}
        # 3) update old formatoptions
        plotfmt.update(fmtnew)
        # -- configure individual formatoptions and plot
        # if plot shall be on two axes in one figure
        if sep and dataorrel == 'all':
            for var in sorted(formatoptions.keys()):
                for key in sorted(formatoptions[var].keys()):
                    if key != 'ax2' and key != 'title':  # title is set below
                        if isinstance(formatoptions[var][key], (unicode, str)):
                            formatoptions[var]['ax2'].setdefault(
                                key, formatoptions[var][key].replace(
                                    var, 'relative ' + var))
                        else:
                            formatoptions[var]['ax2'].setdefault(
                                key, formatoptions[var][key])
                # set title
                formatoptions[var]['ax2']['title'] = \
                    'Relative bias of ' + var + ' for ' + model
                setdefaultfromdict(plotfmt[var], formatoptions[var])
                setdefaultfromdict(plotfmt[var]['ax2'],
                                   formatoptions[var]['ax2'])
            make_barplot(y=plotdata,
                         y2=reldata,
                         output=ofile,
                         formatoptions=plotfmt,
                         ax=(1, 2),
                         output_format=output_format)
        # if plot shall include only data or relative
        elif dataorrel != 'all':
            if dataorrel == 'rel':
                for var in sorted(formatoptions.keys()):
                    formatoptions[var].update(formatoptions[var]['ax2'])
                    formatoptions[var]['title'] = \
                        'Relative bias of ' + var + ' for ' + model
                    set_ylim(plotdata[var], formatoptions[var])
                    # ensure that no wrong default values are used during the
                    # setdefault from ax2
                    for key, val in plotfmt[var].items():
                        if key in keyswithdefault:
                            del plotfmt[var][key]
                    plotfmt[var].update(plotfmt[var].get('ax2', {}))
                    setdefaultfromdict(plotfmt[var], formatoptions[var])
            else:
                for var in sorted(formatoptions.keys()):
                    setdefaultfromdict(plotfmt[var], formatoptions[var])
            make_barplot(y=plotdata,
                         output=ofile,
                         formatoptions=plotfmt,
                         output_format=output_format)
        # if plot shall be on one axes in one figure but with data and relative
        # data
        elif not sep:
            for var in sorted(formatoptions.keys()):
                formatoptions[var].update({
                    'leftcolor': 'blue',
                    'ylabelcolor': 'blue',
                    'rightcolor': 'red',
                    'ylim': [-np.max([abs(np.min(plotdata[var])),
                                     abs(np.max(plotdata[var]))]),
                             np.max([abs(np.min(plotdata[var])),
                                     abs(np.max(plotdata[var]))])],
                    'alpha': 0.5})
                setdefaultfromdict(plotfmt[var], formatoptions[var])
                # set default values for ax2
                formatoptions[var]['ax2'].update({
                    'ylabelcolor': 'red',
                    'ylim': [
                        # use maximum but restrict to 200%
                        -min([np.max([abs(np.min(plotdata[var])),
                                     abs(np.max(plotdata[var]))]),
                              200.0]),
                        min([np.max([abs(np.min(plotdata[var])),
                                     abs(np.max(plotdata[var]))]),
                             200.0])],
                    'alpha': 0.5})
                setdefaultfromdict(plotfmt[var]['ax2'],
                                   formatoptions[var]['ax2'])
            make_barplot(y=plotdata,
                         y2=reldata,
                         output=ofile,
                         formatoptions=plotfmt,
                         output_format=output_format)
    del formatoptions


def make_barplot(y,
                 y2=None,
                 output='lineplots.pdf',
                 x=None,
                 legendlabels=None,
                 formatoptions=None,
                 ax=None,
                 output_format=None):
    # Create a line-plot, optionally with colored x- or y-errorrange
    plotlist = []

    if type(y) is not dict:
        y = {'tmpvar': y}

    # open output class
    if output_format == "pdf":
        pdf = PdfPages(output)

    # loop through variables
    for var in sorted(y.keys()):
        logger.debug("Plotting barplots of %s", var)
        # read plottind data
        ydata = y[var][:]
        if y2 is not None:
            y2data = y2[var][:]
        xdata = [0] * len(ydata)
        if x is None:
            xdata = [range(len(y[var][0]))] * len(ydata)
        else:
            if type(x) is not dict:
                xtmp = x
            else:
                xtmp = x[var]
            try:
                # if x is 2D-array, then choose copy x to xdata, else copy the
                # 1D-array to 2D-array
                len(xtmp[0])
                xdata = xtmp[:]
            except:
                xdata = [xtmp[:] for i in xrange(len(xdata))]
        ticksize = formatoptions.get('ticksize', 'small')
        if legendlabels is None:
            legendlabels = [None] * len(ydata)
        # make plot
        if ax is None:
            fig = plt.figure(1)
            ax1 = plt.subplot(111)
            if y2 is not None:
                ax2 = ax1.twinx()
        else:
            fig, ax0 = plt.subplots(*ax)
            ax1 = ax0[0]
            if y2 is not None:
                ax2 = ax0[1]

        plt.sca(ax1)
        for i in xrange(len(ydata)):
            plot1 = plt.bar(xdata[i], ydata[i], label=legendlabels[i],
                            align='center')
        plt.axhline(c='k', lw=2)
        if y2 is not None:
            plt.sca(ax2)
            plot2 = ax2.bar(xdata[i], y2data[i], label=legendlabels[i],
                            color='red', align='center')
            plt.axhline(c='k', lw=2)
        if legendlabels != [None] * len(ydata):
            ax1.legend(loc='best', fontsize=ticksize)

        # edit format
        if not not formatoptions:
            if var in formatoptions:
                if 'ax2' in formatoptions[var] and y2 is not None:
                    fmt2 = {
                        key: val for key, val in
                        formatoptions[var]['ax2'].items() if key != 'alpha'}
                    for bar in plot2:
                        bar.set_alpha(
                            formatoptions[var]['ax2'].get('alpha', None))
                else:
                    fmt2 = None
                fmt = {
                    key: val for key, val in formatoptions[var].items()
                    if key not in ['ax2', 'alpha']}
                for bar in plot1:
                    bar.set_alpha(formatoptions[var].get('alpha', None))
            else:
                if 'ax2' in formatoptions and y2 is not None:
                    fmt2 = {
                        key: val for key, val in
                        formatoptions['ax2'].items() if key != 'alpha'}
                    plot2[0].set_alpha(formatoptions['ax2'].get('alpha', None))
                else:
                    fmt2 = None
                fmt = {
                    key: val for key, val in formatoptions[var].items()
                    if key not in ['ax2', 'alpha']}
                plot1[0].set_alpha(formatoptions.get('alpha', None))
            if fmt2 is not None:
                xyplotformatter(ax=ax2, formatoptions=fmt2)
                plt.sca(ax1)
            xyplotformatter(ax=ax1, formatoptions=fmt)
        # make output
        if output_format == "pdf":
            pdf.savefig(fig, dpi=300)
        else:
            output_with_var = re.sub("(\.[a-zA-Z]*)$", "-" + var + r"\1", output)
            fig.savefig(output_with_var)
        plt.close(1)
    if output_format == "pdf":
        pdf.close()


def make_coefficient_plots(data,
                           diff2ref,
                           refdata,
                           output_ro_et,
                           output_ro_pr,
                           output_format):
    """Makes the plot of the biases of runoff coeffient (runoff over
    precipitation) against evapotranspiration coefficient (evapotrans. over
    precipiation) for all models in data

    Input:
      - data: input dictionary as got resulted from function computemeans
      - refdata: reference data dictionary
      - output: output file name
    """
    def calc_coeffienct_biases(data, refcoeffs=None):
        """Function to calculate the coefficients and their biases"""
        if refcoeffs is None:
            logger.debug(
                "No reference data specified --> set reference data to 0")
            refcoeffs = {
                catchment: [0, 0] for catchment in data['runoff'].keys()}
        coeff_biases = {
            catchment: [0., 0.] for catchment in data['runoff'].keys()}
        for catchment, val in coeff_biases.items():
            ro = data['runoff'][catchment]['data']
            pr = data['precip'][catchment]['data']
            et = data['ET'][catchment]['data']
            val[0] = et / pr * 100. - refcoeffs[catchment][0]
            val[1] = ro / pr * 100. - refcoeffs[catchment][1]
            logger.debug("Catchment %s:", catchment)
            logger.debug("    Runoff: %s, Precip: %s, reference: %s", ro, pr,
                         refcoeffs[catchment][0])
            logger.debug("        --> Bias: %s", val[0])
            logger.debug("Catchment %s:", catchment)
            logger.debug(
                "    Evapotranspiration: %s, Precip: %s, reference: %s", et,
                pr, refcoeffs[catchment][1])
            logger.debug("        --> Bias: %s", val[1])

        return coeff_biases

    def coefficient_plot(data, xlabel='', ylabel='', ax=None):
        """Calculates coefficients and plots them

        Input:
        - data: dictionary for one model with bias difference
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        marker = cycle(('.', '+', 'o', '*', 'x'))
        # make axis symmetric
        ax.xaxis.set_major_locator(MaxNLocator(symmetric=True))
        ax.yaxis.set_major_locator(MaxNLocator(symmetric=True))
        for catchment, val in sorted(data.items()):
            ax.plot([val[0]], [val[1]], next(marker), label=catchment, mew=2)

        fig.subplots_adjust(bottom=0.3)
        legend = ax.legend(loc='upper left',
                           bbox_to_anchor=(-0.1, -0.1),
                           ncol=3,
                           numpoints=1)
        plt.title(model)
        ax.spines['left'].set_position('center')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('center')
        ax.spines['top'].set_color('none')
        # remove 0 from y and x axis
        for key in ['x', 'y']:
            ticks = list(getattr(ax, 'get_%sticks' % key)())
            try:
                ticks.remove(0)
            except ValueError:
                pass
            getattr(ax, 'set_%sticks' % key)(ticks)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        # xlabel below
        plt.xlabel(xlabel)
        ax.xaxis.set_label_coords(0.5, -0.025)
        # ylabel left
        plt.ylabel(ylabel)
        ax.yaxis.set_label_coords(-0.025, 0.5)

        return fig

    logger.debug("Making scatter plot")
    if output_format == "pdf":
        pdf_ro_et = PdfPages(output_ro_et)
        pdf_ro_pr = PdfPages(output_ro_pr)

    for model, model_dict in data.items():
        missing_keys = [key for key in ['ET', 'runoff', 'precip']
                        if key not in model_dict]
        if missing_keys != []:
            logger.warning(
                "Can not produce scatter plot because %s is missing",
                " and ".join(missing_keys))
            continue
        refcoeffs = calc_coeffienct_biases(refdata[model])
        coefficients = calc_coeffienct_biases(model_dict, refcoeffs)

        ofig = coefficient_plot(coefficients,
                                xlabel="Bias of ET coefficient [%]",
                                ylabel="Bias of runoff coefficient [%]")

        if output_format == "pdf":
            pdf_ro_et.savefig(ofig)
        else:
            ofig.savefig(output_ro_et)

        for catchment, val in coefficients.items():
            val[0] = diff2ref[model]['precip'][catchment]['rel']

        ofig = coefficient_plot(coefficients,
                                xlabel="Relative bias of precipitation [%]",
                                ylabel="Bias of runoff coefficient [%]")

        if output_format == "pdf":
            pdf_ro_pr.savefig(ofig)
        else:
            ofig.savefig(output_ro_pr)

    if output_format == "pdf":
        pdf_ro_et.close()
        pdf_ro_pr.close()


def xyplotformatter(ax, formatoptions):
    """formats title and x- und y-axes of a plot
    Input:
        - ax: matplotlib.axes.AxesSubplot. Axis on which the format
            configuration shall be done
        - formatoptions: dictionary. Contains the format options. Possible keys
            are
            -- 'title': string (Default: None). Gives the title of the plot
            -- 'ylabel': string (Default: None). Defines the y-axis label
            -- 'ylabelcolor': string or color for y-axis label (Default: None)
            -- 'xlabel': string (Default: None). Defines the x-axis label
            -- 'xlabelcolor': string or color for x-axis label (Default: None)
            -- 'yformat': string ('sci' or 'int') (Default: None). Gives the
                    yaxis ticks either in scientific ('sci') or integer ('int')
                    format
            -- 'xformat': string ('sci' or 'int') (Default: None). Gives the
                    xaxis ticks either in scientific ('sci') or integer ('int')
                    format
            -- 'ydeci': integer (Default 0). If 'yformat' is 'sci', ydeci
                    represents the number of decimals
            -- 'xdeci': integer (Default 0). If 'yformat' is 'sci', xdeci
                    represents the number of decimals
            -- 'xrotation': float (Default 0). Degrees between 0 and 360 for
                    which the xticklabels shall be rotated
            -- 'ylim': [ymin,ymax] (Default: None). Specifies the limits of the
                    y-axis. ymin and ymax can be either floats, 'auto', or of
                    the form ['rounded', val]. If 'auto', the value choosed by
                    pyplot is used. if ['rounded', val], the next 0.5-value
                    with respect to its exponent with base 10 is used (i.e.
                    5.4 is rounded to 5.5, 54 is rounded to 55, etc.)
            -- 'xlim': [xmin, xmax] (Default: None). Specifies the limits of
                    the x-axis (same options as for ylim)
            -- 'yscale': string ('linear', 'log' or 'symlog') (Default: None).
                    Sets the scale of the y-axis (See also 'yformat' for
                    scientific axis description)
            -- 'xscale': string ('linear', 'log' or 'symlog') (Default: None).
                    Sets the scale of the x-axis (See also 'xformat' for
                    scientific axis description)
            -- 'yticks': integer or list (Default: None). Defines the y-ticks.
                    If integer i, every i-th tick of the automatically
                    calculated ticks will be used.
            -- 'xticks': integer or list (Default: None). Defines the x-ticks.
                    If integer i, every i-th tick of the automatically
                    calculated ticks will be used.
            -- 'yticklabels': 1D-array (Default: None). Defines the
                    y-axis-ticklabels. If None, the automatically calculated
                    x-ticks will be used.
            -- 'ytickcolor': string or color for y-axis ticklabels (Default:
                    None)
            -- 'xticklabels': 1D-array (Default: None). Defines the
                    x-axis-ticklabels. If None, the automatically calculated
                    x-ticks will be used.
            -- 'xtickcolor': string or color for x-axis ticklabels (Default:
                    None)
            -- 'leftcolor':  string or color for left axis (Default: None)
            -- 'rightcolor':  string or color for right axis (Default: None)
            -- 'topcolor':  string or color for top axis (Default: None)
            -- 'bottomcolor':  string or color for bottom axis (Default: None)
            -- 'fontweight': A numeric value in the range 0-1000 or string.
                    Defines the weight of the labels. Possible strings are one
                    of ‘ultralight’, ‘light’,     ‘normal’, ‘regular’, ‘book’,
                    ‘medium’, ‘roman’, ‘semibold’, ‘demibold’, ‘demi’,
                    ‘bold’, ‘heavy’, ‘extra bold’, ‘black’.
            -- ticksize: string or float. Defines the size of the ticks.
                Strings might be ‘xx-small’, ‘x-small’, ‘small’, ‘medium’,
                ‘large’, ‘x-large’, ‘xx-large’. Floats define the absolute font
                size, e.g., 12
            -- labelsize: string or float. Defines the size of x- and y-axis
                labels (see ticksize for possible values)
            -- titlesize: string or float. Defines the size of the title (see
                ticksize for possible values)
            -- tight: Boolean (Default: False). Make tight_layout if True.
    """
    def get(D, k, d=None):
        """enhancement of dicts get method to set allowed keys
        Input: D: dict instance, k: key, d: default value
        """
        defaultkeys.append(k)
        return D.get(k, d)
    defaultkeys = []
    # read labeling stuff
    title = get(formatoptions, 'title')
    units = get(formatoptions, 'units', '')
    ylabel = get(formatoptions, 'ylabel')
    ylabelcolor = get(formatoptions, 'ylabelcolor')
    leftcolor = get(formatoptions, 'leftcolor')
    rightcolor = get(formatoptions, 'rightcolor')
    topcolor = get(formatoptions, 'topcolor')
    bottomcolor = get(formatoptions, 'bottomcolor')
    xlabel = get(formatoptions, 'xlabel')
    xlabelcolor = get(formatoptions, 'xlabelcolor')
    xformat = get(formatoptions, 'xformat')  # numbering format x-axis
    xdeci = get(formatoptions, 'xdeci', 0)
    yformat = get(formatoptions, 'yformat')  # numbering format y-axis
    ydeci = get(formatoptions, 'ydeci', 0)
    ylim = get(formatoptions, 'ylim')
    xlim = get(formatoptions, 'xlim')
    yscale = get(formatoptions, 'yscale')
    xscale = get(formatoptions, 'xscale')
    yticks = get(formatoptions, 'yticks')
    yticklabels = get(formatoptions, 'yticklabels')
    ytickcolor = get(formatoptions, 'ytickcolor')
    xticks = get(formatoptions, 'xticks')
    xticklabels = get(formatoptions, 'xticklabels')
    xtickcolor = get(formatoptions, 'xtickcolor')
    xrotation = get(formatoptions, 'xrotation', 0)
    fontweight = get(formatoptions, 'fontweight')
    ticksize = get(formatoptions, 'ticksize', 'small')
    labelsize = get(formatoptions, 'labelsize', 'medium')
    titlesize = get(formatoptions, 'titlesize', 'large')
    tight = get(formatoptions, 'tight', False)

    # check whether there are any wrong keys
    for key in formatoptions:
        if key not in defaultkeys:
            similarkeys = get_close_matches(key, defaultkeys)
            if similarkeys == []:
                raise KeyError(
                    'Unknown formatoption keyword ' + key)
            else:
                raise KeyError(
                    'Unknown formatoption keyword ' + key +
                    '! Possible similiar frasings are ' +
                    ', '.join(key for key in similarkeys) + '.')

    if fontweight is None:
        fontformatter = {}
        fontformatter = {label: None for label in
                         ['xticks', 'yticks', 'title', 'ylabel', 'xlabel']}
    elif type(fontweight) is dict:
        fontformatter = {}
        for label in ['xticks', 'yticks', 'title', 'ylabel', 'xlabel']:
            fontformatter[label] = fontweight.get(label, None)
    else:
        fontformatter = {}
        fontformatter = {label: fontweight for label in
                         ['xticks', 'yticks', 'title', 'ylabel', 'xlabel']}

    if yscale:
        ax.set_yscale(yscale)
    if xscale:
        ax.set_xscale(xscale)
    modes = ['min', 'max']
    plotylim = ax.get_ylim()
    if ylim is not None:
        ylim = list(ylim)
        for i, val in enumerate(ylim):
            if val == 'auto':
                ylim[i] = plotylim[i]
            try:
                if val[0] == 'rounded':
                    ylim[i] = round_to_05(val[1], mode=modes[i])
            except (TypeError, IndexError):
                pass
    plt.ylim(ylim)

    plotxlim = ax.get_xlim()
    if xlim is not None:
        xlim = list(xlim)
        for i, val in enumerate(xlim):
            if val == 'auto':
                xlim[i] = plotxlim[i]
            try:
                if val[0] == 'rounded':
                    xlim[i] = round_to_05(val[1], mode=modes[i])
            except TypeError:
                pass
    plt.xlim(xlim)

    # modify yticks
    if not yticks:
        yticks = ax.get_yticks()[:]
    elif type(yticks) is list:
        ax.set_yticks(yticks)
    else:
        yticks = ax.get_yticks()[::yticks]
        ax.set_yticks(yticks)
    if not yticklabels:
        yticklabels = ax.get_yticks()[:]
    ax.set_yticklabels(yticklabels, fontweight=fontformatter['yticks'],
                       fontsize=ticksize)
    if ytickcolor:
        ax.tick_params(axis='y', colors=ylabelcolor)

    # modify xticks
    if not xticks:
        xticks = ax.get_xticks()[:]
    elif type(xticks) is list:
        ax.set_xticks(xticks)
    else:
        xticks = ax.get_xticks()[::xticks]
        ax.set_xticks(xticks)
    if not xticklabels:
        xticklabels = ax.get_xticks()[:]
    ax.set_xticklabels(xticklabels, fontweight=fontformatter['xticks'],
                       fontsize=ticksize, rotation=xrotation)
    if xtickcolor:
        ax.tick_params(axis='x', colors=ylabelcolor)

    # change axis format
    if xformat == 'int':
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.0f'))
    elif xformat == 'sci':
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.xaxis.set_major_formatter(
            FormatStrFormatter('%1.' + str(xdeci) + 'e'))
    if yformat == 'int':
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.0f'))
    elif yformat == 'sci':
        ax.get_yaxis().set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(
            FormatStrFormatter('%1.' + str(ydeci) + 'e'))

    # label stuff
    if title:
        plt.title(title, fontweight=fontformatter['title'], fontsize=titlesize)
    if ylabel:
        plt.ylabel(ylabel, fontweight=fontformatter['ylabel'],
                   fontsize=labelsize)
    if leftcolor:
        ax.spines['left'].set_color(leftcolor)
    if rightcolor:
        ax.spines['right'].set_color(rightcolor)
    if topcolor:
        ax.spines['top'].set_color(topcolor)
    if bottomcolor:
        ax.spines['bottom'].set_color(bottomcolor)
    if ylabelcolor:
        ax.yaxis.label.set_color(ylabelcolor)
    if xlabel:
        plt.xlabel(xlabel, fontweight=fontformatter['xlabel'],
                   fontsize=labelsize)
    if xlabelcolor:
        ax.xaxis.label.set_color(xlabelcolor)
    if tight:
        plt.tight_layout()


#def round_to_05(n, exp=None, func=round):
  #"""Applies the round function specified in func to round n to the
  #next 0.5-value with respect to its exponent with base 10 (i.e.
  #1.3e-4 will be rounded to 1.5e-4) if exp is None or with respect
  #to the given exponent in exp.
  #Input:
    #- n: Float, number to round
    #- exp: Integer. Exponent for rounding
    #- func: Rounding function
  #Output:
    #- Rounded n
  #"""
  #from math import log10, floor
  #if exp is None:
    #exp=floor(log10(abs(n))) # exponent for base 10
  #ntmp=np.abs(n)/10.**exp # mantissa for base 10
  #print exp, ntmp
  #if np.abs(func(ntmp) - ntmp) >= 0.5:
    #return np.sign(n)*(func(ntmp) - 0.5)*10.**exp
  #else:
    #return np.sign(n)*func(ntmp)*10.**exp


def round_to_05(n, exp=None, mode='nearest'):
    """Applies the round function specified in func to round n to the
    next 0.5-value with respect to its exponent with base 10 (i.e.
    1.3e-4 will be rounded to 1.5e-4) if exp is None or with respect
    to the given exponent in exp.
    Input:
        - n: Float, number to round
        - exp: Integer. Exponent for rounding
        - mode: string ('nearest', 'max' or 'min'). If 'nearest', it will be
            rounded to the nearest 0.5-value, if 'max', it will be rounded to
            the next bigger 0.5-value, if 'min', it will be rounded to the next
            smaller 0.5-value
    Output:
        - Rounded n
    """
    from math import log10, floor
    if exp is None:
        exp = floor(log10(abs(n)))  # exponent for base 10
    ntmp = n / 10. ** exp  # mantissa for base 10
    remainder = np.abs(np.floor(ntmp) - ntmp)
    if mode == 'nearest':
        if remainder < 0.25:
            res = np.floor(ntmp) * 10. ** exp
        elif remainder >= 0.75:
            res = np.ceil(ntmp) * 10. ** exp
        else:
            res = (np.floor(ntmp) + 0.5) * 10. ** exp
    elif mode == 'max':
        if remainder > 0.5:
            res = np.ceil(ntmp) * 10. ** exp
        else:
            res = (np.floor(ntmp) + 0.5) * 10. ** exp
    elif mode == 'min':
        if remainder < 0.5:
            res = np.floor(ntmp) * 10. ** exp
        else:
            res = (np.floor(ntmp) + 0.5) * 10. ** exp
    return res


# -----------------------------------------------------------------------------
#
# --------------------reference calculation routines --------------------------
#
# -----------------------------------------------------------------------------

def make_reference_calculations(defaultreffile,
                                PREF,
                                MRROREFCALCULATION,
                                ETREFCALCULATION,
                                PRECIPREFCALCULATIN):
    """Function to make reference calculations for et, mrro and precip"""
    if MRROREFCALCULATION:
        runoff_ref_calculation(
            PREF + 'ref_runoff_catchments_old.txt',
            PREF + 'catchment-areas.txt',
            PREF + 'ref_runoff_catchments.txt')
    elif MRROREFCALCULATION is None:
        runoffrefdata = defaults.runoffrefdata
        if os.path.exists(defaultreffile.replace('<<<VAR>>>', 'runoff')):
            logger.warning(
                'Runoff file already exists! Copy old file to %s~',
                defaultreffile.replace('<<<VAR>>>', 'runoff'))
            spr.call(['cp', defaultreffile.replace('<<<VAR>>>', 'runoff'),
                      defaultreffile.replace('<<<VAR>>>', 'runoff') + '~'])
        writefile(runoffrefdata, defaultreffile.replace('<<<VAR>>>', 'runoff'),
                  'Observations:')
    if ETREFCALCULATION:
        if precfile is None:
            precfile = raw_input(
                ("Please type the filename of the"
                    "file with precipitation data stored in variable"
                    "'Rainf_GPCC' with units in mm/a for the calculation of ET"
                    " reference.\n"))
        ET_ref_calculation(precfile, PREF + 'ref_runoff_catchments.txt',
                           pcatchment, catchments,
                           PREF + 'ref_ET_catchments.txt')
    elif ETREFCALCULATION is None:
        ETrefdata = defaults.ETrefdata
        if os.path.exists(defaultreffile.replace('<<<VAR>>>', 'ET')):
            logger.warning(
                'Evapotranspiration file already exists! Copy old file to %s~',
                defaultreffile.replace('<<<VAR>>>', 'ET'))
            spr.call(['cp', defaultreffile.replace('<<<VAR>>>', 'ET'),
                      defaultreffile.replace('<<<VAR>>>', 'ET') + '~'])
        writefile(ETrefdata, defaultreffile.replace('<<<VAR>>>', 'ET'),
                  'Observations:')
    if PRECIPREFCALCULATIN:
        if precfile is None:
            precfile = raw_input(
                ("Please type the filename of the"
                    "file with precipitation data stored in variable"
                    "'Rainf_GPCC' with units in mm/a for the calculation of ET"
                    " reference.\n"))
        precip_ref_calculation(precfile, pcatchment, catchments,
                               PREF + 'ref_precip_catchments.txt')
    elif PRECIPREFCALCULATIN is None:
        preciprefdata = defaults.preciprefdata
        if os.path.exists(defaultreffile.replace('<<<VAR>>>', 'precip')):
            logger.warning(
                'Evapotranspiration file already exists! Copy old file to %s~',
                defaultreffile.replace('<<<VAR>>>', 'precip'))
            spr.call(['cp', defaultreffile.replace('<<<VAR>>>', 'precip'),
                      defaultreffile.replace('<<<VAR>>>', 'precip') + '~'])
        writefile(preciprefdata, defaultreffile.replace('<<<VAR>>>', 'precip'),
                  'Observations:')


def runoff_ref_calculation(oldfile, catchmentareas, newfile):
    """Runoff reference calculation routines

    transforms reference data from m^3/s to mm/a (and translate from German to
    English)
    Input:
        - oldfile: <string>. Path to to file containing old data in m^3/s,
        - catchmentareas: <string>. Path to file containing data with catchment
            areas
        - newfile: <string. Name of the new file
    """
    logger.info('Calculate runoff reference')
    olddata = readfile(oldfile)
    # change keys to english
    olddata['Nile'] = olddata.pop('Nil')
    olddata['Danube'] = olddata.pop('Donau')
    olddata['Amazon'] = olddata.pop('Amazonas')
    olddata['Yangtze-Kiang'] = olddata.pop('Yangtze')
    olddata['Baltic-Sea-cat'] = olddata.pop('Baltic')
    olddata['Ganges-Brahmaputra'] = olddata.pop('Ganges')
    olddata['6-largest-Arctic-Rivers'] = olddata.pop('ACriv')
    areas = readfile(catchmentareas)
    seconds = 60 * 60 * 24 * 365.25  # seconds per year
    newdata = {}
    for CATCHMENT in sorted(areas.keys()):
        if CATCHMENT in sorted(olddata.keys()):
            newdata[CATCHMENT] = {
                'data': "%.4f" % (
                    olddata[CATCHMENT]['data'] /
                    areas[CATCHMENT]['data'] * 1000 * seconds),
                'unit': 'mm/a'}

    writefile(newdata, newfile, 'Observations:')


def ET_ref_calculation(precfile, runofffile, pcatchment, catchments, newfile,
                       vname='Rainf_GPCC'):
    """Evapotranspiration reference routines

    transform reference data from m^3/s to mm/a
    INPUT:
        - precfile: <string>. Path to file containing grid cell based data on
            observed precipitation (handable for cdos, e.g. netCDF) (units in
            mm/s=kg m-2 s-1).
        - runofffile: <string>. Path to file containing catchment based data on
            observed runoff (units in mm/a).
        - pcatchment: <string>. Path to file containing catchment definitions.
        - catchments: dictionary. Catchment names and ids as used in
            pcatchment.
        - newfile: <string>. Name of the new file
        - vname: <string>. Name of precipitation variable in precfile
    """
    logger.info('Calculate ET reference')
    ifiles = {'ref': {'precip': {
        'file': precfile, 'vname': vname, 'cdo': '-mulc,31471200 ',
        'unit': 'mm/a'}}}
    POUT = os.path.dirname(os.path.abspath(newfile)) + os.path.sep
    refdata = computemean(ifiles, pcatchment, catchments, POUT,
                          False)['ref']['precip']
    runoffdata = readfile(runofffile)
    for catchment in sorted(refdata.keys()):
        refdata[catchment]['data'] = (refdata[catchment]['data'] -
                                      runoffdata[catchment]['data'])
    writefile(refdata, newfile, 'Observations:')


def precip_ref_calculation(precfile, pcatchment, catchments, newfile,
                           vname="Rainf_GPCC"):
    logger.info('Calculate precipitation reference')
    ifiles = {'ref': {'precip': {
        'file': precfile, 'vname': vname, 'cdo': '-mulc,31471200 ',
        'unit': 'mm/a'}}}
    POUT = os.path.dirname(os.path.abspath(newfile)) + os.path.sep
    refdata = computemean(ifiles, pcatchment, catchments, POUT,
                          False)['ref']['precip']
    writefile(refdata, newfile, 'Observations:')


def setdefaultfromdict(D1, D2):
    """function to setdefault from a dictionary instead of key, val pair"""
    for key, val in D2.items():
        D1.setdefault(key, val)


class defaults(object):
    """Class containing default dictionaries

    The properties are used in analysecatchments.
    Properties are
        catchments
        defaultreffile
        runoffrefdata
        ETrefdata
    """
    catchments = {
        # Catchments with name as used in REFFILE as key and the
        # catchment number as used in pcatchment as value
        "Amazon": 94,
        "Parana": 98,
        "Mackenzie": 76,
        "Mississippi": 86,
        "Danube": 14,
        "Congo": 68,
        "Niger": 65,
        "Nile": 60,
        "Lena": 40,
        "Yangtze-Kiang": 52,
        "Ganges-Brahmaputra": 54,
        "Murray": 100
        }

    defaultreffile = "ref_<<<VAR>>>_catchments.txt"

    runoffrefdata = {
        'Amazon': {'data': 1195.4477, 'unit': 'mm/a'},
        'Congo': {'data': 365.6980, 'unit': 'mm/a'},
        'Danube': {'data': 250.9211, 'unit': 'mm/a'},
        'Ganges-Brahmaputra': {'data': 672.5738, 'unit': 'mm/a'},
        'Lena': {'data': 197.3081, 'unit': 'mm/a'},
        'Mackenzie': {'data': 173.9881, 'unit': 'mm/a'},
        'Mississippi': {'data': 182.2420, 'unit': 'mm/a'},
        'Murray': {'data': 8.2041, 'unit': 'mm/a'},
        'Niger': {'data': 31.5160, 'unit': 'mm/a'},
        'Nile': {'data': 48.7528, 'unit': 'mm/a'},
        'Parana': {'data': 203.0060, 'unit': 'mm/a'},
        'Yangtze-Kiang': {'data': 531.6936, 'unit': 'mm/a'}
        }

    preciprefdata = {
        'Amazon': {'data': 2253.61, 'unit': 'mm/a'},
        'Congo': {'data': 1539.98, 'unit': 'mm/a'},
        'Danube': {'data': 809.11, 'unit': 'mm/a'},
        'Ganges-Brahmaputra': {'data': 1387.95, 'unit': 'mm/a'},
        'Lena': {'data': 399.146, 'unit': 'mm/a'},
        'Mackenzie': {'data': 445.342, 'unit': 'mm/a'},
        'Mississippi': {'data': 890.034, 'unit': 'mm/a'},
        'Murray': {'data': 530.441, 'unit': 'mm/a'},
        'Niger': {'data': 436.907, 'unit': 'mm/a'},
        'Nile': {'data': 673.565, 'unit': 'mm/a'},
        'Parana': {'data': 1311.22, 'unit': 'mm/a'},
        'Yangtze-Kiang': {'data': 1032.84, 'unit': 'mm/a'}
        }

    ETrefdata = {
        'Amazon': {'data': 1014.4023, 'unit': 'mm/a'},
        'Congo': {'data': 1203.182, 'unit': 'mm/a'},
        'Danube': {'data': 554.5999, 'unit': 'mm/a'},
        'Ganges-Brahmaputra': {'data': 722.5962, 'unit': 'mm/a'},
        'Lena': {'data': 187.4469, 'unit': 'mm/a'},
        'Mackenzie': {'data': 269.2429, 'unit': 'mm/a'},
        'Mississippi': {'data': 712.192, 'unit': 'mm/a'},
        'Murray': {'data': 465.1909, 'unit': 'mm/a'},
        'Niger': {'data': 402.23, 'unit': 'mm/a'},
        'Nile': {'data': 602.1752, 'unit': 'mm/a'},
        'Parana': {'data': 1085.554, 'unit': 'mm/a'},
        'Yangtze-Kiang': {'data': 538.0664, 'unit': 'mm/a'}
        }
