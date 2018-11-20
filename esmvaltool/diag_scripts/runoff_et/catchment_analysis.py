#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Catchment specific water flux plots

###############################################################
runoff_et/catchment_analysis.py
Authors ESMValToolV1 Version
    Philipp Sommer (philipp.sommer@mpimet.mpg.de)
    Stefan Hagemann (stefan.hagemann@hzg.de)
    Alexander Loew
Port to ESMValTool Version 2
    Tobias Stacke (tobias.stacke@mpimet.mpg.de)
###############################################################

Description
-----------
    Plots temporal and spatial averages of precipitation, runoff and
    evaporation for specific land surface catchments. Additionally,
    relations of runoff coefficient to relative precipitation bias
    and runoff coefficient to evaporation coefficient are computed.

    Default reference data are included in this routine (default class)
    but can be replaced with other datasets. In case a custom catchment
    mask is used, the default class (catchment names, IDs, reference data)
    has to be adapted.

###############################################################

"""
import iris
import numpy as np
import calendar
import os
import logging
import pdb
import esmvaltool.diag_scripts.shared as diag
from itertools import cycle
from esmvaltool.preprocessor._regrid import regrid
from esmvaltool.preprocessor._area_pp import area_average
import matplotlib
matplotlib.use('Agg')

logger = logging.getLogger(os.path.basename(__file__))


class defaults(object):
    """Class containing default dictionaries for predefined catchments

    The properties are used in the routine analysecatchments. Catchments and
    reference values are specific for the default catchment mask. All reference
    values are given in mm a-1. Precip data is based on WFDEI, runoff is based
    on GRDC, ET is derived as the difference of both. The values are updated
    and differ slightly from the ESMValTool 1 version.
    Properties are
        catchments
        mrro
        pr
        evspsbl
    """
    catchments = {
        # Catchments with name as used in make_catchment_plots and
        # associated ID used in the catchment mask netCDF file
        "Amazon": 94,
        "Parana": 98,
        "Mackenzie": 76,
        "Mississippi": 86,
        "Danube": 14,
        "Congo": 68,
        "Niger_Malanville": 65,
        "Nile": 60,
        "Lena": 40,
        "Yangtze-Kiang": 52,
        "Ganges-Brahmaputra": 54,
        "Murray": 100,
    }

    mrro = {
        'Amazon': 1194.63,
        'Congo': 365.45,
        'Danube': 250.75,
        'Ganges-Brahmaputra': 672.11,
        'Lena': 199.61,
        'Mackenzie': 173.87,
        'Mississippi': 182.12,
        'Murray': 8.20,
        'Niger_Malanville': 31.49,
        'Nile': 48.72,
        'Parana': 202.87,
        'Yangtze-Kiang': 531.33,
    }

    pr = {
        'Amazon': 2210.25,
        'Congo': 1571.41,
        'Danube': 808.04,
        'Ganges-Brahmaputra': 1405.84,
        'Lena': 387.01,
        'Mackenzie': 450.16,
        'Mississippi': 897.18,
        'Murray': 474.62,
        'Niger_Malanville': 437.90,
        'Nile': 655.62,
        'Parana': 1314.66,
        'Yangtze-Kiang': 1074.79,
    }

    evspsbl = {
        'Amazon': 1015.62,
        'Congo': 1205.96,
        'Danube': 557.29,
        'Ganges-Brahmaputra': 733.73,
        'Lena': 187.40,
        'Mackenzie': 276.29,
        'Mississippi': 715.06,
        'Murray': 466.42,
        'Niger_Malanville': 406.41,
        'Nile': 606.90,
        'Parana': 1111.80,
        'Yangtze-Kiang': 543.46,
    }


def format_coef_plot(ax):
    """ Moves axis from border to center, adapts ticks and labels accordingly
    Parameters
    ----------
    ax : object
        plot axis object
    """
    # Add infos to axis
    ax.xaxis.set_label_coords(0.5, -0.025)
    ax.yaxis.set_label_coords(-0.025, 0.5)
    # Adapt axis range to center zero
    xmax = np.ceil(
        (np.absolute(np.array(ax.get_xlim())).max() + 5) / 10.0) * 10.0 - 5.0
    ax.set_xlim(xmax * -1, xmax)
    ymax = np.ceil(
        (np.absolute(np.array(ax.get_ylim())).max() + 5) / 10.0) * 10.0 - 5.0
    ax.set_ylim(ymax * -1, ymax)
    # remove 0 from y and x axis
    for key in ['x', 'y']:
        ticks = list(getattr(ax, 'get_%sticks' % key)())
        try:
            ticks.remove(0)
        except ValueError:
            pass
        getattr(ax, 'set_%sticks' % key)(ticks)

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    return ax

def data2file(cfg, filename, title, filedata):
    """ Output data dictionary into ascii file
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    filename : str
        String containing the file name
    title : str
        String containing the file header
    filedata : dict
        Dictionary of catchment averages per river
    """

    # Write experiment data
    filepath = os.path.join(cfg[diag.names.WORK_DIR], filename)
    with open(filepath, 'w') as f:
        f.write(title + '\n\n')
        for river, value in sorted(filedata.items()):
            f.write('{:25} : {:8.2f}\n'.format(river, value))


def write_plotdata(cfg, plotdata, catch_info, reference):
    """ Output catchment averaged values for all datasets.
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    plotdata : dict
        Dictionary containing the catchment averages
    catch_info : object
        Object containing catchment names, IDs, and reference data
    reference : str
        String containing name of the reference dataset
    """

    ref_vars = []
    metric   = "catchment averages"
    unit     = "[mm a-1]"

    for var in plotdata.keys():
        for identifier in plotdata[var].keys():
            # Write experiment data
            filename = '_'.join([var, identifier]) + '.txt'
            title = " ".join(identifier.split(' ') + [var, metric, unit])
            filedata = plotdata[var][identifier]
            data2file(cfg, filename, title, filedata)
            # Write reference data
            if var not in ref_vars:
                filename = '_'.join([var, 'reference']) + '.txt'
                title = " ".join([reference, metric, unit])
                filedata = getattr(catch_info, var)
                data2file(cfg, filename, title, filedata)
                ref_vars.append(var)


def get_expdata(datadict, rivers):
    """ Get list with catchment averages sorted according to
        river list from reference
    Parameters
    ----------
    datadict : dict
        the catchment averages data dictionary
    rivers : list
        list of river catchment names
    """

    expdata = []
    for riv in rivers:
        expdata.append(datadict[riv])
    return np.array([expdata])


def prep_barplot(title, rivers, var):
    """ Prepare barplot
    Parameters
    ----------
    title : str
        multipanel plot title
    rivers : list
        list of river catchment names
    var : str
        short name of the actual variable
    """
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=False)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.35)
    plottitle = ['\nBias for ', '\nRelative bias for ']
    ylabel = [var.upper() + ' [mm a-1]','Relative bias [%]']

    for ia, ax in enumerate(axs.keys()):
        ax.set_title(plottitle[ia] + var.upper())
        ax.set_ylabel(ylabel[ia])
        ax.set_xlabel('Catchment')
        ax.set_xticks(range(len(rivers)))
        ax.set_xticklabels((rivers), fontsize='small')
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        ax.axhline(c='black', lw=2)

    return fig, axs


def prep_scatplot(title, coeftype):
    """ Prepare scatterplot for different coefficients
    Parameters
    ----------
    title : str
        multipanel plot title
    coeftype : str
        string indicting plot type [prbias,etcoef]
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False)
    ax.set_title(title)
    ax.set_ylabel('Bias of runoff coefficient [%]')
    if coeftype == 'prbias':
        ax.set_xlabel('Relative bias of precipitation [%]')
    elif coeftype == 'etcoef':
        ax.set_xlabel('Bias of ET coefficient [%]')
    else:
        raise ValueError('Unexpected coefficient combination in prep_scatplot')

    return fig, ax


def add_legend(fig, rivers, markerlist):
    """ Prepare scatterplot for different coefficients
    Parameters
    ----------
    fig : obj
        plot figure object
    rivers : list
        list of river catchment names
    markerlist : list
        list of marker strings for scatterplot legend
    """

    # Define legend
    fig.subplots_adjust(bottom=0.30)
    marker = cycle(markerlist)
    caxe = fig.add_axes([0.05, 0.01, 0.9, 0.20])
    for i, label in enumerate(rivers):
        caxe.scatter([], [], marker=next(marker), label=label)
    caxe.legend(ncol=3, numpoints=1, loc="lower center", mode="expand")
    caxe.set_axis_off()


def finish_plot(fig, pltdir, name, pdf):
    """ Save actual figure to either png or pdf
    Parameters
    ----------
    fig : obj
        actual figure
    pltdir : str
        target directory to store plots
    name : str
        filename for png output without extension
    pdf : obj
        pdf object collection all pages in case of pdf output
    """
    import matplotlib.pyplot as plt
    if '-bias' in name:
        plt.tight_layout()
    if pdf is None:
        filepath = os.path.join(pltdir, name + ".png")
        fig.savefig(filepath)
    else:
        fig.savefig(pdf, dpi=80, format='pdf')
        plt.close()


def make_catchment_plots(cfg, plotdata, catch_info, reference):
    """ Plot catchment averages for precipitation, evaporation
        runoff as bar plots and relation of derived quantities.
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    plotdata : dict
        Dictionary containing the catchment averages
    catch_info : object
        Object containing catchment names, IDs, and reference data
    reference : str
        String containing name of the reference dataset
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Get colorscheme from recipe
    colorscheme = cfg.get('colorscheme', 'default')
    plt.style.use(colorscheme)
    pltdir = cfg[diag.names.PLOT_DIR]
    markerlist = ('s', '+', 'o', '*', 'x', 'D')

    first = list(plotdata.keys())[0]
    for identifier in plotdata[first].keys():
        if cfg.get('output_file_type', 'png') == 'pdf':
            filepath = os.path.join(pltdir, identifier + ".pdf")
            pdf = PdfPages(filepath)
        else:
            pdf = None

        # Compute diagnostics for plots
        expdata, refdata, absdiff, reldiff = {}, {}, {}, {}
        # 1. Variable biases
        for var in plotdata.keys():
            rivers, refdata[var], expdata[var] = get_expdata(
                plotdata[var][identifier], getattr(catch_info, var))
            absdiff[var] = expdata[var] - refdata[var]
            reldiff[var] = absdiff[var] / refdata[var] * 100
            xrivers = range(len(rivers))
        # 2. Coefficients
        prbias = (expdata['pr'] - refdata['pr']) / refdata['pr'] * 100
        rocoef = (expdata['mrro'] / expdata['pr'] * 100) - (
            refdata['mrro'] / refdata['pr'] * 100)
        etcoef = (expdata['evspsbl'] / expdata['pr'] * 100) - (
            refdata['evspsbl'] / refdata['pr'] * 100)

        # Plot diagnostics
        # 1. Barplots for single variables
        title = identifier.upper() + ' vs ' + reference.upper()
        for var in plotdata.keys():
            fig, ax = prep_barplot(title, rivers, var)
            # 1a. Plot absolut bias for every catchment
            ax[0].bar(xrivers, absdiff[var], color="C{}".format(0))
            # 1b. Plot relative bias for every catchment
            ax[1].bar(xrivers, reldiff[var], color="C{}".format(1))
            finish_plot(fig, pltdir, identifier + '_' + var + '-bias', pdf)

        # 2. Runoff coefficient vs relative precipitation bias
        title = identifier.upper() + ' vs ' + reference.upper()
        fig, ax = prep_scatplot(title, 'prbias')
        marker = cycle(markerlist)
        for i, label in enumerate(rivers):
            ax.scatter(prbias[i], rocoef[i], marker=next(marker))
        format_coef_plot(ax)
        add_legend(fig, rivers, markerlist)
        finish_plot(fig, pltdir, identifier + '_pr-vs-ro', pdf)

        # 3. Runoff coefficient vs evaporation coefficient bias
        title = identifier.upper() + ' vs ' + reference.upper()
        fig, ax = prep_scatplot(title, 'etcoef')
        marker = cycle(markerlist)
        for i, label in enumerate(rivers):
            ax.scatter(etcoef[i], rocoef[i], marker=next(marker))
        format_coef_plot(ax)
        add_legend(fig, rivers, markerlist)
        finish_plot(fig, pltdir, identifier + '_et-vs-ro', pdf)

        # Finish pdf if it is the chosen output
        if pdf is not None:
            pdf.close()


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
    Configuration dictionary of the recipe.
    """

    # Get dataset and variable information
    datasets = diag.Datasets(cfg)
    logging.debug("Found datasets in recipe:\n%s", datasets)
    varlist = diag.Variables(cfg)
    logging.debug("Found variables in recipe:\n%s", varlist)

    # Check for correct variables
    if not varlist.vars_available('pr', 'mrro', 'evspsbl'):
        raise ValueError(
            "Diagnostic requires precipitation, runoff and evaporation data")

    # Read catchmentmask
    # to check: Correct way to read auxillary data using recipes?
    catchment_filepath = cfg.get('catchmentmask')
    catchment_cube = iris.load_cube(catchment_filepath)
    if catchment_cube.coord('latitude').bounds is None:
        catchment_cube.coord('latitude').guess_bounds()
    if catchment_cube.coord('longitude').bounds is None:
        catchment_cube.coord('longitude').guess_bounds()
    catchment_areas = iris.analysis.cartography.area_weights(catchment_cube)

    catch_info = defaults()
    reference = 'default'

    # Read data and compute long term means
    # to check: Shouldn't this be part of preprocessing?
    # to check: How to regrid onto catchment_cube grid
    #           with preproc recipe statements
    #           instead of using regrid here?
    allcubes = {}
    plotdata = {}
    for dataset_path in datasets:
        # Prepare data dictionary
        # to check: what is a smart way to do this in python3?
        datainfo = datasets.get_dataset_info(path=dataset_path)
        dset, dexp, dens, dvar = datainfo['dataset'], datainfo[
            'exp'], datainfo['ensemble'], datainfo['short_name']
        if dset not in allcubes.keys():
            allcubes[dset] = []
        # Load data into iris cube
        new_cube = iris.load(dataset_path, varlist.standard_names())[0]
        # Check for expected unit
        if new_cube.units != 'kg m-2 s-1':
            raise ValueError('Unit [kg m-2 s-1] is expected for ',
                             new_cube.long_name.lower(), ' flux')
        # Convert to unit mm per month
        timelist = new_cube.coord('time')
        daypermonth = []
        for mydate in timelist.units.num2date(timelist.points):
            daypermonth.append(
                calendar.monthrange(mydate.year, mydate.month)[1])
        new_cube.data *= 86400.0
        for i, days in enumerate(daypermonth):
            new_cube.data[i] *= days
        # Aggregate over year --> unit mm per year
        year_cube = new_cube.aggregated_by('year', iris.analysis.SUM)
        year_cube.units = "mm a-1"
        # Compute long term mean
        mean_cube = year_cube.collapsed([diag.names.TIME], iris.analysis.MEAN)
        # Regrid to catchment data grid --> maybe use area_weighted instead?
        if mean_cube.coord('latitude').bounds is None:
            mean_cube.coord('latitude').guess_bounds()
        if mean_cube.coord('longitude').bounds is None:
            mean_cube.coord('longitude').guess_bounds()
        # mean_cube_regrid = mean_cube.regrid(catchment_cube,
        #                                     iris.analysis.Linear())
        mean_cube_regrid = mean_cube.regrid(catchment_cube,
                                            iris.analysis.AreaWeighted())
        # Get catchment area means
        rivervalues = {}
        for river, rid in catch_info.catchments.items():
            data_catch = np.ma.masked_where(
                catchment_cube.data.astype(np.int) != rid,
                mean_cube_regrid.data)
            area_catch = np.ma.masked_where(
                catchment_cube.data.astype(np.int) != rid,
                catchment_areas.data)
            rivervalues[river] = (
                data_catch * (area_catch / area_catch.sum())).sum()
        if dset == datainfo.get('reference_dataset', None):
            if reference == 'default':
                reference = datainfo['reference_dataset']
            elif reference != datainfo['reference_dataset']:
                raise ValueError(
                    'Reference must be the same for all variables!')
            setattr(catch_info, dvar, rivervalues)
        else:
            identifier = "_".join([dset.upper(), dexp, dens])
            if dvar not in plotdata.keys():
                plotdata[dvar] = {}
            if identifier in plotdata[dvar].keys():
                raise StandardError(
                    'Variable', dvar,
                    'already exists in plot dictionary --> check script')
            else:
                plotdata[dvar][identifier] = rivervalues

        # Update data for dataset
        # to check: necessary at all? dataset not used later...
        datasets.set_data(mean_cube_regrid.data, dataset_path)
        # Append to cubelist for temporary output
        allcubes[dset].append(mean_cube_regrid)

    # Write regridded and temporal aggregated netCDF data files (one per model)
    # to do: update attributes
    for model in allcubes.keys():
        filepath = os.path.join(cfg[diag.names.WORK_DIR],
                                '_'.join(['postproc', model]) + '.nc')
        if cfg[diag.names.WRITE_NETCDF]:
            iris.save(allcubes[model], filepath)
            logger.info("Writing %s", filepath)

    # Write plotdata as ascii files for user information
    write_plotdata(cfg, plotdata, catch_info, reference)

    # Plot catchment data
    make_catchment_plots(cfg, plotdata, catch_info, reference)


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
