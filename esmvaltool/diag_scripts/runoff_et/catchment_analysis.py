#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Catchment specific water flux plots.

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
import calendar
import logging
import os
from itertools import cycle

import iris
import matplotlib
import numpy as np

import esmvaltool.diag_scripts.shared as diag

matplotlib.use('Agg')

logger = logging.getLogger(os.path.basename(__file__))


def get_defaults():
    """Return default reference values for predefined catchments.

    The entries are used in the routine analysecatchments. Catchments and
    reference values are specific for the default catchment mask. All reference
    values are given in mm a-1. Precip data is based on WFDEI, runoff is based
    on GRDC, ET is derived as the difference of both. The values are updated
    and differ slightly from the ESMValTool 1 version.
    Dictionary entries are
        catchments
        mrro
        pr
        evspsbl
    """
    defaults = {
        'catchments': {
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
        },
        'mrro': {
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
        },
        'pr': {
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
        },
        'evspsbl': {
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
    }

    return defaults


def format_coef_plot(my_ax):
    """Move axis from border to center, adapts ticks and labels accordingly.

    Parameters
    ----------
    my_ax : object
        plot axis object
    """
    # Add infos to axis
    my_ax.xaxis.set_label_coords(0.5, -0.025)
    my_ax.yaxis.set_label_coords(-0.025, 0.5)
    # Adapt axis range to center zero
    xmax = np.ceil(
        (np.absolute(np.array(my_ax.get_xlim())).max() + 5) / 10.0) * 10 - 5
    my_ax.set_xlim(xmax * -1, xmax)
    ymax = np.ceil(
        (np.absolute(np.array(my_ax.get_ylim())).max() + 5) / 10.0) * 10 - 5
    my_ax.set_ylim(ymax * -1, ymax)
    # remove 0 from y and x axis
    for key in ['x', 'y']:
        ticks = list(getattr(my_ax, 'get_%sticks' % key)())
        try:
            ticks.remove(0)
        except ValueError:
            pass
        getattr(my_ax, 'set_%sticks' % key)(ticks)

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    my_ax.spines['left'].set_position('center')
    my_ax.spines['bottom'].set_position('center')
    # Eliminate upper and right axes
    my_ax.spines['right'].set_color('none')
    my_ax.spines['top'].set_color('none')
    # Show ticks in the left and lower axes only
    my_ax.xaxis.set_ticks_position('bottom')
    my_ax.yaxis.set_ticks_position('left')


def data2file(cfg, filename, title, filedata):
    """Write data dictionary into ascii file.

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
    with open(filepath, 'w') as out:
        out.write(title + '\n\n')
        for river, value in sorted(filedata.items()):
            out.write('{:25} : {:8.2f}\n'.format(river, value))


def write_plotdata(cfg, plotdata, catchments):
    """Write catchment averaged values for all datasets.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    plotdata : dict
        Dictionary containing the catchment averages
    catchments : dict
        Dictionary containing infomation about catchment mask,
        grid cell size, and reference values
    """
    ref_vars = []
    metric = "catchment averages"
    unit = "[mm a-1]"

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
                title = " ".join([catchments['refname'], metric, unit])
                filedata = catchments[var]
                data2file(cfg, filename, title, filedata)
                ref_vars.append(var)


def get_expdata(expdict, refdict):
    """Get list with catchment averages for experiment and reference.

    Parameters
    ----------
    expdict : dict
        the catchment averages experiments dictionary
    refdict : dict
        the catchment averages reference dictionary
    """
    expdata, refdata, rivers = [], [], []
    for riv, ref in sorted(refdict.items()):
        rivers.append(riv)
        refdata.append(ref)
    for riv in rivers:
        expdata.append(expdict[riv])
    return rivers, np.array(refdata), np.array(expdata)


def prep_barplot(title, rivers, var):
    """Prepare barplot.

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

    fig, my_axs = plt.subplots(nrows=1, ncols=2, sharex=False)
    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.35)
    plottitle = ['\nBias for ', '\nRelative bias for ']
    ylabel = [var.upper() + ' [mm a-1]', 'Relative bias [%]']

    for iax, axs in enumerate(my_axs.tolist()):
        axs.set_title(plottitle[iax] + var.upper())
        axs.set_ylabel(ylabel[iax])
        axs.set_xlabel('Catchment')
        axs.set_xticks(range(len(rivers)))
        axs.set_xticklabels((rivers), fontsize='small')
        for tick in axs.get_xticklabels():
            tick.set_rotation(90)
        axs.axhline(c='black', lw=2)

    return fig, my_axs


def prep_scatplot(title, coeftype):
    """Prepare scatterplot for different coefficients.

    Parameters
    ----------
    title : str
        multipanel plot title
    coeftype : str
        string indicting plot type [prbias,etcoef]
    """
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(nrows=1, ncols=1, sharex=False)
    axs.set_title(title)
    axs.set_ylabel('Bias of runoff coefficient [%]')
    if coeftype == 'prbias':
        axs.set_xlabel('Relative bias of precipitation [%]')
    elif coeftype == 'etcoef':
        axs.set_xlabel('Bias of ET coefficient [%]')
    else:
        raise ValueError('Unexpected coefficient combination in prep_scatplot')

    return fig, axs


def add_legend(fig, rivers, markerlist):
    """Add scatter plot legend with separate axis.

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
    for label in rivers:
        caxe.scatter([], [], marker=next(marker), label=label)
    caxe.legend(ncol=3, numpoints=1, loc="lower center", mode="expand")
    caxe.set_axis_off()


def finish_plot(fig, pltdir, name, pdf):
    """Save actual figure to either png or pdf.

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


def make_catchment_plots(cfg, plotdata, catchments):
    """Plot catchment averages for different metrics.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    plotdata : dict
        Dictionary containing the catchment averages
    catchments : dict
        Dictionary containing infomation about catchment mask,
        grid cell size, and reference values
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
                plotdata[var][identifier], catchments[var])
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
        title = identifier.upper() + ' vs ' + catchments['refname'].upper()
        for var in plotdata.keys():
            fig, axs = prep_barplot(title, rivers, var)
            # 1a. Plot absolut bias for every catchment
            axs[0].bar(xrivers, absdiff[var], color="C{}".format(0))
            # 1b. Plot relative bias for every catchment
            axs[1].bar(xrivers, reldiff[var], color="C{}".format(1))
            finish_plot(fig, pltdir, identifier + '_' + var + '-bias', pdf)

        # 2. Runoff coefficient vs relative precipitation bias
        title = identifier.upper() + ' vs ' + catchments['refname'].upper()
        fig, axs = prep_scatplot(title, 'prbias')
        marker = cycle(markerlist)
        for iriv in range(len(rivers)):
            axs.scatter(prbias[iriv], rocoef[iriv], marker=next(marker))
        format_coef_plot(axs)
        add_legend(fig, rivers, markerlist)
        finish_plot(fig, pltdir, identifier + '_pr-vs-ro', pdf)

        # 3. Runoff coefficient vs evaporation coefficient bias
        title = identifier.upper() + ' vs ' + catchments['refname'].upper()
        fig, axs = prep_scatplot(title, 'etcoef')
        marker = cycle(markerlist)
        for iriv in range(len(rivers)):
            axs.scatter(etcoef[iriv], rocoef[iriv], marker=next(marker))
        format_coef_plot(axs)
        add_legend(fig, rivers, markerlist)
        finish_plot(fig, pltdir, identifier + '_et-vs-ro', pdf)

        # Finish pdf if it is the chosen output
        if pdf is not None:
            pdf.close()


def get_catchment_data(cfg):
    """Read and prepare catchment mask.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    """
    catchments = get_defaults()
    catchments['refname'] = 'default'
    catchment_filepath = cfg.get('catchmentmask')
    catchments['cube'] = iris.load_cube(catchment_filepath)
    if catchments['cube'].coord('latitude').bounds is None:
        catchments['cube'].coord('latitude').guess_bounds()
    if catchments['cube'].coord('longitude').bounds is None:
        catchments['cube'].coord('longitude').guess_bounds()
    catchments['area'] = iris.analysis.cartography.area_weights(catchments['cube'])

    return catchments


def get_sim_data(cfg, datapath, catchment_cube):
    """Read and postprocess netcdf data from experiments.

    Check units, aggregate to long term mean yearly sum and
    regrid to resolution of catchment mask.
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.
    dataset_path : str
        Path to the netcdf file
    catchment_cube : obj
        iris cube object containing simulation data
    """
    datainfo = diag.Datasets(cfg).get_dataset_info(path=datapath)
    identifier = "_".join(
        [datainfo['dataset'].upper(), datainfo['exp'], datainfo['ensemble']])
    # Load data into iris cube
    new_cube = iris.load(datapath, diag.Variables(cfg).standard_names())[0]
    # Check for expected unit
    if new_cube.units != 'kg m-2 s-1':
        raise ValueError('Unit [kg m-2 s-1] is expected for ',
                         new_cube.long_name.lower(), ' flux')
    # Convert to unit mm per month
    timelist = new_cube.coord('time')
    daypermonth = []
    for mydate in timelist.units.num2date(timelist.points):
        daypermonth.append(calendar.monthrange(mydate.year, mydate.month)[1])
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
    m_grid = [iris.analysis.Linear(), iris.analysis.AreaWeighted()]
    mean_cube_regrid = mean_cube.regrid(catchment_cube, m_grid[1])

    return datainfo['short_name'], identifier, mean_cube_regrid


def get_catch_avg(catchments, sim_cube):
    """Compute area weighted averages for river catchments.

    Parameters
    ----------
    catchments : dict
        Dictionary containing infomation about catchment mask,
        grid cell size, and reference values
    sim_cube : obj
        iris cube object containing the simulation data
    """
    avg = {}
    for river, rid in catchments['catchments'].items():
        data_catch = np.ma.masked_where(
            catchments['cube'].data.astype(np.int) != rid, sim_cube.data)
        area_catch = np.ma.masked_where(
            catchments['cube'].data.astype(np.int) != rid, catchments['area'].data)
        avg[river] = (data_catch * (area_catch / area_catch.sum())).sum()
    return avg


def update_reference(catchments, model, rivervalues, var):
    """Update reference catchment averages.

    Parameters
    ----------
    catchments : dict
        Dictionary containing infomation about catchment mask,
        grid cell size, and reference values
    model : str
        name of the data set
    rivervalues : dict
        dictionary of river catchment averages
    var : str
        short name of the variable
    """
    if catchments['refname'] != model and catchments['refname'] != 'default':
        raise ValueError('Reference must be the same for all variables!')
    catchments[var] = rivervalues
    catchments['refname'] = model


def update_plotdata(identifier, plotdata, rivervalues, var):
    """Update simulation catchment averages.

    Parameters
    ----------
    identifier : str
        string consisting of dataset, experiment and ensemble information
    plotdata : dict
        river catchment averages for different variables and datasets
    rivervalues : dict
        river catchment averages for different variables
    var : str
        short name of the variable
    """
    if var not in plotdata.keys():
        plotdata[var] = {}
    if identifier in plotdata[var].keys():
        raise ValueError('Variable', var, 'already exists in plot dict')
    else:
        plotdata[var][identifier] = rivervalues


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.
    """
    # Get dataset and variable information
    logging.debug("Found datasets in recipe:\n%s", diag.Datasets(cfg))
    logging.debug("Found variables in recipe:\n%s", diag.Variables(cfg))

    # Check for correct variables
    if not diag.Variables(cfg).vars_available('pr', 'mrro', 'evspsbl'):
        raise ValueError(
            "Diagnostic requires precipitation, runoff and evaporation data")

    # Read catchmentmask
    # to check: Correct way to read auxillary data using recipes?
    my_catch = get_catchment_data(cfg)

    # Read data, convert units and compute long term means
    # to check: Shouldn't this be part of preprocessing?
    # to check: How to regrid onto catchment_cube grid
    #           with preproc recipe statements
    #           instead of using regrid here?
    allcubes = {}
    plotdata = {}
    for datapath in diag.Datasets(cfg):
        # Get simulation data
        var, identifier, cube = get_sim_data(cfg, datapath, my_catch['cube'])
        # Get river catchment averages
        rivervalues = get_catch_avg(my_catch, cube)
        # Sort into data dictionaries
        datainfo = diag.Datasets(cfg).get_dataset_info(path=datapath)
        model = datainfo['dataset']
        if model == datainfo.get('reference_dataset', None):
            update_reference(my_catch, model, rivervalues, var)
        else:
            update_plotdata(identifier, plotdata, rivervalues, var)

        # Append to cubelist for temporary output
        if model not in allcubes.keys():
            allcubes[model] = []
        allcubes[model].append(cube)

    # Write regridded and temporal aggregated netCDF data files (one per model)
    # to do: update attributes, something fishy with unlimited dimension
    for model, mcube in allcubes.items():
        filepath = os.path.join(cfg[diag.names.WORK_DIR],
                                '_'.join(['postproc', model]) + '.nc')
        if cfg[diag.names.WRITE_NETCDF]:
            iris.save(mcube, filepath)
            logger.info("Writing %s", filepath)

    # Write plotdata as ascii files for user information
    write_plotdata(cfg, plotdata, my_catch)

    # Plot catchment data
    make_catchment_plots(cfg, plotdata, my_catch)


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
