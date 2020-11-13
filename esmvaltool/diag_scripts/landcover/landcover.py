#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Landcover analysis plots.

###############################################################
landcover/landcover.py
Authors ESMValToolV1 Version
    Stefan Hagemann (stefan.hagemann@hzg.de)
    Alexander Loew
    Benjamin Mueller (b.mueller@iggf.geo.uni-muenchen.de)
Port to ESMValTool Version 2
    Tobias Stacke (tobias.stacke@mpimet.mpg.de)
###############################################################

Description
-----------
    Computes accumulated and fractional extent for major land
    cover types (bare soil, crops, grasses, shrubs and trees)
    for the whole globe as well as separated into regions
    (tropics, northern extratropics and southern extratropics).
    The fractions are compared to ESA-CCI land cover data.

    ESA-CCI land cover data needs to be downloaded separately
    by the user and converted to netCDF files containing the
    grid cell fractions for the respective cover type.
    The data and a conversion tool are available at
    https://maps.elie.ucl.ac.be/CCI/viewer/ upon registration.
    Detailed instructions for the installation and use of the
    CCI-LC user tools is available on the same page.

    Note, that all experiments will be regridded onto the
    grid of the ESA-CCI data, thus it is recommended to
    download the coarses resolution which is sufficient for
    the planned study. For testing, ESA-CCI data on 0.5 degree
    resolution was used.
"""

import logging
import os
import numpy as np

import iris
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import esmvaltool.diag_scripts.shared as diag
from esmvaltool.diag_scripts.shared import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


def write_plotdata(infos, modnam, values):
    """Write region values for all datasets of one variable.

    Parameters
    ----------
    infos : list
        contains infos about configuration, regions and provenance
    modnam : dict
        containing list of dataset names for specific metrics
    values : dict
        dictionary of nested list containing the keys
        area --> region sums in 1.0e+6 km2
        frac --> region average fractions in %
    """
    cfg, regnam, prov_rec, var = infos
    # Header information for different metrics
    filehead = {
        'area':
        'Accumulated land coverage for ' + var +
        ' in different regions [1.0e+6 km2]',
        'frac':
        'Average land cover fraction for ' + var + ' in different regions [%]',
        'bias':
        'Bias in average land cover fraction for ' + var +
        ' compared to reference [%]'
    }
    # Write experiment data
    for metric in values.keys():
        filepath = os.path.join(cfg[diag.names.WORK_DIR],
                                '_'.join([metric, var]) + '.txt')
        ncol = len(regnam)
        with open(filepath, 'w') as fout:
            header = '{:35} ' + ncol * ' {:>12}' + '\n'
            body = '{:35} ' + ncol * ' {:12.4f}' + '\n'
            line = [
                ' ',
            ] + regnam
            fout.write(filehead[metric] + '\n\n')
            fout.write(header.format(*line))
            for irow, row in enumerate(values[metric]):
                line = [modnam[metric][irow]] + row
                fout.write(body.format(*line))

    # provenance tracking, only if comparison == variable
    if prov_rec is not None:
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(filepath, prov_rec[var])


def init_plot(cfg, var):
    """Prepare plot and set defaults.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    var : str
        variable short name
    """
    if cfg.get('output_file_type', 'png') == 'pdf':
        filepath = os.path.join(cfg[diag.names.PLOT_DIR],
                                '_'.join(['metrics', var]) + ".pdf")
        pdf = PdfPages(filepath)
    else:
        pdf = None

    nicename = {
        'baresoilFrac': 'bare soil covered',
        'treeFrac': 'tree covered',
        'grassFrac': 'grass covered',
        'cropFrac': 'crop covered',
        'shrubFrac': 'shrub covered'
    }

    info = {
        # Plot titles
        'pt': {
            'area': ' '.join(['Accumulated',
                              nicename.get(var, var), 'area']),
            'frac': ' '.join(['Average',
                              nicename.get(var, var), 'fraction']),
            'bias':
            ' '.join(['Average',
                      nicename.get(var, var), 'fraction bias'])
        },
        # Labels for y axis
        'yl': {
            'area': r'Area [$10^6$ km$^2$]',
            'frac': r'Fraction [%]',
            'bias': r'Bias [%]'
        },
        # Plot directory
        'pd': cfg[diag.names.PLOT_DIR]
    }

    return pdf, info


def plot_bars(info, metric, data, regnam):
    """Add legend and save plot to either png or pdf.

    Parameters
    ----------
    info : dict
        compilation of plot properties
    metric : str
        plot type [area, fraction or bias]
    data : list
        list of floats for plotting
    regnam : list
        list containing the region names
    """
    fig, axs = plt.subplots(nrows=1, ncols=1, sharex=False)
    axs.set_title(info['pt'][metric])
    axs.set_ylabel(info['yl'][metric])
    nbar, ncat = np.array(data).shape
    index = np.arange(0, (nbar + 1) * ncat, nbar + 1)
    xticks = np.linspace((nbar + 1) / 2.0,
                         (nbar + 1) * ncat - (nbar + 1) / 2.0, ncat) - 1.0
    axs.set_xticklabels(regnam)
    axs.set_xticks(xticks)
    for irow, row in enumerate(data):
        axs.bar(index + irow, row)

    return fig


def finish_plot(fig, labels, pltdir, name, pdf):
    """Add legend and save plot to either png or pdf.

    Parameters
    ----------
    fig : obj
        actual figure
    labels : list
        list of plot labels
    pltdir : str
        target directory to store plots
    name : str
        filename for png output without extension
    pdf : obj
        pdf object collection all pages in case of pdf output
    """
    fig.subplots_adjust(bottom=0.20)
    caxe = fig.add_axes([0.05, 0.01, 0.9, 0.20])
    for lbl in labels:
        caxe.plot([], [], lw=4, label=lbl)
    caxe.legend(ncol=2, loc="lower center", fontsize='small')
    caxe.set_axis_off()

    if pdf is None:
        filepath = os.path.join(pltdir, name + ".png")
        fig.savefig(filepath)
    else:
        fig.savefig(pdf, dpi=80, format='pdf')
        plt.close()


def make_landcover_bars(cfg, regnam, modnam, values, var):
    """Make bar plots for regional values.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    regnam : list
        list containing the region names
    modnam : dict
        containing list of dataset names for specific metrics
    values : dict
        dictionary of nested list containing the keys
        area --> region sums in 1.0e+6 km2
        frac --> region average fractions in %
    var : str
        variable short name
    """
    # Get colorscheme from recipe
    plt.style.use(cfg.get('colorscheme', 'seaborn'))

    # Set up plot
    pdf, info = init_plot(cfg, var)

    # Loop over metrices
    for metr in values.keys():
        # Plot plot with bars
        fig = plot_bars(info, metr, values[metr], regnam)
        # Add legend and finish plot
        finish_plot(fig, modnam[metr], info['pd'], '_'.join([metr, var]), pdf)

    if pdf is not None:
        pdf.close()


def sel_lats(latlist, bounds):
    """Return subset of latitudes within bounds.

    Parameters
    ----------
    latlist : numpy array
        contains all latitudes for the cube
    bounds : list
        bounds for latitude selection
    """
    subset = []
    for lat in latlist.tolist():
        if min(bounds) < lat < max(bounds):
            subset.append(lat)

    return subset


def get_timmeans(attr, cubes, refset, prov_rec):
    """Return time averaged data cubes.

    Parameters
    ----------
    attr : dict
        contains metadata for dataset.
    cubes : dict
        collection of iris data cubes.
    refset : dict
        reference dataset names for all variables.
    prov_rec : dict
        contains information for provenance tracking.
    """
    # Get dataset information
    var = attr['short_name']
    # Store name of reference data for given variable
    if var not in refset.keys():
        refset[var] = attr.get('reference_dataset', None)
    # Load data into iris cube
    new_cube = iris.load_cube(attr['filename'])
    # Check for expected unit
    if new_cube.units != '%':
        raise ValueError('Unit % is expected for ' +
                         new_cube.long_name.lower() + ' area fraction')
    # Compute long term mean
    mean_cube = new_cube.collapsed([diag.names.TIME], iris.analysis.MEAN)
    # Rename variable in cube
    mean_cube.var_name = "_".join([
        attr.get('cmor_table', ''),
        attr.get('dataset', ''),
        attr.get('exp', ''),
        attr.get('ensemble', '')
    ]).replace('__', '_').strip("_")
    mean_cube.long_name = " ".join([var, 'for dataset', attr['dataset']])
    # Append to cubelist for temporary output
    if attr['dataset'] == refset[var]:
        cubes['ref'][var].append(mean_cube)
    else:
        cubes['exp'][var].append(mean_cube)
    # Add information to provenance record
    if prov_rec[var] == {}:
        caption = ("Mean land cover fraction for {long_name} between "
                   "{start_year} and {end_year} for different datasets".format(
                       **attr))
        prov_rec[var] = {
            'caption': caption,
            'statistics': ['mean'],
            'domains': ['global'],
            'plot_type': 'regional averages',
            'authors': [
                'hagemann_stefan',
                'loew_alexander',
                'mueller_benjamin',
                'stacke_tobias',
            ],
            'references': [
                'acknow_project',
            ],
            'ancestors': [attr['filename']]
        }
    else:
        prov_rec[var]['ancestors'].append(attr['filename'])


def write_data(cfg, cubes, var, prov_rec):
    """Write intermediate datafield for one variable.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.
    cubes : dict
        collection of iris data cubes.
    var : str
        variable short name
    prov_rec : dict
        contains information for provenance tracking.
    """
    # Compile output path
    filepath = os.path.join(cfg[diag.names.WORK_DIR],
                            '_'.join(['postproc', var]) + '.nc')

    # Join cubes in one list with ref being the last entry
    outcubes = cubes['exp'][var] + cubes['ref'][var]
    if cfg[diag.names.WRITE_NETCDF]:
        iris.save(outcubes, filepath)
        logger.info("Writing %s", filepath)

    # provenance tracking
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filepath, prov_rec[var])


def compute_landcover(var, lcdata, cubes):
    """Return aggregated and averaged land cover values.

    Parameters
    ----------
    var : str
        variable short name
    lcdata : dict
        collection of land cover values per region
    cubes : dict
        collection of time averaged iris data cubes.
    """
    # Define regions
    regdef = {
        'Global': None,
        'Tropics': [-30, 30],
        'North. Hem.': [30, 90],
        'South. Hem.': [-90, -30]
    }

    values = {'area': [], 'frac': [], 'bias': []}
    modnam = {'area': [], 'frac': [], 'bias': []}
    # Compute metrices for all datasets of a given variable
    for sub_cube in cubes:
        modnam['area'].append(sub_cube.var_name)
        modnam['frac'].append(sub_cube.var_name)
        cellarea = sub_cube.copy()
        cellarea.name = 'cellarea'
        cellarea.data = iris.analysis.cartography.area_weights(cubes[0])
        row = {'area': [], 'frac': []}
        # Compute land cover area in million km2:
        # area = Percentage * 0.01 * area [m2]
        #      / 1.0e+6 [km2]
        #      / 1.0e+6 [1.0e+6 km2]
        coverarea = sub_cube.copy()
        coverarea.data *= (0.01 * cellarea.data / 1.0E+6 / 1.0e+6)
        # Sum over area for different regions
        for reg in regdef:
            if regdef[reg] is not None:
                zone = iris.Constraint(
                    latitude=sel_lats(
                        sub_cube.coord('latitude').points, regdef[reg]))
                row['area'].append(
                    coverarea.extract(zone).collapsed(
                        ['longitude', 'latitude'],
                        iris.analysis.SUM).data.tolist())
                row['frac'].append(
                    sub_cube.extract(zone).collapsed(
                        ['longitude', 'latitude'],
                        iris.analysis.MEAN,
                        weights=cellarea.extract(zone).data).data.tolist())

            else:
                row['area'].append(
                    coverarea.collapsed(['longitude', 'latitude'],
                                        iris.analysis.SUM).data.tolist())
                row['frac'].append(
                    sub_cube.collapsed(['longitude', 'latitude'],
                                       iris.analysis.MEAN,
                                       weights=cellarea.data).data.tolist())
        values['area'].append(row['area'])
        values['frac'].append(row['frac'])
    # Compute relative bias in average fractions compared to reference
    reffrac = np.array(values['frac'][-1])
    for imod, modfrac in enumerate(values['frac'][:-1]):
        values['bias'].append(
            ((np.array(modfrac) - reffrac) / reffrac * 100.0).tolist())
        modnam['bias'].append(modnam['frac'][imod])

    lcdata[var] = {'values': values, 'groups': modnam}

    return list(regdef.keys())


def focus2model(cfg, lcdata, refset):
    """Resort lcdata for model focus.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.
    lcdata : dict
        collection of land cover values per region
    refset : dict
        reference dataset names for all variables.
    """
    var = diag.Variables(cfg).short_names()[0]
    shuffle = {key: {} for key in lcdata[var]['groups']['area']}
    for dset in shuffle.keys():
        ids = lcdata[var]['groups']['area'].index(dset)
        if refset[var] in dset:
            shuffle[dset] = {
                'groups': {
                    'area': [],
                    'frac': []
                },
                'values': {
                    'area': [],
                    'frac': []
                }
            }
        else:
            shuffle[dset] = {
                'groups': {
                    'area': [],
                    'frac': [],
                    'bias': []
                },
                'values': {
                    'area': [],
                    'frac': [],
                    'bias': []
                }
            }
        for var in sorted(diag.Variables(cfg).short_names()):
            for metric in shuffle[dset]['groups'].keys():
                shuffle[dset]['groups'][metric].append(var)
                shuffle[dset]['values'][metric].append(
                    lcdata[var]['values'][metric][ids])
    lcdata = shuffle


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.
    """
    # Print dataset and variable information
    logging.debug("Found datasets in recipe:\n%s", diag.Datasets(cfg))
    logging.debug("Found variables in recipe:\n%s", diag.Variables(cfg))

    # Get metadata information
    grouped_input_data = diag.group_metadata(
        cfg['input_data'].values(), 'standard_name', sort='dataset')

    # Prepare dictionaries
    timcubes = {
        'exp': {key: []
                for key in diag.Variables(cfg).short_names()},
        'ref': {key: []
                for key in diag.Variables(cfg).short_names()}
    }
    lcdata = {key: {} for key in diag.Variables(cfg).short_names()}
    refset = {}
    prov_rec = {key: {} for key in diag.Variables(cfg).short_names()}

    # Read data and compute long term means
    for standard_name in grouped_input_data:
        for attributes in grouped_input_data[standard_name]:
            get_timmeans(attributes, timcubes, refset, prov_rec)

    for var in diag.Variables(cfg).short_names():
        # Write regridded and temporal aggregated netCDF data files
        write_data(cfg, timcubes, var, prov_rec)
        # Compute aggregated and fraction average land cover
        regnam = compute_landcover(var, lcdata,
                                   timcubes['exp'][var] + timcubes['ref'][var])

    # Reshuffle data if models are the comparison target
    if cfg.get('comparison', 'variable') == 'model':
        focus2model(cfg, lcdata, refset)
        prov_rec = None
    elif cfg.get('comparison', 'variable') != 'variable':
        raise ValueError('Only variable or model are valid comparison targets')

    # Output ascii files and plots
    for target in lcdata.keys():
        # Write plotdata as ascii files for user information
        infos = [cfg, regnam, prov_rec, target]
        write_plotdata(infos, lcdata[target]['groups'],
                       lcdata[target]['values'])

        # Plot area values
        make_landcover_bars(cfg, regnam, lcdata[target]['groups'],
                            lcdata[target]['values'], target)


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
