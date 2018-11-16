"""Python example diagnostic."""
import logging
import os
import numpy as np
import pdb

import iris
import esmvaltool.diag_scripts.shared as diag
import matplotlib
matplotlib.use('Agg')


logger = logging.getLogger(os.path.basename(__file__))


def write_plotdata(cfg, regnam, modnam, values, var):
    """ Output region sums for all datasets of one variable.
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    regname : list
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

    # Header information for different metrics
    filehead = {'area': 'Accumulated land coverage for '+var+' in different regions [1.0e+6 km2]',
                'frac': 'Average land cover fraction for '+var+' in different regions [%]',
                'bias': 'Bias in average land cover fraction for '+var+' compared to reference [%]'}
    # Write experiment data
    for metric in values.keys():
        filepath = os.path.join(cfg[diag.names.WORK_DIR], '_'.join([metric,var]) + '.txt')
        ncol = len(regnam)
        with open(filepath, 'w') as f:
            header = '{:35} ' + ncol * ' {:>12}' + '\n'
            body = '{:35} ' + ncol * ' {:12.4f}' + '\n'
            line = [' ',] + regnam
            f.write(filehead[metric]+'\n\n')
            f.write(header.format(*line))
            for ir, row in enumerate(values[metric]):
                line = [modnam[metric][ir]] + row
                f.write(body.format(*line))


def make_landcover_bars(cfg, regnam, modnam, values, var):
    """ makes bar plots for the accumulated land cover area
        for the different regions and datasets
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    regname : list
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
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Get colorscheme from recipe
    colorscheme = cfg.get('colorscheme', 'seaborn')
    plt.style.use(colorscheme)

    outtype = cfg.get('output_file_type', 'png')
    logger.info('Generating plots for filetype: ' + outtype)

    nicename = {'baresoilFrac': 'bare soil covered', 'treeFrac': 'tree covered', 'grassFrac': 'grass covered'}
    plottitle = {'area': ' '.join(['Accumulated',nicename.get(var,var),'area']),
                 'frac': ' '.join(['Average',nicename.get(var,var),'fraction']),
                 'bias': ' '.join(['Average',nicename.get(var,var),'fraction bias'])}
    ylabel = {'area': r'Area [$10^6$ km$^2$]',
              'frac': r'Fraction [%]',
              'bias': r'Bias [%]'}

    # Create pdf in case of pdf output
    if outtype == 'pdf':
        filepath = os.path.join(cfg[diag.names.PLOT_DIR],'_'.join(['metrics',var]) + '.' + outtype)
        pdf = PdfPages(filepath)

    # Loop over metrices
    for metric in values.keys():
        filepath = os.path.join(cfg[diag.names.PLOT_DIR],'_'.join([metric,var]) + '.' + outtype)
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False)
        ax.set_title(plottitle[metric])
        ax.set_ylabel(ylabel[metric])
        nbar, ncat = np.array(values[metric]).shape
        index = np.arange(0, (nbar+1)*ncat, nbar+1)
        xticks = np.linspace((nbar+1) / 2.0, (nbar+1)*ncat - (nbar+1) / 2.0, ncat) - 1.0
        ax.set_xticklabels(regnam)
        ax.set_xticks(xticks)
        for ir, row in enumerate(values[metric]):
            ax.bar(index+ir, row)

        fig.subplots_adjust(bottom=0.20)
        caxe = fig.add_axes([0.05, 0.01, 0.9, 0.20])
        for i, label in enumerate(modnam[metric]):
            caxe.plot([], [], lw=4, label=label)
        caxe.legend(ncol=2, loc="lower center", fontsize='small')
        caxe.set_axis_off()

        if outtype == "pdf":
            fig.savefig(pdf, dpi=80, format='pdf')
            plt.close()
        else:
            fig.savefig(filepath)

    if outtype == "pdf":
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

    # Get target for comparison [variable, model]
    comparison = cfg.get('comparison', 'variable')
    if comparison not in ['variable', 'model']:
        raise ValueError('Only variable or model are valid comparison targets')

    # Read data and compute long term means
    # to check: Shouldn't this be part of preprocessing?
    expcubes = {key: [] for key in varlist.short_names()}
    refcubes = {key: [] for key in varlist.short_names()}
    refset = {}
    for dataset_path in sorted(datasets):
        # Get dataset information
        datainfo = datasets.get_dataset_info(path=dataset_path)
        ds = datainfo['dataset']
        var = datainfo['short_name']
        # Store name of reference data for given variable
        if var not in refset.keys():
            refset[var] = datainfo.get('reference_dataset', None)
        # Load data into iris cube
        new_cube = iris.load(dataset_path, varlist.standard_names())[0]
        # Check for expected unit
        if new_cube.units != '%':
            raise ValueError('Unit % is expected for ',
                             new_cube.long_name.lower(), ' area fraction')
        # Compute long term mean
        mean_cube = new_cube.collapsed([diag.names.TIME], iris.analysis.MEAN)
        # Rename variable in cube
        mean_cube._var_name = "_".join([datainfo.get('cmor_table',''),datainfo.get('dataset',''),
                              datainfo.get('exp',''),datainfo.get('ensemble','')]).replace('__','_').strip("_")
        mean_cube.long_name = " ".join([var,'for dataset',ds])
        # Update data for dataset
        datasets.set_data(mean_cube.data, dataset_path)
        # Append to cubelist for temporary output
        if ds == refset[var]:
            refcubes[var].append(mean_cube)
        else:
            expcubes[var].append(mean_cube)

    # Write regridded and temporal aggregated netCDF data files (one per model)
    # to do: update attributes
    allcubes = {key: [] for key in varlist.short_names()}
    for var in expcubes.keys():
        filepath = os.path.join(cfg[diag.names.WORK_DIR],
                                '_'.join(['postproc', var]) + '.nc')
        # Join cubes in one list with ref being the last entry
        allcubes[var] = expcubes[var] + refcubes[var]
        if cfg[diag.names.WRITE_NETCDF]:
            iris.save(allcubes[var], filepath)
            logger.info("Writing %s", filepath)

    # Compute aggregated area and average fractional coverage for different regions
    regdef = {'Global': None, 'Tropics': [-30,30], 'North. Hem.': [30,90], 'South. Hem.': [-90,-30]}
    regnam = list(regdef.keys())
    mydata = {key: {} for key in varlist.short_names()}
    for var in allcubes.keys():
        values = {'area': [], 'frac': [], 'bias': []}
        modnam = {'area': [], 'frac': [], 'bias': []}
        # Compute metrices for all datasets of a given variable
        for sub_cube in allcubes[var]:
            dataset_name = sub_cube._var_name
            modnam['area'].append(dataset_name)
            modnam['frac'].append(dataset_name)
            cellarea = sub_cube.copy()
            cellarea.name = 'cellarea'
            cellarea.data = iris.analysis.cartography.area_weights(allcubes[var][0], normalize=False)
            row = {'area': [], 'frac': []}
            # Compute land cover area in million km2:
            # area = Percentage * 0.01 * area [m2]
            #      / 1.0e+6 [km2]
            #      / 1.0e+6 [1.0e+6 km2]
            coverarea = sub_cube.copy()
            coverarea.data *= (0.01 * cellarea.data / 1.0E+6 / 1.0e+6)
            # Sum over area for different regions
            for reg in regnam:
                if regdef[reg] is not None:
                    zone = iris.Constraint(latitude=lambda v: regdef[reg][0] < v < regdef[reg][1])
                    row['area'].append(coverarea.extract(zone).collapsed(['longitude', 'latitude'], iris.analysis.SUM).data.tolist())
                    row['frac'].append(sub_cube.extract(zone).collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=cellarea.extract(zone).data).data.tolist())

                else:
                    row['area'].append(coverarea.collapsed(['longitude', 'latitude'], iris.analysis.SUM).data.tolist())
                    row['frac'].append(sub_cube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=cellarea.data).data.tolist())
            values['area'].append(row['area'])
            values['frac'].append(row['frac'])
        # Compute relative bias in average fractions compared to reference
        reffrac = np.array(values['frac'][-1])
        for im, modfrac in enumerate(values['frac'][:-1]):
            row = ((np.array(modfrac) - reffrac) /  reffrac * 100.0).tolist()
            values['bias'].append(row)
            modnam['bias'].append(modnam['frac'][im])
        mydata[var] = {'values': values, 'groups': modnam}

    # Reshuffle data if models are the comparison target
    if comparison == 'model':
        shuffle = {key: {} for key in mydata[var]['groups']['area']}
        for ds in shuffle.keys():
            ids = mydata[var]['groups']['area'].index(ds)
            if refset[var] in ds:
                shuffle[ds] = {'groups': {'area': [], 'frac': []},
                               'values': {'area': [], 'frac': []}}
            else:
                shuffle[ds] = {'groups': {'area': [], 'frac': [], 'bias': []},
                               'values': {'area': [], 'frac': [], 'bias': []}}
            for var in sorted(allcubes.keys()):
                for metric in shuffle[ds]['groups'].keys():
                    shuffle[ds]['groups'][metric].append(var)
                    shuffle[ds]['values'][metric].append(mydata[var]['values'][metric][ids])
        mydata = shuffle


    # Output ascii files and plots
    for target in mydata.keys():
        # Write plotdata as ascii files for user information
        write_plotdata(cfg, regnam, mydata[target]['groups'], mydata[target]['values'], target)

        # Plot area values
        make_landcover_bars(cfg, regnam, mydata[target]['groups'], mydata[target]['values'], target)


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
