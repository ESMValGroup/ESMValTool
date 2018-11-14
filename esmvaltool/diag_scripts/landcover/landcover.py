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


def write_plotdata(cfg, regnam, modnam, values, datacont):
    """ Output region sums for all datasets of one variable.
    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe
    regname : list
        list containing the region names
    modnam : list
        list containing the dataset names
    value : list
        nested list containing the region sums in 1.0e+6 km2
    datacont : str
        str containing a data identifier used for the file name
    """

    # Write experiment data
    filepath = os.path.join(cfg[diag.names.WORK_DIR], datacont + '.txt')
    ncol = len(regnam)
    with open(filepath, 'w') as f:
        header = '{:25} ' + ncol * ' {:>12}' + '\n'
        body = '{:25} ' + ncol * ' {:12.4f}' + '\n'
        line = [' ',] + [*regnam]
        f.write('Accumulated land coverage for different regions [1.0e+6 km2]\n\n')
        f.write(header.format(*line))
        for ir, row in enumerate(values):
            line = [modnam[ir]] + row
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
    modnam : list
        list containing the dataset names
    value : list
        nested list containing the region sums in 1.0e+6 km2
    var : str
        variable short name
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    outtype = cfg.get('output_file_type', 'png')
    logger.info('Generating plots for filetype: ' + outtype)
    nicename = {'baresoilFrac': 'bare soil covered', 'treeFrac': 'tree covered', 'grassFrac': 'grass covered'}

    # Create pdf in case of pdf output
    if outtype == 'pdf':
        filepath = os.path.join(cfg[diag.names.PLOT_DIR],'coversum_' + var + '.' + outtype)
        pdf = PdfPages(filepath)

    # Loop over metrices
    filepath = os.path.join(cfg[diag.names.PLOT_DIR],'metric_' + var + '.' + outtype)
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False)
    ax.set_title(' '.join(['Accumulated',nicename.get(var,var),'area']))
    ax.set_ylabel(r'Area [$10^6$ km$^2$]')
    nbar, ncat = np.array(values).shape
    index = np.arange(0, (nbar+1)*ncat, nbar+1)
    xticks = np.linspace((nbar+1) / 2.0, (nbar+1)*ncat - (nbar+1) / 2.0, ncat) - 1.0
    ax.set_xticklabels([*regnam])
    ax.set_xticks(xticks)
    for ir, row in enumerate(values):
        ax.bar(index+ir, row)

    fig.subplots_adjust(bottom=0.30)
    caxe = fig.add_axes([0.05, 0.01, 0.9, 0.20])
    for i, label in enumerate(modnam):
        caxe.plot([], [], lw=4, label=label)
    caxe.legend(ncol=2, loc="lower center", mode="expand")
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

    # Read data and compute long term means
    # to check: Shouldn't this be part of preprocessing?
    allcubes = {key: [] for key in varlist.short_names()}
    for dataset_path in datasets:
        # Get dataset information
        datainfo = datasets.get_dataset_info(path=dataset_path)
        ds = datainfo['dataset']
        var = datainfo['short_name']
        # Load data into iris cube
        new_cube = iris.load(dataset_path, varlist.standard_names())[0]
        # Check for expected unit
        if new_cube.units != '%':
            raise ValueError('Unit % is expected for ',
                             new_cube.long_name.lower(), ' area fraction')
        # Compute long term mean
        mean_cube = new_cube.collapsed([diag.names.TIME], iris.analysis.MEAN)
        # Rename variable in cube
        mean_cube._var_name = "_".join([var,ds])
        mean_cube.long_name = " ".join([var,'for dataset',ds])
        # Update data for dataset
        datasets.set_data(mean_cube.data, dataset_path)
        # Append to cubelist for temporary output
        allcubes[var].append(mean_cube)

    # Write regridded and temporal aggregated netCDF data files (one per model)
    # to do: update attributes
    for var in allcubes.keys():
        filepath = os.path.join(cfg[diag.names.WORK_DIR],
                                '_'.join(['postproc', var]) + '.nc')
        if cfg[diag.names.WRITE_NETCDF]:
            iris.save(allcubes[var], filepath)
            logger.info("Writing %s", filepath)

    # Compute land cover fraction areas
    regdef = {'Global': None, 'Tropics': [-30,30], 'North. Hem.': [30,90], 'South. Hem.': [-90,-30]}
    regnam  = regdef.keys()
    for var in allcubes.keys():
        values = []
        modnam = []
        cellarea = iris.analysis.cartography.area_weights(allcubes[var][0], normalize=False)
        for sub_cube in allcubes[var]:
            modnam.append('_'.join(sub_cube._var_name.split('_')[1:]))
            coverarea = sub_cube.copy()
            row = []
            # Compute land cover area in million km2:
            # area = Percentage * 0.01 * area [m2]
            #      / 1.0e+6 [km2]
            #      / 1.0e+6 [1.0e+6 km2]
            coverarea.data *= (0.01 * cellarea / 1.0E+6 / 1.0e+6)
            # Sum over area for different regions
            for reg in regnam:
                if regdef[reg] is not None:
                    zone = iris.Constraint(latitude=lambda v: regdef[reg][0] <= v <= regdef[reg][1])
                    row.append(coverarea.extract(zone).collapsed(['longitude', 'latitude'], iris.analysis.SUM).data.tolist())

                else:
                    row.append(coverarea.collapsed(['longitude', 'latitude'], iris.analysis.SUM).data.tolist())
            values.append(row)

        # Write plotdata as ascii files for user information
        write_plotdata(cfg, regnam, modnam, values, '_'.join(['area',var]))

        # Plot area values
        make_landcover_bars(cfg, regnam, modnam, values, var)


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
