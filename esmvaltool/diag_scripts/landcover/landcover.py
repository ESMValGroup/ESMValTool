"""Python example diagnostic."""
import logging
import os
import pdb

import iris
import esmvaltool.diag_scripts.shared as diag


logger = logging.getLogger(os.path.basename(__file__))


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


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
