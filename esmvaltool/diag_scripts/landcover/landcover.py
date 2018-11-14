"""Python example diagnostic."""
import logging
import os

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
    allcubes = {}
    for dataset_path in datasets:
        # Get dataset information
        datainfo = datasets.get_dataset_info(path=dataset_path)
        dset = datainfo['dataset']
        if dset not in allcubes.keys():
            allcubes[dset] = []
        # Load data into iris cube
        new_cube = iris.load(dataset_path, varlist.standard_names())[0]
        # Check for expected unit
        if new_cube.units != '%':
            raise ValueError('Unit % is expected for ',
                             new_cube.long_name.lower(), ' area fraction')
        # Compute long term mean
        mean_cube = new_cube.collapsed([diag.names.TIME], iris.analysis.MEAN)
        # Update data for dataset
        datasets.set_data(mean_cube.data, dataset_path)
        # Append to cubelist for temporary output
        allcubes[dset].append(mean_cube)

    # Write regridded and temporal aggregated netCDF data files (one per model)
    # to do: update attributes
    for model in allcubes.keys():
        filepath = os.path.join(cfg[diag.names.WORK_DIR],
                                '_'.join(['postproc', model]) + '.nc')
        if cfg[diag.names.WRITE_NETCDF]:
            iris.save(allcubes[model], filepath)
            logger.info("Writing %s", filepath)


if __name__ == '__main__':

    with diag.run_diagnostic() as config:
        main(config)
