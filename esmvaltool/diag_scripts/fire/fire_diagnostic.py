"""Running a diagnostic about the climate drivers of fire.

Diagnostic is based on:
- Jones et al., State of wildfires (2023-2024), Earth Syst. Sci. Data, 16,
    3601-3685, https://doi.org/10.5194/essd-16-3601-2024, 2024.
- https://github.com/douglask3/Bayesian_fire_models/tree/AR7_REF
    maintainer: Douglas Kelley 
        https://orcid.org/0000-0003-1413-4969 https://github.com/douglask3

"""
import logging
from pathlib import Path
import numpy as np
import glob
import os

import iris

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    group_metadata,
    get_diagnostic_filename,
    save_figure,
)

from esmvaltool.diag_scripts.fire.diagnostic_run_ConFire import (
    diagnostic_run_ConFire
)

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(ancestors):
    record = {
        'ancestors': ancestors
    }
    return record


def setup_basename_file(filename, var, table):
    """Generate the file basename for a newly computed variable using as
    template the filename of an ancestor variable.
    
    Example: If the ancestor filename is 
    CMIP6_MPI-ESM1-2-LR_Amon_historical-ssp370_r1i1p1f1_pr_gn_2010-2020.nc
    then the returned basename will have the following structure
    CMIP6_MPI-ESM1-2-LR_{table}_historical-ssp370_r1i1p1f1_{var}_gn_2010-2020.
    
    Parameters
    ----------
    filename : str
        The full filename to an ancestor variable.
    var : str
        Short name of the variable.
    table : str
        Reference table for the variable (if it exists).
    
    Returns
    -------
    basename : str
        Basename for the file in which the processed variable will be saved.
    """
    # Recover only the filename without path and type then split elements
    file = filename.split('/')[-1].split('.')[0].split('_')
    # Fill in variable and table
    file[2] = table
    file[-3] = var
    return '_'.join(file)


def compute_vpd(config, tas, hurs, provenance):
    """Compute the Vapor Pressure Deficit (VPD) using tas and hurs following
    equations 4 and 5 in:
    Bjarke, N., Barsugli, J. & Livneh, B. Ensemble of CMIP6 derived reference
    and potential evapotranspiration with radiative and advective components.
    Sci Data 10, 417 (2023). https://doi.org/10.1038/s41597-023-02290-0

    The resulting cube is then saved in the work directory of the recipe
    alongside the other variables.
    
    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    tas : iris.cube.Cube
        The cube containing the surface temperature.
    hurs : iris.cube.Cube
        The cube containing the surface relative humidity.
    provenance: dict
        Dictionary containing the attributes from tas/hurs and
        the corresponding files used in the computation.
    """
    # Convert temperature to Celsius for the operation
    tas.convert_units('degrees_C')
    # Compute vpd
    e_s = 0.6108 * np.exp(np.divide(
        17.2694 * tas.data, tas.data + 237.3
    ))
    data = 10 * np.multiply(1 - 0.01 * hurs.data, e_s)  # Convert kPa to hPa
    data = np.maximum(data, 0.)
    # Transfer coordinates
    dim_coords = [
        (coord, idx) for idx, coord in enumerate(
            tas.coords()
        ) if coord.var_name != 'height'
    ]
    # Apply guess_bounds to each coordinate if it doesn't have bounds
    for idx, (coord, _) in enumerate(dim_coords):
        if not coord.has_bounds():
            coord = coord.copy()
            coord.guess_bounds()
            dim_coords[idx] = (coord, idx)
    # Create cube and save it
    vpd = iris.cube.Cube(
        data,
        dim_coords_and_dims=dim_coords,
        long_name='vapor_pressure_deficit',
        var_name='vpd',
        units='hPa'
    )
    vpd.attributes['ancestors'] = provenance['ancestors']
    logger.debug(f'vapor_pressure_deficit iris cube {vpd}')
    basename = setup_basename_file(provenance['ancestors'][0], 'vpd', 'Amon')
    logger.info(
        f'Saving vapor_pressure_deficit in {config['work_dir']}/{basename}'
    )
    # save_data(basename, provenance, config, vpd)
    filename = get_diagnostic_filename(basename, config, extension='nc')
    iris.save(vpd, filename)
    return filename


def main(config):
    """Produce the diagnostics regarding the climate drivers of fire.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")
    logger.info(config)

    # Add diagnostic directory to the config
    diag_dir = os.path.abspath(os.path.dirname(__file__))
    logger.info(f'Diagnostic directory {diag_dir}')
    config['diag_dir'] = diag_dir
    # Add the namelists to the config
    # config['training_namelist'] = f'{diag_dir}/{config_namelist}'
    # config['config_namelist'] = f'{diag_dir}/{config_namelist}'

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info('Processing data for %s', model_dataset)
        logger.info(group)

        vars_file = {}
        for i, attributes in enumerate(group):
            logger.info(f'Variable {attributes['short_name']}')
            vars_file[attributes['short_name']] = attributes
            # Save model information for output plot name
            if i == 0:
                plot_file_info = '_'.join([
                    model_dataset,
                    '_'.join(attributes["activity"]),
                    attributes["exp"],
                    str(attributes["start_year"]),
                    str(attributes["end_year"]),
                ])
        logger.info(vars_file.keys())

        # Compute Vapor Pressure Deficit (VPD)
        logger.info(f'Processing vapor_pressure_deficit for {model_dataset}')
        tas = iris.load_cube(vars_file['tas']['filename'])
        hurs = iris.load_cube(vars_file['hurs']['filename'])
        provenance_record = get_provenance_record(
            [vars_file['tas']['filename'], vars_file['hurs']['filename']]
        )
        filename_vpd = compute_vpd(config, tas, hurs, provenance_record)
        # Log vpd filename for ConFire model run
        vars_file['vpd'] = {}
        vars_file['vpd']['filename'] = filename_vpd

        # Run ConFire model evaluationuld be ordered follo
        # Add the list of input filenames
        # WARNING = they should be ordered following the model run config
        config['files_input'] = [
            [vars_file[v]['filename'], v] for v in config['var_order']
        ]
        logger.info(
            f'Input files used for diagnostic {config['files_input']}'
        )
        logger.info('Running diagnostic model ConFire.')
        figure = diagnostic_run_ConFire(config)

        # Save output figure
        save_figure(
            basename=f"ConFire_results_{plot_file_info}",
            provenance=get_provenance_record(
                [vars_file[v]['filename'] for v in config['var_order'] \
                    if v != 'vpd']
            ),
            cfg=config,
            figure=figure,
            close=True
        )

        # Remove VPD files after diagnostic run = NOT WORKING?
        # only the .nc file is visible, othert files created later?
        if config['remove_vpd_files']:
            logger.info(f'Removing VPD files in {config['work_dir']}')  
            f_not_removed = []
            for f in glob.glob(f'{config['work_dir']}/*vpd*.*'):
                logger.info(f'Removing {f.split('/')[-1]}')
                try:
                    os.remove(f)
                    logger.info(f'Removed {f.split('/')[-1]}')
                except OSError as e:
                    logger.debug(f"Error removing {f.split('/')[-1]}: {e}")
                    f_not_removed.append(f.split('/')[-1])
            logger.info(f'Files not removed: {f_not_removed}')
        
        # Remove ConFire files after diagnostic run
        if config['remove_confire_files']:
            logger.info(
                f'Removing files in {config['work_dir']}/ConFire_outputs'
            )  
            f_not_removed = []
            for f in glob.glob(f'{config['work_dir']}/ConFire_outputs/*.nc'):
                logger.info(f'Removing {f.split('/')[-1]}')
                try:
                    os.remove(f)
                    logger.info(f'Removed {f.split('/')[-1]}')
                except OSError as e:
                    logger.debug(f"Error removing {f.split('/')[-1]}: {e}")
                    f_not_removed.append(f.split('/')[-1])
            logger.info(f'Files not removed: {f_not_removed}')
            if len(f_not_removed) == 0:
                os.rmdir(f'{config['work_dir']}/ConFire_outputs')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)                
