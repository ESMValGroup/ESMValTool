"""Running a diagnostic about the climate drivers of fire.

Diagnostic is based on:
- Jones et al., State of wildfires (2023-2024), Earth Syst. Sci. Data, 16,
    3601-3685, https://doi.org/10.5194/essd-16-3601-2024, 2024.
- https://github.com/douglask3/Bayesian_fire_models/tree/AR7_REF
    maintainer: Douglas Kelley, kelley_douglas
        https://orcid.org/0000-0003-1413-4969
        https://github.com/douglask3

authors:
- lenhardt_julien
- kelley_douglas

"""
import logging
from pathlib import Path
import numpy as np
import glob
import os
import requests

import iris

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    group_metadata,
    get_diagnostic_filename,
    get_plot_filename,
)

from esmvaltool.diag_scripts.fire.diagnostic_run_confire import (
    diagnostic_run_confire
)

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(ancestors):
    """Get record of ancestors.

    Parameters
    ----------
    ancestors : str, list, dict

    Returns
    -------
    record : dict
        Dictionary containing the ancestors.
    """
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


def download_files_from_zenodo(zenodo_link, output_dir):
    """Download files from a Zenodo archive.
    Only the relevant files will be downloaded i.e.:
        - trace.nc
        - none-trace-params.txt
        - scalers.csv
    These files are relevant to run the diagnostic (ConFire model).

    Parameters
    ----------
    zenodo_link : str
        URL link to the Zenodo archive.
    output_dir : str
        Output directory where to store the downloaded files.

    Returns
    -------
    """
    # Extract the record ID from the Zenodo link
    record_id = zenodo_link.split('/')[-1]
    try:
        int(record_id)
    except ValueError:
        raise ValueError(
            f"Extracted Zenodo id '{record_id}' is not a valid integer.")
    api_url = f"https://zenodo.org/api/records/{record_id}"

    # Fetch the record metadata
    response = requests.get(api_url)
    if response.status_code != 200:
        raise ValueError("Failed to fetch Zenodo record metadata.")

    record_metadata = response.json()

    # Download relevant files in the record
    # - trace.nc
    # - none_trace-params.txt
    # - scalers.csv
    files_to_download = [
        ('trace', 'nc'), ('none_trace-params', 'txt'), ('scalers', 'csv')]
    for file_info in record_metadata['files']:
        file_url = file_info['links']['self']
        file_name = file_info['key']
        if any([(f[0] in file_name) and
                (f[1] == file_name.split('.')[-1])
                for f in files_to_download]):
            file_path = os.path.join(output_dir, file_name)
            file_response = requests.get(file_url)
            if file_response.status_code == 200:
                with open(file_path, 'wb') as file:
                    file.write(file_response.content)
                logger.info(f"Downloaded {file_name} to {output_dir}.")
            else:
                raise ValueError(
                    f'Failed to download {file_name} to {output_dir}.')


def get_parameter_directory(config):
    """Get the path to the parameter directory.
    The parameter directory can either consist of:
        - an existing directory which must contain the necessary files
        - a Zenodo URL link to an existing archive from which the necessary
        files will be downloaded from.

    Parameters
    ----------
    config : dict
        ESMValTool configuration dictionary from the recipe.

    Returns
    -------
    param_dir : str
        Path to the directory containing the relevant files.
    """
    param_dir = None
    # If confire_param is a local directory = return this directory
    if os.path.isdir(os.path.expanduser(config['confire_param'])):
        param_dir = os.path.expanduser(config['confire_param'])
        if param_dir[-1] != '/':
            param_dir += '/'
        logger.info(
            'Using ConFire model parameter files from directory ' +
            param_dir
        )
        # Check if the needed files are present in the directory
        # trace.nc file
        file_trace = glob.glob(param_dir + "trace*.nc")
        if not file_trace:
            raise ValueError(
                f'{param_dir} should contain a trace file named trace*.nc.'
            )
        # none_trace-params.txt file
        file_none_trace = glob.glob(param_dir + "none_trace-params*.txt")
        if not file_none_trace:
            raise ValueError(
                f'{param_dir} should contain a trace parameter file named ' +
                'none_trace-params*.txt.'
            )
        # scalers.csv file
        scale_file = glob.glob(param_dir + "scalers*.csv")
        if not scale_file:
            raise ValueError(
                f'{param_dir} should contain a scalers file named ' +
                'scalers*.csv.'
            )
    # If confire_param is a Zenodo URL = download relevant files
    # Zenodo URL should follow the layout:
    # https://zenodo.org/records/{record_id}
    # where record_id is a number typically corresponding to
    # the attribute zenodo.{record_id} in the DOI.
    elif "https://zenodo.org/records/" in config['confire_param']:
        logger.info(
            'Retrieving ConFire model parameter files from ' +
            config['confire_param']
        )
        param_dir = config['work_dir'] + '/ConFire_parameter_files/'
        os.makedirs(param_dir, exist_ok=True)
        download_files_from_zenodo(
            zenodo_link=config['confire_param'],
            output_dir=param_dir
        )
    # Otherwise raise an error
    else:
        raise ValueError(
            'Recipe input confire_param should either be a directory or a ' +
            'Zenodo URL. However, parameter is set to ' +
            f'{config['confire_param']}.'
        )
    return param_dir


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

    Returns
    ----------
    filename: str
        Filename containing the processed VPD output.
    """
    # Convert temperature to Celsius for the operation
    tas.convert_units('degrees_C')
    # Compute vpd
    e_s = 0.6108 * np.exp(np.divide(
        17.2694 * tas.data, tas.data + 237.3
    ))
    # Convert kPa to Pa
    data = 1000. * np.multiply(1 - 0.01 * hurs.data, e_s)
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
        units='Pa'
    )
    vpd.attributes['ancestors'] = provenance['ancestors']
    logger.debug(f'vapor_pressure_deficit iris cube {vpd}')
    basename = setup_basename_file(provenance['ancestors'][0], 'vpd', 'Amon')
    logger.info(
        f'Saving vapor_pressure_deficit in {config['work_dir']}/{basename}'
    )
    filename = get_diagnostic_filename(basename, config, extension='nc')
    iris.save(vpd, filename)
    return filename


def main(config):
    """Produce the diagnostics for the climate drivers of fire.

    Parameters
    ----------
    config : dict
        ESMValTool configuration from the recipe.
    """
    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")
    logger.info(config)

    # Add diagnostic directory to the config
    diag_dir = os.path.abspath(os.path.dirname(__file__))
    logger.info(f'Diagnostic directory {diag_dir}')
    config['diag_dir'] = diag_dir

    # ConFire model parameter files
    # - download files from a Zenodo archive
    # - use a local path to files
    config['confire_param_dir'] = get_parameter_directory(config)

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info('Processing data for %s', model_dataset)
        logger.info(group)
        timerange = group[0]['timerange']

        vars_file = {}
        for i, attributes in enumerate(group):
            logger.info(f'Variable {attributes['short_name']}')
            vars_file[attributes['short_name']] = attributes
            # Save model information for output plot name
            if i == 0:
                plot_file_info = '_'.join([
                    model_dataset,
                    attributes["exp"],
                    str(attributes["start_year"]),
                    str(attributes["end_year"]),
                ])
        logger.info(vars_file.keys())

        if 'vpd' in config['var_order']:
            # Compute Vapor Pressure Deficit (VPD)
            logger.info(
                f'Processing vapor_pressure_deficit for {model_dataset}'
            )
            tas = iris.load_cube(vars_file['tas']['filename'])
            hurs = iris.load_cube(vars_file['hurs']['filename'])
            provenance_record = get_provenance_record(
                [vars_file['tas']['filename'], vars_file['hurs']['filename']]
            )
            filename_vpd = compute_vpd(config, tas, hurs, provenance_record)
            # Log vpd filename for ConFire model run
            vars_file['vpd'] = {}
            vars_file['vpd']['filename'] = filename_vpd

        # Run ConFire model evaluation
        # Add the list of input filenames
        # WARNING = they should be ordered following the model run config
        config['files_input'] = [
            [vars_file[v]['filename'], v] for v in config['var_order']
        ]
        logger.info(
            f'Input files used for diagnostic {config['files_input']}'
        )
        config['filenames_out'] = [
            'burnt_fraction',
            'fire_weather_control',
            'fuel_load_continuity_control',
            'stochastic_control'
        ]
        logger.info('Running diagnostic model ConFire.')
        figures = diagnostic_run_confire(
            config, model_name=model_dataset, timerange=timerange
        )

        # Save output figures
        output_file = f"{plot_file_info}"
        for i, f in enumerate(config['filenames_out']):
            output_path = get_plot_filename(f'{f}_{output_file}', config)
            figures[i].savefig(output_path, bbox_inches='tight', dpi=300)

        # Remove VPD files after diagnostic run
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
