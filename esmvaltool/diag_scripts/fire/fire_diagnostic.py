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

from __future__ import annotations

import logging
from pathlib import Path

import dask.array as da
import iris
import requests

from esmvaltool.diag_scripts.fire.diagnostic_run_confire import (
    diagnostic_run_confire,
    get_provenance_record,
)
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


def setup_basename_file(
    filename: str,
    var: str,
    table: str,
) -> str:
    """Generate the file basename for a newly computed variable.

    the filename of an ancestor variable is used as template.

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
    file = Path(filename).stem.split("_")
    # Fill in variable and table
    file[2] = table
    file[-3] = var
    return "_".join(file)


def download_files_from_zenodo(
    zenodo_link: str,
    output_dir: str,
) -> None:
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
    """
    # Extract the record ID from the Zenodo link
    record_id = zenodo_link.split("/")[-1]
    try:
        int(record_id)
    except ValueError as err:
        msg = f"Extracted Zenodo id '{record_id}' is not a valid integer."
        raise ValueError(msg) from err
    api_url = f"https://zenodo.org/api/records/{record_id}"
    status_request_success = 200

    # Fetch the record metadata
    response = requests.get(api_url, timeout=60)
    if response.status_code != status_request_success:
        msg = "Failed to fetch Zenodo record metadata."
        raise ValueError(msg)

    record_metadata = response.json()

    # Download relevant files in the record
    # - trace.nc
    # - none_trace-params.txt
    # - scalers.csv
    files_to_download = [
        ("trace", "nc"),
        ("none_trace-params", "txt"),
        ("scalers", "csv"),
    ]
    for file_info in record_metadata["files"]:
        file_url = file_info["links"]["self"]
        file_name = file_info["key"]
        if any(
            (f[0] in file_name) and (f[1] == file_name.split(".")[-1])
            for f in files_to_download
        ):
            file_path = Path(output_dir) / file_name
            file_response = requests.get(file_url, timeout=60)
            if file_response.status_code == status_request_success:
                with Path.open(file_path, mode="wb") as file:
                    file.write(file_response.content)
                logger.info("Downloaded %s to %s.", file_name, output_dir)
            else:
                msg = f"Failed to download {file_name} to {output_dir}."
                raise ValueError(msg)


def get_parameter_directory(cfg: dict) -> dict:
    """Get the path to the parameter directory.

    The parameter directory can either consist of:
        - an existing directory which must contain the necessary files
        - a Zenodo URL link to an existing archive from which the necessary
        files will be downloaded from.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary from the recipe.

    Returns
    -------
    param_dir : str
        Path to the directory containing the relevant files.
    """
    param_dir = None
    # If confire_param is a directory = return this directory
    if Path(cfg["confire_param"]).is_dir():
        param_dir = str(cfg["confire_param"])
        if param_dir[-1] != "/":
            param_dir += "/"
        msg = f"Using ConFire model parameter files from directory {param_dir}"
        logger.info(msg)
        # Check if the needed files are present in the directory
        # trace.nc file
        file_trace = list(Path(param_dir).glob("trace*.nc"))
        if not file_trace:
            msg = f"{param_dir} should contain a trace file named trace*.nc."
            raise ValueError(msg)
        # none_trace-params.txt file
        file_none_trace = list(Path(param_dir).glob("none_trace-params*.txt"))
        if not file_none_trace:
            msg = (
                f"{param_dir} should contain a trace parameter file named "
                "none_trace-params*.txt."
            )
            raise ValueError(msg)
        # scalers.csv file
        scale_file = list(Path(param_dir).glob("scalers*.csv"))
        if not scale_file:
            msg = (
                f"{param_dir} should contain a scalers file named "
                "scalers*.csv."
            )
            raise ValueError(msg)
    # If confire_param is a Zenodo URL = download relevant files
    # Zenodo URL should follow the layout:
    # https://zenodo.org/records/{record_id}
    # where record_id is a number typically corresponding to
    # the attribute zenodo.{record_id} in the DOI.
    elif "https://zenodo.org/records/" in cfg["confire_param"]:
        msg = (
            "Retrieving ConFire model parameter files from "
            + cfg["confire_param"]
        )
        logger.info(msg)
        param_dir = cfg["work_dir"] + "/ConFire_parameter_files/"
        Path.mkdir(param_dir, parents=True, exist_ok=True)
        download_files_from_zenodo(
            zenodo_link=cfg["confire_param"],
            output_dir=param_dir,
        )
    # Otherwise raise an error
    else:
        msg = (
            "Recipe input confire_param should either be a directory or a "
            "Zenodo URL. However, parameter is set to "
            f"{cfg['confire_param']}."
        )
        raise ValueError(msg)
    return param_dir


def compute_vpd(
    cfg: dict,
    tas: iris.cube.Cube,
    hurs: iris.cube.Cube,
    provenance: dict,
) -> str:
    """Compute the Vapor Pressure Deficit (VPD).

    This function uses tas and hurs following equations 4 and 5 in:
    Bjarke, N., Barsugli, J. & Livneh, B. Ensemble of CMIP6 derived reference
    and potential evapotranspiration with radiative and advective components.
    Sci Data 10, 417 (2023). https://doi.org/10.1038/s41597-023-02290-0

    The resulting cube is then saved in the work directory of the recipe
    alongside the other variables.

    Parameters
    ----------
    cfg : dict
        The ESMValTool configuration.
    tas : iris.cube.Cube
        The cube containing the surface temperature.
    hurs : iris.cube.Cube
        The cube containing the surface relative humidity.
    provenance: dict
        Dictionary containing the attributes from tas/hurs and
        the corresponding files used in the computation.

    Returns
    -------
    filename: str
        Filename containing the processed VPD output.
    """
    # Convert temperature to Celsius for the operation
    tas.convert_units("degrees_C")
    # Compute vpd
    e_s = 0.6108 * da.exp(
        da.divide(17.2694 * tas.core_data(), tas.core_data() + 237.3),
    )
    # Convert kPa to Pa
    data = 1000.0 * da.multiply(1 - 0.01 * hurs.core_data(), e_s)
    data = da.maximum(data, 0.0)
    # Transfer coordinates
    dim_coords = [
        (coord, idx)
        for idx, coord in enumerate(tas.coords())
        if coord.var_name != "height"
    ]
    # Apply guess_bounds to each coordinate if it doesn't have bounds
    for idx, (coord, _) in enumerate(dim_coords):
        if not coord.has_bounds():
            coord_cp = coord.copy()
            coord_cp.guess_bounds()
            dim_coords[idx] = (coord_cp, idx)
    # Create cube and save it
    vpd = iris.cube.Cube(
        data,
        dim_coords_and_dims=dim_coords,
        long_name="vapor_pressure_deficit",
        var_name="vpd",
        units="Pa",
    )
    vpd.attributes["ancestors"] = provenance["ancestors"]
    debug_vpd = f"vapor_pressure_deficit iris cube {vpd}"
    logger.debug(debug_vpd)
    basename = setup_basename_file(provenance["ancestors"][0], "vpd", "Amon")
    msg = f"Saving vapor_pressure_deficit in {cfg['work_dir']}/{basename}"
    logger.info(msg)
    filename = get_diagnostic_filename(basename, cfg, extension="nc")
    iris.save(vpd, filename)
    if not cfg["remove_vpd_files"]:
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(filename, provenance)
    return filename


def main(cfg: dict) -> None:
    """Produce the diagnostics for the climate drivers of fire.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration from the recipe.
    """
    input_data = cfg["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")
    logger.info(cfg)

    # Add diagnostic directory to the config
    diag_dir = Path(__file__).parent
    cfg["diag_dir"] = diag_dir

    # ConFire model parameter files
    # - user-defined path to files
    # - download files from a Zenodo archive
    # - defaults to the path to files from ESMValTool
    if "confire_param" in cfg:
        if "zenodo" not in cfg["confire_param"]:
            cfg["param_confire"] = (
                Path(cfg["auxiliary_data_dir"]) / cfg["confire_param"]
            )
    else:
        cfg["confire_param"] = "https://zenodo.org/records/14917245"
    cfg["confire_param_dir"] = get_parameter_directory(cfg)

    # Removing or not the computed vapor pressure deficit files
    # default is False
    if "remove_vpd_files" not in cfg:
        cfg["remove_vpd_files"] = False
    # Removing or not the files produced during the ConFire model evaluation
    # default is False
    if "remove_confire_files" not in cfg:
        cfg["remove_confire_files"] = False

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        logger.info(group)
        timerange = group[0]["timerange"]
        project = list({gr["project"] for gr in group})
        experiment = list({gr["exp"] for gr in group})
        plot_file_info = None

        vars_file = {}
        for i, attributes in enumerate(group):
            logger.info("Variable %s", attributes["short_name"])
            vars_file[attributes["short_name"]] = attributes
            # Save model information for output plot name
            if i == 0:
                plot_file_info = "_".join(
                    [
                        model_dataset,
                        attributes["exp"],
                        str(attributes["start_year"]),
                        str(attributes["end_year"]),
                    ],
                )
        logger.info(vars_file.keys())

        if "vpd" in cfg["var_order"]:
            # Compute Vapor Pressure Deficit (VPD)
            logger.info(
                "Processing vapor_pressure_deficit for %s",
                model_dataset,
            )
            tas = iris.load_cube(vars_file["tas"]["filename"])
            hurs = iris.load_cube(vars_file["hurs"]["filename"])
            provenance_record_vpd = get_provenance_record(
                ancestors=[
                    vars_file["tas"]["filename"],
                    vars_file["hurs"]["filename"],
                ],
                model_name=model_dataset,
                project=project,
                experiment=experiment,
                timerange=timerange,
            )
            filename_vpd = compute_vpd(
                cfg,
                tas,
                hurs,
                provenance_record_vpd,
            )
            # Log vpd filename for ConFire model run
            vars_file["vpd"] = {}
            vars_file["vpd"]["filename"] = filename_vpd
            cfg["provenance_record_vpd"] = provenance_record_vpd

        # Run ConFire model evaluation
        # Add the list of input filenames
        # WARNING = they should be ordered following the model run config
        cfg["files_input"] = [
            [vars_file[v]["filename"], v] for v in cfg["var_order"]
        ]
        logger.info(
            "Input files used for diagnostic %s",
            cfg["files_input"],
        )
        cfg["filenames_out"] = [
            "burnt_fraction",
            "fire_weather_control",
            "fuel_load_continuity_control",
        ]
        logger.info("Running diagnostic model ConFire.")
        figures = diagnostic_run_confire(
            cfg,
            model_name=model_dataset,
            timerange=timerange,
            project=project,
            experiment=experiment,
        )

        # Save output figures
        output_file = f"{plot_file_info}"
        ancestors = [
            cfg["files_input"][i][0]
            for i in range(len(cfg["files_input"]))
            if cfg["files_input"][i][1] != "vpd"
        ]
        if "vpd" in cfg["var_order"]:
            ancestors += provenance_record_vpd["ancestors"]
        for i, f in enumerate(cfg["filenames_out"]):
            provenance = get_provenance_record(
                ancestors=ancestors,
                model_name=model_dataset,
                project=group[0]["project"],
                experiment=group[0]["exp"],
                timerange=timerange,
                var=f,
            )
            output_path = get_plot_filename(f"{f}_{output_file}", cfg)
            figures[i].savefig(output_path, bbox_inches="tight", dpi=300)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_path, provenance)

        # Remove or not VPD files after diagnostic run
        if cfg["remove_vpd_files"]:
            logger.info("Removing VPD files in %s", cfg["work_dir"])
            f_not_removed = []
            for f in Path(f"{cfg['work_dir']}/").glob("*vpd*.*"):
                logger.info("Removing %s", f.name)
                try:
                    Path(f).unlink()
                    logger.info("Removed %s", f.name)
                except OSError as e:
                    logger.debug("Error removing %s: %s", f.name, e)
                    f_not_removed.append(f.name)
            logger.info("Files not removed: %s", f_not_removed)

        # Remove or not ConFire files after diagnostic run
        if cfg["remove_confire_files"]:
            logger.info(
                "Removing files in %s/ConFire_outputs",
                cfg["work_dir"],
            )
            f_not_removed = []
            for f in Path(f"{cfg['work_dir']}/ConFire_outputs/").glob("*.nc"):
                logger.info("Removing %s", f.name)
                try:
                    Path(f).unlink()
                    logger.info("Removed %s", f.name)
                except OSError as e:
                    logger.debug("Error removing %s: %s", f.name, e)
                    f_not_removed.append(f.name)
            logger.info("Files not removed: %s", f_not_removed)
            if len(f_not_removed) == 0:
                Path.rmdir(f"{cfg['work_dir']}/ConFire_outputs")


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
