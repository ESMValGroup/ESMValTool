'''
Main diagnostic script for the computation of climate change hotspots.
A comparison between the global and hotspot region tas and pr can be done
Output:
 - scatter plots relating large vs regional scale changes
 - fields of the hotspot for DJF and JJA in CMIP5 and CMIP6
'''
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import iris
from core._ploting_functions import (
    hotspot_fields_plot,
    timeseries_scatter_plot,
)
from esmvalcore.preprocessor import (
    annual_statistics,
    anomalies,
    area_statistics,
    climate_statistics,
    extract_region,
    extract_season,
    extract_time,
    mask_landsea,
    seasonal_statistics,
)

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    names,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


class HotspotDiag(object):
    """class that reads, postprocesses and calls the plotting functions
    necessari to obtain the hotspot figures."""
    def __init__(self, config):
        """config is a dictionary containing metadata regarding input files and
        overall, as the name suggests, configuration options."""
        self.cfg = config
        self.baseline_period = self.cfg["baseline_period"]
        self.scenario = {}
        self.variables = ["tas", "pr"]
        self.periods = ["2041-2060", "2081-2100"]
        self.seasons = ["djf", "jja", "annual"]

    def main(self):
        self.retrieve_data()
        self.compute_10yr_means()
        self.compute_hotspot_fields()

        # delete unused cubes
        for exp in self.experiment_dict:
            for position, ls in enumerate(self.experiment_dict[exp]):
                del self.experiment_dict[exp][position]["mediterranean"]
                del self.experiment_dict[exp][position]["large_scale"]
                logger.info(f"Experiment {exp} for {ls['short_name']}: {ls}")

        # Scatter plots
        for season in self.seasons:
            for var_combination in ["tas:tas", "pr:pr", "pr:tas"]:
                timeseries_scatter_plot(self.experiment_dict, season,
                                        var_combination, self.cfg, names)

        # Hotspot field plots
        for scen in ["26", "45", "85"]:
            for variable in self.variables:
                plot_input = {
                    key: var_meta
                    for key, exp_meta in self.experiment_dict.items()
                    for var_meta in exp_meta if
                    var_meta["short_name"] == variable and key.endswith(scen)
                }
                logger.info(f"plotting {variable.upper()} for scenario {scen}")
                logger.info(f"{plot_input}")
                # plot hostpot field results
                hotspot_fields_plot(
                    plot_input,
                    self.cfg,
                    names,
                    {
                        "variables":
                        self.variables,
                        "seasons":
                        self.seasons,
                        "periods":
                        self.periods,
                        "baseline":
                        f"{self.baseline_period[0]}-{self.baseline_period[1]}",
                    },
                )

    def retrieve_data(self):
        """read data from the preprocessor and extract season and regions."""
        # self.experiment_dict contains the datasets metadata
        # and will contain the computed results.
        self.experiment_dict = group_metadata(self.cfg["input_data"].values(),
                                              "exp")
        process_dict = {
        }  # dict where parallel task results are going to be stored
        # Start parallel task pool
        with ProcessPoolExecutor() as executor:
            for exp in self.experiment_dict:
                for position, dataset_meta in enumerate(
                        self.experiment_dict[exp]):
                    # Gather metadata
                    variable_file = dataset_meta["filename"]
                    project = dataset_meta["project"]
                    variable = dataset_meta["short_name"]
                    logging.info(f"Reading {variable} for {project} "
                                 f"ensemble and {exp} scenario.")
                    # Submit self.load_and_extract_regions
                    # method as a parallel task
                    process_dict[executor.submit(
                        self.load_and_extract_regions, variable_file,
                        variable)] = f"{exp}..{position}"  # Result identifier
        # Gather parallel task results once completed
        for future in as_completed(process_dict):
            identifier = process_dict[future]
            exp = identifier.split("..")[0]
            position = int(identifier.split("..")[-1])
            (
                self.experiment_dict[exp][position]["mediterranean"],
                self.experiment_dict[exp][position]["large_scale"],
            ) = future.result()
        logger.info("Data retrieved")

    def load_and_extract_regions(self, variable_file, variable):
        """Load the cube from variable_file and extract the large scale and
        mediterranean regions."""
        cube = iris.load_cube(variable_file)
        mediterranean_cube, large_scale_cube = self.extract_regions(
            cube, variable)
        mediterranean_cube = self.extract_seasons(mediterranean_cube, variable)
        large_scale_cube = self.extract_seasons(large_scale_cube, variable)
        return [mediterranean_cube, large_scale_cube]

    def extract_regions(self, cube, variable):
        """extract mediterranean and large scale region from cube.

        - If the variable is pr, the 30ºN-45ºN latitudinal band
        is used as large scale region.
        - If the variable is tas, the whole globe is used.

        A sea mask is applied to the Mediterranean region
        """
        if variable == "pr":
            large_scale_cube = extract_region(cube, -180, 180,
                                              self.cfg["region"][2],
                                              self.cfg["region"][3])
        elif variable == "tas":
            large_scale_cube = cube
        mediterranean_cube = extract_region(
            cube,
            self.cfg["region"][0],  # start_lon
            self.cfg["region"][1],  # end_lon
            self.cfg["region"][2],  # start_lat
            self.cfg["region"][3],  # end_lat
        )
        mediterranean_cube = mask_landsea(mediterranean_cube, "sea", True)
        return mediterranean_cube, large_scale_cube

    def extract_seasons(self, cube, variable):
        """The seasons listed in self.seasons are extracted."""
        ensemble = {}
        for season in self.seasons:
            if season == "annual":
                ensemble["annual"] = annual_statistics(cube)
            else:
                for coord in [names.TIME, names.LAT, names.LON]:
                    if not cube.coord(coord).has_bounds():
                        cube.coord(coord).guess_bounds()
                cube = seasonal_statistics(cube)
                ensemble[season] = extract_season(cube, season)
        return ensemble

    def compute_hotspot_fields(self):
        """Compute the hotspot fields with the cubes obtained from
        self.retireve_data.

        the resulting cubes are saved inside self.experiment_dict
        """

        for exp in self.experiment_dict:
            for position, dataset_meta in enumerate(self.experiment_dict[exp]):
                logger.info(f"Computing Med and large scale changes for "
                            f"{dataset_meta['short_name']} "
                            f"{dataset_meta['project']} {exp}")
                self.compute_changes_field(exp, position)
                logger.info(
                    f"Computing hotspot for {dataset_meta['short_name']} "
                    f"{dataset_meta['project']} {exp}")
                logger.info(self.experiment_dict[exp][position])

    def compute_changes_field(self, exp, position):
        """Compute anomaly/change fields with respect to the baseline
        period."""
        variable = self.experiment_dict[exp][position]["short_name"]
        for period in self.periods:
            self.experiment_dict[exp][position][f"hs_field_{period}"] = {}
        med_dict = self.experiment_dict[exp][position]["mediterranean"].items()
        large_dict = self.experiment_dict[exp][position]["large_scale"].items()
        for (season, cube_med), (_, cube_large) in zip(med_dict, large_dict):
            anomaly_med = self.compute_anomaly(cube_med)
            anomaly_large = self.compute_anomaly(cube_large)
            if variable == "pr":
                anomaly_med = self.relative_change(anomaly_med, cube_med)
                anomaly_large = self.relative_change(anomaly_large, cube_large)
            hotspot = anomaly_med - area_statistics(anomaly_large, "mean")
            for period in self.periods:
                start_year = int(period.split("-")[0])
                end_year = int(period.split("-")[-1])
                period_anomaly = extract_time(hotspot, start_year, 1, 1,
                                              end_year, 12, 31)
                period_anomaly = climate_statistics(period_anomaly)
                self.experiment_dict[exp][position][f"hs_field_{period}"][
                    season] = period_anomaly

    def compute_10yr_means(self):
        """Compute the 10yr rolling mean anomaly/change timeseries with the
        cubes obtained from self.retireve_data.

        the resulting cubes are saved inside self.experiment_dict
        """
        for exp in self.experiment_dict:
            for position, dataset_meta in enumerate(self.experiment_dict[exp]):

                for region in ["mediterranean", "large_scale"]:
                    self.log_message(region, dataset_meta, exp)
                    self.compute_timeseries(exp, position, region)

    def compute_timeseries(self, exp, position, region):
        """Compute anomaly/change fields with respect to the baseline
        period."""
        variable = self.experiment_dict[exp][position]["short_name"]
        self.experiment_dict[exp][position][f"{region}_timeseries"] = {}
        for season, cube in self.experiment_dict[exp][position][region].items(
        ):
            area_mean = area_statistics(cube, 'mean')
            anomaly = self.compute_anomaly(area_mean)
            if variable == "pr":
                anomaly = self.relative_change(anomaly, area_mean)
            self.experiment_dict[exp][position][f"{region}_timeseries"][
                season] = anomaly.rolling_window("time", iris.analysis.MEAN,
                                                 10)

    def compute_anomaly(self, cube):
        """Compute anomaly of the cube with respect to the baseline period.

        If the variable is precipitation, the % anomaly is computed with
        respect to the baseline climatology.
        """
        anomaly_reference = {
            "start_year": int(self.baseline_period[0]),
            "start_month": 1,
            "start_day": 1,
            "end_year": int(self.baseline_period[1]),
            "end_month": 12,
            "end_day": 31,
        }
        anomaly = anomalies(cube, "full", reference=anomaly_reference)
        return anomaly

    def relative_change(self, anomaly, cube):
        baseline_cube = extract_time(cube, int(self.baseline_period[0]), 1, 1,
                                     int(self.baseline_period[1]), 12, 31)
        baseline_cube = climate_statistics(baseline_cube)
        anomaly /= baseline_cube
        return anomaly * 100

    @staticmethod
    def log_message(region, dataset_meta, exp):
        logger.info(f"Computing {region.capitalize()} timeseries for "
                    f"{dataset_meta['short_name']} "
                    f"{dataset_meta['project']} {exp}")


def main():
    with run_diagnostic() as config:
        HotspotDiag(config).main()


if __name__ == "__main__":
    main()
