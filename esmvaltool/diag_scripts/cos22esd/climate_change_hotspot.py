"""Main diagnostic script for the computation of climate change hotspots.

A comparison between the global and hotspot region tas and pr can be done
Output:
 - scatter plots relating large vs regional-scale changes.
 - fields of the hotspot for DJF and JJA in CMIP5 and CMIP6.
"""
import iris
from esmvalcore.preprocessor import (
    annual_statistics,
    anomalies,
    area_statistics,
    climate_statistics,
    extract_region,
    extract_season,
    extract_time,
    mask_landsea,
)

from esmvaltool.diag_scripts.shared import run_diagnostic, save_data


class HotspotDiag:
    """Hotspot diagnostics' class.

    Class that reads, post-processes and calls the plotting functions
    necessary to obtain the hotspot figures.
    """

    def __init__(self, config):
        """Variables definition.

        config is a dictionary containing metadata regarding input files
        and overall, as the name suggests, configuration options.
        """
        self.cfg = config
        self.anomaly_reference = {
            "start_year": int(self.cfg["baseline_period"][0]),
            "start_month": 1,
            "start_day": 1,
            "end_year": int(self.cfg["baseline_period"][1]),
            "end_month": 12,
            "end_day": 31,
        }
        self.regions = {
            "tas": {
                "large_scale": "global",
                "regional": self.cfg["region_name"]
            },
            "pr": {
                "large_scale": "latitudinal belt",
                "regional": self.cfg["region_name"]
            },
        }

    def compute(self):
        """Compute Hotspot diagnostics."""
        input_data_dict = list(self.cfg["input_data"].values())[0]
        filename = input_data_dict['filename']
        large_scale_cube = iris.load_cube(filename)
        regional_cube = self.extract_regional(large_scale_cube)
        if large_scale_cube.var_name == "pr":
            large_scale_cube = extract_region(large_scale_cube, -180, 180,
                                              self.cfg["region"][2],
                                              self.cfg["region"][3])

        for season in ["jja", "djf", "annual"]:
            time_extr_large_scale_cube = self.extract_annual_or_season(
                large_scale_cube, season)
            time_extr_regional_cube = self.extract_annual_or_season(
                regional_cube, season)
            time_subset_cubes = [
                time_extr_large_scale_cube, time_extr_regional_cube
            ]
            self.compute_hotspot_fields(time_subset_cubes, season, filename,
                                        input_data_dict)
            for key, cube in zip(
                    ["large_scale", "regional"],
                    [time_extr_large_scale_cube, time_extr_regional_cube]):
                rolling_mean_cube = self.compute_rolling_means(cube)
                basename = f"rolling_mean_{key}_{season}"
                region = self.regions[regional_cube.var_name][key]
                provenance_record = self.get_rolling_mean_record(
                    input_data_dict, season, region, filename)
                save_data(basename, provenance_record, self.cfg,
                          rolling_mean_cube)

    def extract_regional(self, cube):
        """Extract the hotspot region and mask its sea."""
        regional_cube = extract_region(
            cube,
            self.cfg["region"][0],  # start_lon
            self.cfg["region"][1],  # end_lon
            self.cfg["region"][2],  # start_lat
            self.cfg["region"][3],  # end_lat
        )
        regional_cube = mask_landsea(regional_cube, "sea", True)

        return regional_cube

    @staticmethod
    def extract_annual_or_season(cube, season):
        """Compute the statistics of the cube.

        Either annual, or extract a season (djf, mam, jja or son).
        """
        if season == "annual":
            cube = annual_statistics(cube)
        else:
            cube = extract_season(cube, season)
        return cube

    def compute_hotspot_fields(self, cubes, season, filename, input_data_dict):
        """Compute the hotspot fields with self.retireve_data cubes.

        the resulting cubes are saved inside self.experiment_dict.
        """
        regional_anomaly = anomalies(cubes[1],
                                     "full",
                                     reference=self.anomaly_reference)
        large_scale_anomaly = anomalies(cubes[0],
                                        "full",
                                        reference=self.anomaly_reference)
        if cubes[1].var_name == "pr":
            regional_anomaly = self.relative_change(regional_anomaly, cubes[1])
            large_scale_anomaly = self.relative_change(large_scale_anomaly,
                                                       cubes[0])
        hotspot = regional_anomaly - \
            area_statistics(large_scale_anomaly, "mean")
        hotspot.var_name = f"{regional_anomaly.var_name}_hotspot"
        for period in self.cfg["future_periods"]:
            start_year = int(period.split("-")[0])
            end_year = int(period.split("-")[-1])
            period_anomaly = extract_time(hotspot, start_year, 1, 1, end_year,
                                          12, 31)
            period_anomaly = climate_statistics(period_anomaly)

            basename = f"hotspot_{season}_{period}"
            provenance_record = self.get_hotspot_provenance_record(
                input_data_dict, period, season, filename)
            save_data(basename, provenance_record, self.cfg, period_anomaly)

    def compute_rolling_means(self, cube):
        """Compute the 10yr rolling mean anomaly.

        A timeseries with respect to the baseline period is obtained for
        the input cube.
        """
        area_mean = area_statistics(cube, 'mean')
        anomaly = anomalies(area_mean,
                            "full",
                            reference=self.anomaly_reference)
        if cube.var_name == "pr":
            anomaly = self.relative_change(anomaly, area_mean)
        timeseries_cube = anomaly.rolling_window("time", iris.analysis.MEAN,
                                                 10)
        return timeseries_cube

    def relative_change(self, anomaly, cube):
        """Compute relative anomaly."""
        baseline_cube = extract_time(cube,
                                     int(self.cfg["baseline_period"][0]), 1, 1,
                                     int(self.cfg["baseline_period"][1]), 12,
                                     31)
        baseline_cube = climate_statistics(baseline_cube)
        anomaly /= baseline_cube
        return anomaly * 100

    def get_hotspot_provenance_record(self, attributes, period, season,
                                      ancestor_files):
        """Create a provenance record.

        It describes the hotspot fields diagnostic data.
        """
        baseline = self.cfg["baseline_period"]
        caption = (f"{attributes['project']} {season.upper()} "
                   f"{attributes['long_name']} differences between global "
                   f"and {self.cfg['region_name']} anomalies of period "
                   f"{period}, against the baseline period "
                   f"{baseline[0]}-{baseline[1]}, according to.")

        record = {
            'caption': caption,
            'statistics': ['mean', 'anomaly', 'diff'],
            'domains': ['global', 'reg'],
            'authors': [
                'cos_josep',
            ],
            'references': [
                'cos22esd',
            ],
            'ancestors': [ancestor_files],
        }
        return record

    def get_rolling_mean_record(self, attributes, season, region,
                                ancestor_files):
        """Create a provenance record.

        It describes the rolling mean diagnostic data.
        """
        baseline = self.cfg["baseline_period"]
        caption = (f"{attributes['project']} {season.upper()} 10 year "
                   f"{attributes['long_name']} rolling mean between "
                   f"{attributes['start_year']} and {attributes['end_year']}, "
                   f"against the baseline period {baseline[0]}-{baseline[1]}, "
                   f"according to the {region} mean")

        domain = "reg"
        if region == "global":
            domain = "global"

        record = {
            'caption': caption,
            'statistics': ['anomaly', "other"],
            'domains': [domain],
            'authors': [
                'cos_josep',
            ],
            'references': [
                'cos22esd',
            ],
            'ancestors': [ancestor_files],
        }
        return record


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        HotspotDiag(cfg).compute()
