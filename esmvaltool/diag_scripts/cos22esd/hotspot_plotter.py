"""Functions to plot output from ancestor cos22esd/climate_change_hotspot.py.

The plots produced reproduce Figs. 2, 3, S1, S2, S4 from Cos et al.
2022.
"""
import os

import cartopy.crs as ccrs
import cartopy.feature as cf
import iris
import iris.plot as iplt
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    get_cfg,
    group_metadata,
    io,
    names,
    run_diagnostic,
    save_figure,
)


class HotspotPlot:
    """class that plots the results."""
    def __init__(self, config):
        """Variable definition.

        Config is a dictionary containing metadata regarding input files
        and overall, as the name suggests, configuration options.
        """
        self.cfg = config
        self.var_combinations = ["tas:tas", "pr:pr", "pr:tas"]
        self.seasons = ["jja", "djf", "annual"]
        self.projects = ["cmip5", "cmip6"]
        self.variables = ["tas", "pr"]
        self.scenarios = ["26", "45", "85"]

        # generate list of candidate bound limits
        small = np.arange(0.1, 1, 0.1)
        medium = np.arange(1, 11)
        high = np.arange(20, 100, 10)
        v_high = np.arange(150, 400, 50)
        self.bound_candidates = np.concatenate(
            (small, medium, high, v_high)) * 5 / 4

    def compute(self):
        """Collect datasets and call the plotting functions."""
        # find number of models inside every multi-model mean
        metadata_files = [
            file for file in self.cfg["input_files"]
            if "tas/metadata.yml" in file
        ]
        self.cfg["N"] = {}
        for meta_file in metadata_files:
            n_identifyer = meta_file.split("/tas/")[0].split("/tas_")[-1]
            metadata = group_metadata(get_cfg(meta_file).values(), "dataset")
            self.cfg["N"][n_identifyer] = len(metadata.keys()) - 1
        # call hotspot field plots
        for scenario in self.scenarios:
            fields_dict = {}
            ancestor_files = []
            for filename in io.get_all_ancestor_files(self.cfg,
                                                      pattern='hotspot_*.nc'):
                key = os.path.basename(os.path.dirname(filename))
                splitname = os.path.basename(filename).split("_")
                if key.split("_")[-1] == scenario:
                    fields_dict[(
                        f"{splitname[-1].split('.nc')[0]}_"
                        f"{splitname[1]}_{key}")] = iris.load_cube(filename)
                    ancestor_files.append(filename)
            self.hotspot_fields_plot(fields_dict, scenario, ancestor_files)

        # call scatter plots
        for season in self.seasons:
            timeseries_dict = {"large_scale": {}, "regional": {}}
            ancestors_dict = {"large_scale": {}, "regional": {}}
            for region in timeseries_dict:
                for filename in io.get_all_ancestor_files(
                        self.cfg,
                        pattern=f'rolling_mean_{region}_{season}.nc'):
                    key = os.path.basename(os.path.dirname(filename))
                    timeseries_dict[region][key] = iris.load_cube(filename)
                    ancestors_dict[region][key] = filename
            for var_combination in self.var_combinations:
                self.timeseries_scatter_plot(timeseries_dict, season,
                                             var_combination, ancestors_dict)

    def hotspot_fields_plot(self,
                            results_dict,
                            scenario,
                            ancestor_files_var,
                            tas_bound=None,
                            pr_bound=None):
        """Regional climate change hotspot maps for TAS and PR.

        Local temperature and precipitation change differences
        with respect to the mean global temperature change and the
        mean latitudinal belt precipitation change.
        The changes are computed with respect to the 1986-2005 mean
        for the mid-term and long-term periods.

        The differences are shown for the CMIP5 and CMIP6 winter,
        summer and annual mean projections.

        N indicates the number of models included in the ensemble mean.
        """
        top, bottom, left, right = 0.75, 0.2, 0.02, 0.99
        wspace, hspace = 0.005, 0.005
        sorted_keys = [
            f"{period}_{season}_{variable}_{project}_{scenario}"
            for variable in self.variables
            for period in self.cfg["future_periods"]
            for project in self.projects for season in self.seasons
        ]
        sorted_keys = [
            sorted_keys[:len(sorted_keys) // 2],
            sorted_keys[len(sorted_keys) // 2:]
        ]
        ancestor_files_var = [[
            ancestor_file for ancestor_file in ancestor_files_var
            if f"/{var}_" in ancestor_file
        ] for var in self.variables]
        for ancestor_files, keys, variable in zip(ancestor_files_var,
                                                  sorted_keys, self.variables):
            fig = plt.figure(figsize=(14.4, 3.4),
                             constrained_layout=True,
                             dpi=300)
            plt.gcf().subplots_adjust()
            proj, projection, path_ext, plotextend = self.define_projection(
                self.cfg["region"])
            gspec = GridSpec(
                len(self.cfg["future_periods"]),
                len(self.seasons) * 2,
                figure=fig,
                hspace=hspace,
                wspace=wspace,
                top=top,
                bottom=bottom,
                left=left,
                right=right,
            )

            # bound colorbar to abs(max) value on the map
            # or use input bounds {"tas": [bounds], "pr": [bounds]}
            if variable == "tas":
                if tas_bound:
                    bound_limit = tas_bound
                else:
                    bound_limit = self.find_abs_bound_range(results_dict, keys)
                cmap = plt.cm.RdBu_r
            else:
                if pr_bound:
                    bound_limit = pr_bound
                else:
                    bound_limit = self.find_abs_bound_range(results_dict,
                                                            keys,
                                                            avg_over=25)
                cmap = plt.cm.BrBG
            bounds = np.linspace(-1 * bound_limit, bound_limit, 11)

            # plot each panel
            for i, key in enumerate(keys):
                if i < 6:
                    axes = fig.add_subplot(gspec[0, i], projection=proj)
                    plt.title(f"{key.split('_')[1].upper()}")
                else:
                    axes = fig.add_subplot(gspec[1, i - 6], projection=proj)

                if projection == "LambertConformal":
                    proj_to_data = (
                        ccrs.PlateCarree()._as_mpl_transform(axes) -
                        axes.transData)
                    rect_in_target = proj_to_data.transform_path(path_ext)
                    axes.set_boundary(rect_in_target, use_as_clip_path=True)
                    axes.set_facecolor("silver")
                axes.set_extent(plotextend, crs=ccrs.PlateCarree())
                axes.coastlines("50m", linewidth=0.8)
                axes.add_feature(cf.BORDERS, alpha=0.4)

                norm = colors.BoundaryNorm(boundaries=bounds,
                                           ncolors=256,
                                           extend="both")
                cube = self.regrid_longitude_coord(results_dict[key])
                fill = iplt.pcolormesh(
                    cube,
                    norm=norm,
                    coords=(names.LON, names.LAT),
                    cmap=cmap,
                )

            for p_ind, project in enumerate(self.projects):
                n_models = self.cfg["N"][f"{project}_{scenario}"]
                plt.figtext(
                    left + 0.18 + p_ind * (right - left) / 2,
                    0.85,
                    (f"{self.formatter(project.upper())} "
                     f"{self.formatter(f'{project}-{scenario}')} "
                     f"(N={n_models})"),
                    fontsize="large",
                )
            for row, period in enumerate(self.cfg["future_periods"]):
                ypos = top - (top - bottom) / 2 * (1 + row * 1.1) + 0.05
                plt.figtext(
                    0.005,
                    ypos,
                    period,
                    rotation="vertical",
                    fontsize="11",
                )
            mid = left + (right - left) / 2
            line = plt.Line2D((mid, mid), (bottom, 0.9),
                              color="k",
                              linewidth=1)
            fig.add_artist(line)

            colorbar_axes = plt.gcf().add_axes([0.25, 0.125, 0.5, 0.04])
            cbar = plt.colorbar(fill,
                                colorbar_axes,
                                orientation="horizontal",
                                extend="both")
            if variable == "pr":
                cbar.set_label("%")
            else:
                cbar.set_label(
                    self.formatter(str(results_dict[keys[-1]].units)))
            if variable == "tas":
                against_region = "global"
            elif variable == "pr":
                against_region = (
                    f"{self.cfg['region'][2]}$^o$ N-"
                    f"{self.cfg['region'][3]}$^o$ N latitudinal belt")

            suptitle = (f"{self.cfg['region_name']} {variable.upper()} "
                        f"change against mean {against_region} future "
                        f"climatology. Baseline period: "
                        f"{self.cfg['baseline_period'][0]}-"
                        f"{self.cfg['baseline_period'][1]}")
            plt.suptitle(suptitle, fontsize=13)

            basename = f"{variable}_{scenario}"
            provenance_record = self.get_hotspot_provenance(
                suptitle, scenario, ancestor_files)
            save_figure(basename, provenance_record, self.cfg)

    def timeseries_scatter_plot(self, results_dict, season, var_combination,
                                ancestors_dict):
        """Regional vs large scale changes for three scenarios.

        Computed for different seasons for the CMIP5 and CMIP6 ensemble means.
        Each dot in the plot represents a 10 year mean change beginning from
        1960-1969 (light coloring) until 2091-2100 (opaque coloring).
        The changes are computed with 1986-2005 as baseline.

        An ordinary least squares linear regression is computed and the
        slope and rvalue are shown. N indicates the number of models
        included in the ensemble mean.
        """
        base_colors = {"cmip5": "#2161A6", "cmip6": "#BB3437"}
        legend_elements = {}
        fig = plt.figure(figsize=(12, 4), constrained_layout=True, dpi=300)
        gspec = fig.add_gridspec(1, 3)
        rvalue = {}
        slope = {}
        min_range, max_range = {}, {}
        min_glob, max_glob = [], []
        axes = []
        ancestor_files = []
        for panel, scen in enumerate(self.scenarios):
            legend_elements = {scen: []}
            axes.append(fig.add_subplot(gspec[0, panel]))
            regional_var = var_combination.split(":")[0]
            regional_keys = [
                f"{regional_var}_{project}_{scen}" for project in self.projects
            ]
            large_scale_var = var_combination.split(":")[1]
            large_scale_keys = [
                f"{large_scale_var}_{project}_{scen}"
                for project in self.projects
            ]
            for regional_key, large_scale_key in zip(regional_keys,
                                                     large_scale_keys):
                project = regional_key.split("_")[1]
                ls_cube = results_dict["large_scale"][large_scale_key]
                large_scale_signal_ts = ls_cube.data
                r_cube = results_dict["regional"][regional_key]
                regional_signal_ts = r_cube.data

                res = stats.linregress(large_scale_signal_ts,
                                       regional_signal_ts)
                y_values = res.intercept + res.slope * \
                    np.array(large_scale_signal_ts)
                rvalue[project] = res.rvalue
                slope[project] = res.slope
                timesteps = np.linspace(0, 1, len(large_scale_signal_ts))
                if project == "cmip6":
                    cb_colors = plt.cm.Reds(
                        np.linspace(0, 1, len(regional_signal_ts)))
                if project == "cmip5":
                    cb_colors = plt.cm.Blues(
                        np.linspace(0, 1, len(regional_signal_ts)))
                cb_colors[:, -1] = timesteps

                if var_combination == "pr:tas":
                    max_range[project] = max(regional_signal_ts)
                    min_range[project] = min(regional_signal_ts)
                else:
                    if max(regional_signal_ts) > max(large_scale_signal_ts):
                        max_range[project] = max(regional_signal_ts)
                    else:
                        max_range[project] = max(large_scale_signal_ts)
                    if min(regional_signal_ts) < min(large_scale_signal_ts):
                        min_range[project] = min(regional_signal_ts)
                    else:
                        min_range[project] = min(large_scale_signal_ts)

                title_format = {
                    "26": "RCP2.6/SSP1-2.6",
                    "45": "RCP4.5/SSP2-4.5",
                    "60": "RCP6.0/SSP4-6.0",
                    "85": "RCP8.5/SSP5-8.5",
                }

                axes[panel].scatter(
                    large_scale_signal_ts,
                    regional_signal_ts,
                    facecolors="none",
                    linewidths=0.8,
                    s=70,
                    color=cb_colors,
                    label=self.formatter(project.upper()),
                )
                max_glob.append(max(large_scale_signal_ts))
                min_glob.append(min(large_scale_signal_ts))
                axes[panel].plot(large_scale_signal_ts,
                                 y_values,
                                 color=base_colors[project])
                if len(legend_elements[scen]) < 2:
                    n_models = self.cfg["N"][f"{project}_{scen}"]
                    legend_elements[scen].append(
                        Patch(
                            facecolor=base_colors[project],
                            edgecolor=base_colors[project],
                            label=(f"{self.formatter(project.upper())} "
                                   f"(N={n_models})"),
                        ))

                # collect used ancestor files
                ancestor_files.append(
                    ancestors_dict["large_scale"][regional_key])
                ancestor_files.append(ancestors_dict["regional"][regional_key])

            if var_combination.partition(":")[-1] == "tas":
                against_region = "Global"
            else:
                against_region = (
                    f"{self.cfg['region'][2]}$^o$ N-{self.cfg['region'][3]}"
                    f"$^o$ N latitudinal belt")
            large_scale_units = self.formatter(
                str(results_dict['large_scale'][large_scale_keys[-1]].units))
            regional_units = self.formatter(
                str(results_dict['regional'][regional_keys[-1]].units))
            xlabel = (f"{against_region} "
                      f"{var_combination.partition(':')[-1].upper()} "
                      f"[{large_scale_units}]")
            axes[panel].set_xlabel(xlabel)
            ylabel = (f"{self.cfg['region_name']} "
                      f"{var_combination.partition(':')[0].upper()} "
                      f"[{regional_units}]")
            axes[panel].set_ylabel(ylabel)

            max_lim = max(max_range.values())
            min_lim = min(min_range.values())

            delta_range = max_lim - min_lim
            min_lim -= delta_range * 0.1
            max_lim += delta_range * 0.1

            axes[panel].axvline(
                x=0,
                ymin=-1000,
                ymax=1000,
                color="grey",
                linestyle="dotted",
                alpha=0.6,
            )
            axes[panel].axhline(
                y_values=0,
                xmin=-1000,
                xmax=1000,
                color="grey",
                linestyle="dotted",
                alpha=0.6,
            )
            axes[panel].set_title(
                f"Scenario: {title_format[scen]} \n CMIP5: rval="
                f"{rvalue['cmip5']:.3f}; slope={slope['cmip5']:.3f} "
                f"\n CMIP6: rval={rvalue['cmip6']:.3f}; "
                f"slope={slope['cmip6']:.3f}")
            axes[panel].legend(handles=legend_elements[scen])

            long_name_dict = {"pr": "precipitation", "tas": "temperature"}
            if var_combination == "pr:tas":
                suptitle = (f"{self.cfg['region_name']} {season.upper()} "
                            f"precipitation vs global {season.upper()} "
                            f"temperature.\n 10yr rolling means 1960-2100, "
                            f"Baseline: 1986-2005")
                plt.suptitle(suptitle)
            else:
                y_combination = var_combination.partition(':')[0]
                suptitle = (f"{self.cfg['region_name']} vs {against_region} "
                            f"{season.upper()} "
                            f"{long_name_dict[y_combination]}"
                            f".\n 10yr rolling means 1960-2100, "
                            f"Baseline: 1986-2005")
                plt.suptitle(suptitle)
        for box in range(3):
            axes[box].set_ylim(min_lim, max_lim)
            if var_combination == "pr:tas":
                min_l = min(min_glob) - (max(max_glob) - min(min_glob)) * 0.1
                max_l = max(max_glob) + (max(max_glob) - min(min_glob)) * 0.1
                axes[box].set_xlim(min_l, max_l)
            else:
                axes[box].set_xlim(min_lim, max_lim)

            if (slope["cmip5"] + slope["cmip6"]) >= 0:
                axes[box].plot(
                    [-1000, 1000],
                    [-1000, 1000],
                    color="gray",
                    alpha=0.6,
                )
            else:
                axes[box].plot(
                    [-1000, 1000],
                    [1000, -1000],
                    color="gray",
                    alpha=0.6,
                )
        basename = f"scenario_combination_{var_combination}_{season}"
        provenance_record = self.get_rolling_mean_provenance(
            suptitle, ancestor_files)
        save_figure(basename, provenance_record, self.cfg)

    @staticmethod
    def formatter(text):
        """Text definitions to format strings."""
        repl_map = {
            "degC": "$^o$C",
            "K": "$^o$C",
            "month-1": "month$^{{-1}}$",
            "day-1": "day$^{{-1}}$",
            "d-1": "day$^{{-1}}$",
            "decade-1": "decade$^{{-1}}$",
            "year-1": "year$^{{-1}}$",
            "rcp85": "RCP8.5",
            "rcp45": "RCP4.5",
            "rcp26": "RCP2.6",
            "RCP85": "RCP8.5",
            "RCP45": "RCP4.5",
            "RCP26": "RCP2.6",
            "cmip5-85": "RCP8.5",
            "cmip5-60": "RCP6.0",
            "cmip5-45": "RCP4.5",
            "cmip5-26": "RCP2.6",
            "ssp585": "SSP5-8.5",
            "ssp245": "SSP2-4.5",
            "ssp126": "SSP1-2.6",
            "SSP585": "SSP5-8.5",
            "SSP245": "SSP2-4.5",
            "SSP126": "SSP1-2.6",
            "cmip6-85": "SSP5-8.5",
            "cmip6-70": "SSP3-7.0",
            "cmip6-60": "SSP4-6.0",
            "cmip6-34": "SSP4-3.4",
            "cmip6-45": "SSP2-4.5",
            "cmip6-26": "SSP1-2.6",
            "cmip6-19": "SSP1-1.9",
            "1": "%",
            "era5": "ERA5",
            "gpcc025x025_v8": "GPCC",
            "cru": "CRU",
            "jra55": "JRA55",
            "HIGHRESMIP": "HighResMIP",
            " ": "",
        }
        for key, val in repl_map.items():
            if key in text:
                text = text.replace(key, val)
                break
        return text

    def find_abs_bound_range(self, results_dict, keys, avg_over=5):
        """Find suitable bounds for the colorbar.

        It takes into account the absolute maximum value from all the
        panels.
        """
        max_averages = []
        min_averages = []
        for key in keys:
            result_data = results_dict[key].data
            # compress to remove masked values
            sorted_data = np.sort(result_data.compressed())
            # select the "avg_over" extreme values from the array
            # and find it's average value
            max_average_data = np.average(sorted_data[-avg_over:])
            min_average_data = np.average(sorted_data[:avg_over])
            max_averages.append(max_average_data)
            min_averages.append(min_average_data)

        # the maximum absolute value for the bound
        abs_max = np.abs(np.max(max_averages))
        abs_min = np.abs(np.min(min_averages))
        max_bound = np.max([abs_min, abs_max])

        # find the bound candidate suited for the bound range
        index = np.argwhere(self.bound_candidates - max_bound > 0)[0, 0]

        return self.bound_candidates[index]

    @staticmethod
    def region_to_square(region, dimension):
        """Definition of the region polygon."""
        if dimension == "latitude":
            boundaries = [
                region["start_latitude"],
                region["start_latitude"],
                region["end_latitude"],
                region["end_latitude"],
                region["start_latitude"],
            ]
        elif dimension == "longitude":
            boundaries = [
                region["start_longitude"],
                region["end_longitude"],
                region["end_longitude"],
                region["start_longitude"],
                region["start_longitude"],
            ]
        return boundaries

    def define_projection(self, region):
        """Projection definition to get LambertConformal borders."""
        region = {
            "start_longitude": region[0],
            "end_longitude": region[1],
            "start_latitude": region[2],
            "end_latitude": region[3],
        }
        projection = "LambertConformal"
        plotextend = [
            region["start_longitude"],
            region["end_longitude"],
            region["start_latitude"],
            region["end_latitude"],
        ]
        if projection == "LambertConformal":
            # plotextend has to be a little larger so everything is on there
            plotextend = [
                plotextend[0] - 1.0,
                plotextend[1] + 1.0,
                plotextend[2] - 1.0,
                plotextend[3] + 1.0,
            ]
            # path to cut out is exact though
            lons = self.region_to_square(region, "longitude")
            lats = self.region_to_square(region, "latitude")
            path_ext = [[lon, lat] for lon, lat in zip(lons, lats)]
            path_ext = mpath.Path(path_ext).interpolated(20)
        # South Hemisfere
        if region["start_latitude"] <= 0 and region["end_latitude"] <= 0:
            proj = ccrs.LambertConformal(
                central_longitude=np.sum(plotextend[:2]) / 2.0,
                central_latitude=np.sum(plotextend[2:]) / 2.0,
                cutoff=+30,
                standard_parallels=(-33, -45),
            )
        # North Hemisphere
        else:
            proj = ccrs.LambertConformal(
                central_longitude=np.sum(plotextend[:2]) / 2.0,
                central_latitude=np.sum(plotextend[2:]) / 2.0,
            )
        return proj, projection, path_ext, plotextend

    @staticmethod
    def sorted_dim(cube, coord="longitude"):
        """Sorts the cube data according to the longitude coordinate values.

        example: 180/-180 --> -180/180
        """
        coord_to_sort = cube.coord(coord)
        assert coord_to_sort.ndim == 1, "Coord should be 1-dimensional."
        (dim, ) = cube.coord_dims(coord_to_sort)
        index = [slice(None)] * cube.ndim
        index[dim] = np.argsort(coord_to_sort.points)
        cube = cube[tuple(index)]
        coord = cube.coord(coord)
        iris.util.promote_aux_coord_to_dim_coord(cube, "longitude")
        return cube

    def regrid_longitude_coord(self, cube):
        """Sorts the longitudes of the cubes from 0/360 degrees to -180/180."""
        # make a list with the 'longitude' coord in the form: 0/180/-180/0
        neg_lons = ((cube.coord("longitude").points + 180) % 360) - 180
        # interpolates the cube data to the new 'longitude' dimensions
        cube = cube.interpolate([("longitude", neg_lons)],
                                iris.analysis.Linear())
        sorted_cube = self.sorted_dim(cube)
        return sorted_cube

    def get_hotspot_provenance(self, suptitle, scenario, ancestor_files):
        """Create a provenance record describing the hotspot fields plots."""
        caption = (f"{suptitle}. Calculated for seasons "
                   f"{self.seasons[0].upper()}, "
                   f"{self.seasons[1].upper()} and {self.seasons[2].upper()} "
                   f"in the future periods {self.cfg['future_periods'][0]} "
                   f"and {self.cfg['future_periods'][1]} "
                   f"for CMIP5 {self.formatter(f'cmip5-{scenario}')} "
                   f"and CMIP6 {self.formatter(f'cmip6-{scenario}')}")

        record = {
            'caption': caption,
            'statistics': ['anomaly', 'diff'],
            'domains': ['reg'],
            'plot_types': ['map'],
            'authors': [
                'loosveldt-tomas_saskia',
                # 'cos_josep',
            ],
            'references': [
                'cos22esd',
            ],
            'ancestors': ancestor_files,
        }
        return record

    def get_rolling_mean_provenance(self, suptitle, ancestor_files):
        """Create a provenance record with the rolling mean diagnostic data."""
        suptitle = suptitle.replace('\n', '')
        caption = (f"{suptitle}. For CMIP5 ("
                   f"{self.formatter(f'cmip5-{self.scenarios[0]}')}, "
                   f"{self.formatter(f'cmip5-{self.scenarios[1]}')} and "
                   f"{self.formatter(f'cmip5-{self.scenarios[2]}')}) "
                   f"and CMIP6 "
                   f"({self.formatter(f'cmip6-{self.scenarios[0]}')}, "
                   f"{self.formatter(f'cmip6-{self.scenarios[1]}')} and "
                   f"{self.formatter(f'cmip6-{self.scenarios[2]}')})")

        record = {
            'caption': caption,
            'statistics': ['anomaly', "other"],
            'domains': ['reg', 'global'],
            'plot_types': ['scatter', 'line', 'times'],
            'authors': [
                'loosveldt-tomas_saskia',
                # 'cos_josep',
            ],
            'references': [
                'cos22esd',
            ],
            'ancestors': ancestor_files,
        }
        return record


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        HotspotPlot(cfg).compute()
