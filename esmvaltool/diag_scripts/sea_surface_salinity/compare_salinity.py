import logging
import os
import string
from datetime import datetime
from functools import reduce

import cf_units
import iris
import iris.quickplot as qplot
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import regrid_time, climate_statistics
from iris.coord_categorisation import add_month_number, add_year
from iris.util import unify_time_units
from matplotlib.legend import Legend
from matplotlib.legend_handler import HandlerBase
from matplotlib.text import Text

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.shared import group_metadata, names
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(__name__)


class CompareSalinity(object):
    def __init__(self, config):
        self.cfg = config

    def compute(self):
        data = group_metadata(self.cfg[names.INPUT_DATA].values(),
                              names.SHORT_NAME)
        for short_name in data:
            logger.info("Processing variable %s", short_name)
            variables = group_metadata(data[short_name], names.ALIAS)
            ref_alias = list(variables.values())[0][0]['reference_dataset']
            reference_dataset = variables.pop(ref_alias)[0]
            reference = iris.load_cube(reference_dataset[names.FILENAME])
            reference_ancestor = reference_dataset[names.FILENAME]
            logger.debug("Info reference dataset:")
            logger.debug(reference)
            for alias, dataset_info in variables.items():
                logger.info("Plotting dataset %s", alias)
                dataset_info = dataset_info[0]
                dataset = iris.load_cube(dataset_info[names.FILENAME])
                time_coord = dataset.coord('time')
                if time_coord.units.calendar == 'proleptic_gregorian':
                    time_coord.units = cf_units.Unit(
                        time_coord.units.name,
                        calendar='gregorian',
                    )
                unify_time_units((reference, dataset))
                logger.debug("Info dataset %s:", alias)
                logger.debug(dataset)
                ancestors = (dataset_info[names.FILENAME], reference_ancestor)
                for region_slice in dataset.slices_over('shape_id'):
                    region = region_slice.coord('shape_id').points[0]
                    self.create_timeseries_plot(region, region_slice,
                                                reference, ref_alias,
                                                dataset_info, ancestors)
                self.create_radar_plot(dataset_info, dataset, reference,
                                       ref_alias, ancestors)

    def create_timeseries_plot(self, region, data, reference, reference_alias,
                               dataset_info, ancestors):
        alias = dataset_info[names.ALIAS]
        qplot.plot(data, label=alias)
        qplot.plot(reference.extract(iris.Constraint(shape_id=region)),
                   label=reference_alias)
        plt.legend()
        plt.title(f"{dataset_info[names.LONG_NAME]} ({region})")
        plt.tight_layout()
        plt.savefig(f'test_timeseries_{region}.png')
        plot_path = os.path.join(
            self.cfg[names.PLOT_DIR],
            f"{dataset_info[names.SHORT_NAME]}_{region.replace(' ', '')}"
            f"_{alias}.{self.cfg[names.OUTPUT_FILE_TYPE]}")
        plt.savefig(plot_path)
        plt.close()
        caption = (f"{dataset_info[names.SHORT_NAME]} mean in {region} for "
                   f"{alias} and {reference_alias}")
        self._create_prov_record(plot_path, caption, ancestors)

    def create_radar_plot(self, data_info, data, reference, reference_alias,
                          ancestors):
        interval = self._get_overlap([data, reference])
        indices = self._slice_cube(data, interval[0], interval[1])
        data = data[indices[0]:indices[1] + 1]
        indices = self._slice_cube(reference, interval[0], interval[1])
        reference = reference[indices[0]:indices[1] + 1]

        add_month_number(data, 'time')
        add_year(data, 'time')

        add_month_number(reference, 'time')
        add_year(reference, 'time')

        data_alias = data_info[names.ALIAS]
        climat_ref = climate_statistics(reference, 'mean')
        climat_data = climate_statistics(data, 'mean')
        bias_data = np.absolute(climat_ref.data) - climat_data.data
        bias = climat_ref.copy(bias_data)
        angles = np.linspace(0, 2 * np.pi, bias.shape[0] + 1)
        # Initialise the spider plot
        ax = plt.subplot(111, polar=True)
        for spine in ax.spines.values():
            spine.set_color('grey')

        # Draw one axe per variable + add labels labels yet
        letters = [string.ascii_uppercase[i] for i in range(0, bias.shape[0])]
        plt.xticks(angles[:-1],
                   letters,
                   color='grey',
                   size=8,
                   rotation=angles[:-1])

        # Draw ylabels
        ax.set_rlabel_position(0)
        plt.yticks([3, 6, 9], ["3", "6", "9"],
                   color="grey",
                   size=7)
        plt.ylim(0, 1)

        data = np.append(bias.data, bias.data[0])
        more_angles = np.linspace(0, 2 * np.pi, bias.shape[0] * 20 + 1)
        interp_data = np.interp(more_angles, angles, data)

        # Plot data
        ax.plot(more_angles, interp_data, linewidth=1, linestyle='solid')
        ax.fill(more_angles, interp_data, 'b', alpha=0.1)
        ax.legend(letters,
                  bias.coord('shape_id').points,
                  loc='upper center',
                  ncol=2,
                  frameon=False,
                  bbox_to_anchor=(0.5, -0.1),
                  borderaxespad=0.)
        plt.title(
            f'{data_info[names.SHORT_NAME]} bias\n'
            f'{data_alias} vs {reference_alias}',
            pad=20)
        plt.tight_layout()
        plot_path = os.path.join(
            self.cfg[names.PLOT_DIR],
            f"{data_info[names.SHORT_NAME]}_comparison_{data_alias}_"
            f"{reference_alias}.{self.cfg[names.OUTPUT_FILE_TYPE]}")
        plt.savefig(plot_path)
        plt.close()
        caption = (f"Absolute bias comparison in diferent regions for "
                   f"{data_alias} and {reference_alias}")
        self._create_prov_record(plot_path, caption, ancestors)

    def _create_prov_record(self, filepath, caption, ancestors):
        record = {
            'caption': caption,
            'domains': [
                'global',
            ],
            'autors': ['vegas-regidor_javier'],
            'references': ['acknow_author'],
            'ancestors': ancestors
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(filepath, record)

    def _get_time_offset(self, time_unit):
        """Return a datetime object equivalent to tunit."""
        # tunit e.g. 'day since 1950-01-01 00:00:00.0000000 UTC'
        cfunit = cf_units.Unit(time_unit, calendar=cf_units.CALENDAR_STANDARD)
        time_offset = cfunit.num2date(0)
        return time_offset

    def _align_yearly_axes(self, cube):
        """
        Perform a time-regridding operation to align time axes for yr data.
        """
        years = [cell.point.year for cell in cube.coord('time').cells()]
        # be extra sure that the first point is not in the previous year
        if 0 not in np.diff(years):
            return regrid_time(cube, 'yr')
        return cube

    def _datetime_to_int_days(self, cube):
        """Return list of int(days) converted from cube datetime cells."""
        cube = self._align_yearly_axes(cube)
        time_cells = [cell.point for cell in cube.coord('time').cells()]

        # extract date info
        real_dates = []
        for date_obj in time_cells:
            # real_date resets the actual data point day
            # to the 1st of the month so that there are no
            # wrong overlap indices
            real_date = datetime(date_obj.year, date_obj.month, 1, 0, 0, 0)
            real_dates.append(real_date)

        # get the number of days starting from the reference unit
        time_unit = cube.coord('time').units.name
        time_offset = self._get_time_offset(time_unit)
        days = [(date_obj - time_offset).days for date_obj in real_dates]

        return days

    def _get_overlap(self, cubes):
        """
        Get discrete time overlaps.
        This method gets the bounds of coord time
        from the cube and assembles a continuous time
        axis with smallest unit 1; then it finds the
        overlaps by doing a 1-dim intersect;
        takes the floor of first date and
        ceil of last date.
        """
        all_times = []
        for cube in cubes:
            span = self._datetime_to_int_days(cube)
            start, stop = span[0], span[-1]
            all_times.append([start, stop])
        bounds = [range(b[0], b[-1] + 1) for b in all_times]
        time_pts = reduce(np.intersect1d, bounds)
        if len(time_pts) > 1:
            time_bounds_list = [time_pts[0], time_pts[-1]]
            return time_bounds_list

    def _slice_cube(self, cube, t_1, t_2):
        """
        Efficient slicer.
        Simple cube data slicer on indices
        of common time-data elements.
        """
        time_pts = [t for t in cube.coord('time').points]
        converted_t = self._datetime_to_int_days(cube)
        idxs = sorted([
            time_pts.index(ii) for ii, jj in zip(time_pts, converted_t)
            if t_1 <= jj <= t_2
        ])
        return [idxs[0], idxs[-1]]


class TextHandler(HandlerBase):
    def create_artists(self, legend, text, xdescent, ydescent, width, height,
                       fontsize, trans):
        tx = Text(width / 2.,
                  height / 2,
                  text,
                  fontsize=fontsize,
                  ha="center",
                  va="center",
                  fontweight="bold")
        return [tx]


Legend.update_default_handler_map({str: TextHandler()})


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        CompareSalinity(config).compute()


if __name__ == "__main__":
    main()
