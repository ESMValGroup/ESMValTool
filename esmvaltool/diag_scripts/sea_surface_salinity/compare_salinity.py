import logging
import string
import os
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase
from matplotlib.text import Text
from matplotlib.legend import Legend
import numpy as np
import cf_units
import iris
from iris.util import unify_time_units
import iris.quickplot as qplot

from iris.analysis.stats import pearsonr
from iris.coord_categorisation import add_month_number, add_year

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata
from esmvalcore.preprocessor._multimodel import _get_overlap, _slice_cube
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger


logger = logging.getLogger(__name__)


class CompareSalinity(object):
    def __init__(self, config):
        self.cfg = config

    def compute(self):
        data = group_metadata(self.cfg[n.INPUT_DATA].values(), n.SHORT_NAME)
        for short_name in data:
            logger.info("Processing variable %s", short_name)
            variables = group_metadata(data[short_name], n.ALIAS)
            ref_alias = list(variables.values())[0][0]['reference_dataset']
            reference_dataset = variables.pop(ref_alias)[0]
            reference = iris.load_cube(reference_dataset[n.FILENAME])
            reference_ancestor = reference_dataset[n.FILENAME]
            logger.debug(reference)
            for alias, dataset_info in variables.items():
                logger.info("Plotting dataset %s", alias)
                dataset_info = dataset_info[0]
                dataset = iris.load_cube(dataset_info[n.FILENAME])
                time_coord = dataset.coord('time')
                if time_coord.units.calendar == 'proleptic_gregorian':
                    time_coord.units = cf_units.Unit(
                        time_coord.units.name,
                        calendar='gregorian',
                    )
                unify_time_units((reference, dataset))
                logger.debug(dataset)
                ancestors = (dataset_info[n.FILENAME], reference_ancestor)
                for region_slice in dataset.slices_over('shape_id'):
                    region = region_slice.coord('shape_id').points[0]
                    self.create_timeseries_plot(
                        region, region_slice, reference, ref_alias,
                        dataset_info, ancestors
                    )
                self.create_radar_plot(
                    dataset_info, dataset, reference, ref_alias,
                    ancestors)

    def create_timeseries_plot(self, region, data, reference, reference_alias,
                               dataset_info, ancestors):
        alias = dataset_info[n.ALIAS]
        qplot.plot(data, label=alias)
        qplot.plot(
            reference.extract(iris.Constraint(shape_id=region)),
            label=reference_alias
        )
        plt.legend()
        plt.title(f"{dataset_info[n.LONG_NAME]} ({region})")
        plt.tight_layout()
        plt.savefig(f'test_timeseries_{region}.png')
        plot_path = os.path.join(
            self.cfg[n.PLOT_DIR],
            f"{dataset_info[n.SHORT_NAME]}_{region.replace(' ', '')}_{alias}"
            f".{self.cfg[n.OUTPUT_FILE_TYPE]}"
        )
        plt.savefig(plot_path)
        plt.close()
        caption = (
            f"{dataset_info[n.SHORT_NAME]} mean in {region} for "
            f"{alias} and {reference_alias}"
        )
        self._create_prov_record(plot_path, caption, ancestors)

    def create_radar_plot(self, data_info, data, reference, reference_alias,
                          ancestors):
        interval = _get_overlap([data, reference])
        indices = _slice_cube(data, interval[0], interval[1])
        data = data[indices[0]:indices[1] + 1]
        indices = _slice_cube(reference, interval[0], interval[1])
        reference = reference[indices[0]:indices[1] + 1]

        add_month_number(data, 'time')
        add_year(data, 'time')
        data.remove_coord('time')

        add_month_number(reference, 'time')
        add_year(reference, 'time')
        reference.remove_coord('time')

        data_alias = data_info[n.ALIAS]
        corr = pearsonr(data, reference, ('month_number', 'year'))
        angles = np.linspace(0, 2 * np.pi, corr.shape[0] + 1)
        # Initialise the spider plot
        ax = plt.subplot(111, polar=True)
        for spine in ax.spines.values():
            spine.set_color('grey')

        # Draw one axe per variable + add labels labels yet
        letters = [string.ascii_uppercase[i] for i in range(0, corr.shape[0])]
        plt.xticks(
            angles[:-1], letters, color='grey', size=8, rotation=angles[:-1])

        # Draw ylabels
        ax.set_rlabel_position(0)
        plt.yticks(
            [0.25, 0.5, 0.75], ["0.25", "0.5", "0.75"], color="grey", size=7)
        plt.ylim(0, 1)

        data = np.append(corr.data, corr.data[0])

        more_angles = np.linspace(0, 2 * np.pi, corr.shape[0] * 20 + 1)
        interp_data = np.interp(more_angles, angles, data)

        # Plot data
        ax.plot(more_angles, interp_data, linewidth=1, linestyle='solid')
        ax.fill(more_angles, interp_data, 'b', alpha=0.1)
        ax.legend(
            letters, corr.coord('shape_id').points,
            loc='upper center', ncol=2, frameon=False,
            bbox_to_anchor=(0.5, -0.1), borderaxespad=0.
        )
        plt.title(
            f'{data_info[n.SHORT_NAME]} correlation\n'
            f'{data_alias} vs {reference_alias}',
            pad=20
        )
        plt.tight_layout()
        plot_path = os.path.join(
            self.cfg[n.PLOT_DIR],
            f"{data_info[n.SHORT_NAME]}_comparison_{data_alias}_"
            f"{reference_alias}.{self.cfg[n.OUTPUT_FILE_TYPE]}"
        )
        plt.savefig(plot_path)
        plt.close()
        caption = (
            f"Correlation comparison in diferent regions for "
            f"{data_alias} and {reference_alias}"
        )
        self._create_prov_record(plot_path, caption, ancestors)

    def _create_prov_record(self, filepath, caption, ancestors):
        record = {
            'caption': caption,
            'domains': ['global', ],
            'autors': ['vegas-regidor_javier'],
            'references': ['contact_author'],
            'ancestors': ancestors
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(filepath, record)


class TextHandler(HandlerBase):
    def create_artists(self, legend, text, xdescent, ydescent,
                       width, height, fontsize, trans):
        tx = Text(width/2., height/2, text, fontsize=fontsize,
                  ha="center", va="center", fontweight="bold")
        return [tx]


Legend.update_default_handler_map({str: TextHandler()})


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        CompareSalinity(config).compute()


if __name__ == "__main__":
    main()
