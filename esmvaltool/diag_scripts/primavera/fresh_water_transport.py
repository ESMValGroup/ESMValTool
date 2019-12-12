import os
import logging
import collections
import calendar
import warnings

import numpy as np
import matplotlib.pyplot as plt

import iris
import iris.analysis
import iris.plot as iplot

from geopy.distance import great_circle

from esmvalcore.preprocessor import extract_region, mask_above_threshold, \
    mask_outside_range

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(os.path.basename(__file__))


class FreshWaterTransport(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.output_name = None
        self.target_grid = self.cfg.get('target_grid')
        self.grid_cube = None
        self.sftlf = None

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        straits = self.cfg.get(
            'straits',
            {
                'BERINGSTRAIT': {'latitude': (63.9, 67.3), 'longitude': (-174, -165)},
                'FRAMSTRAIT': {'latitude': (77, 78), 'longitude': (-20, 14)},
                'GIN': {'latitude': (65, 82), 'longitude': (-20, 14)},
                'LABRADOR': {'latitude': (45, 72), 'longitude': (-70, -40)},
                'BARENTS': {'latitude': (79.7, 80.3), 'longitude': (17.5, 52)},
                'GREEN-ICE': {'latitude': (64.6, 66), 'longitude': (-41, -21)},
                'ICE-NORWAY': {'latitude': (64.6, 65), 'longitude': (-15.7, 12.4)},
                'NORTHBAFFIN': {'latitude': (78, 79), 'longitude': (-78.4, -63)},
            }
        )
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            logger.info('Processing %s', alias)
            for strait, config in straits.items():
                logger.info('Computing %s', strait)
                so = iris.load_cube(var['so'][0]['filename'])
                vo = iris.load_cube(var['vo'][0]['filename'])
                self._plot_salinity_profile(alias, so, vo, strait, config)

    def _plot_salinity_profile(self, alias, so, vo, strait, config):
        so = extract_region(so, config['longitude'][0], config['longitude'][1],
                            config['latitude'][0], config['latitude'][1])
        vo = extract_region(vo, config['longitude'][0], config['longitude'][1],
                            config['latitude'][0], config['latitude'][1])

        latitude_index, longitude_index = self._get_lat_lon_index(so)
        so = so[:, :, latitude_index, longitude_index]
        vo = vo[:, :, latitude_index, longitude_index]
        lat_coord = so.coord('latitude')
        lon_coord = so.coord('longitude')

        distance = np.zeros(len(longitude_index))
        for x in range(len(longitude_index)):
            point0 = (lat_coord.bounds[x][0], lon_coord.bounds[x][0])
            point1 = (lat_coord.bounds[x][1], lon_coord.bounds[x][1])
            point2 = (lat_coord.bounds[x][2], lon_coord.bounds[x][2])
            point3 = (lat_coord.bounds[x][3], lon_coord.bounds[x][3])
            distance[x] = (
                great_circle(point0, point1).m +
                great_circle(point2, point3).m
            ) / 2.

        so = mask_outside_range(so, 10, 45)
        vo = mask_above_threshold(vo, 30)

        transport = so.copy(so.lazy_data() * vo.lazy_data())
        depth = transport.coord('depth')
        gradient = np.diff(depth.points, prepend=0)
        area = np.outer(gradient, distance)
        salt_ts = self._get_time_series(transport, area)
        velocity_ts = self._get_time_series(vo, area)
        fwt = (velocity_ts - (1/34.8) * salt_ts.lazy_data()) / 1000.
        fwt.standard_name = None
        fwt.var_name = 'fwt'
        fwt.long_name = 'Fresh water transport'
        fwt.units = 'mSv'
        iris.save(
            fwt,
            os.path.join(self.cfg[n.WORK_DIR], f'{alias}_{strait}_fwt.nc')
        )
        self._plot(
            fwt.aggregated_by('year', iris.analysis.SUM),
            f'{alias} (Yearly)',
            f'{alias}_{strait}_yearly_fwt'
        )
        monthly_fwt = fwt.aggregated_by(
            ('month_number', 'year'), iris.analysis.SUM)

        for month in range(1, 13):
            cube = monthly_fwt.extract(iris.Constraint(month_number=month))
            title = f'{alias} ({calendar.month_name[month]})'
            filename = f'{alias}_{strait}_{calendar.month_abbr[month]}_fwt'
            self._plot(cube, title, filename)

    def _get_lat_lon_index(self, so):
        data_mask = so[0, 0, ...].data.mask
        ind_lat, _ = np.where(data_mask == False)
        counter = collections.Counter(ind_lat)
        most_com = counter.most_common(1)
        ss = most_com[0]
        latitude_index = counter.most_common(1)[0][0]
        data_mask = data_mask[latitude_index, ...]
        longitude_index = np.where(data_mask == False)[0]
        return latitude_index, longitude_index

    def _plot(self, cube, title, filename):
        if not self.cfg[n.WRITE_PLOTS]:
            return
        iplot.plot(cube.coord('year'), cube)
        plt.title(title)
        plt.ylabel(f'{cube.long_name} ({cube.units})')
        plt.xlabel('Year')
        plt.grid(True)

        plt.savefig(os.path.join(
            self.cfg[n.PLOT_DIR],
            f'{filename}.{self.cfg[n.OUTPUT_FILE_TYPE]}'
        ))
        plt.close()

    @staticmethod
    def _get_time_series(cube, area):
        cube = cube.copy(cube.lazy_data() * area)
        cube.coord('latitude').bounds = None
        cube.coord('longitude').bounds = None
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Collapsing spatial coordinate*")
            warnings.filterwarnings(
                "ignore", message="Collapsing a non-contiguous coordinate.*")
            cube = cube.collapsed(
                ('latitude', 'longitude', 'depth'),
                iris.analysis.SUM
            )
        return cube


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        FreshWaterTransport(config).compute()


if __name__ == "__main__":
    main()
