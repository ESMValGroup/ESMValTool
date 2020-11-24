import logging
import os
import sys

import cartopy.crs as ccrs
import iris
import iris.analysis
import iris.cube
import iris.quickplot as qplt
import iris.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from dask import array as da

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvalcore.preprocessor import (annual_statistics, extract_levels, regrid,
                                     seasonal_statistics)
from esmvaltool.diag_scripts.shared import ProvenanceLogger, group_metadata

matplotlib.use('agg')

logger = logging.getLogger(os.path.basename(__file__))


class EadyGrowthRate(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.fill_value = 1e20
        self.ref_p = 1000.0
        self.g = 9.80665
        self.con = 0.3098
        self.omega = 7.292e-5
        self.time_statistic = self.cfg['time_statistic']

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            ta = iris.load_cube(var['ta'][0]['filename'])
            plev = ta.coord('air_pressure')

            theta = self.potential_temperature(ta, plev)

            del ta

            zg = iris.load_cube(var['zg'][0]['filename'])

            brunt = self.brunt_vaisala_frq(theta, zg)

            lats = zg.coord('latitude')
            fcor = self.coriolis(lats, zg.shape)

            ua = iris.load_cube(var['ua'][0]['filename'])
            if ua.shape is not zg.shape:
                ua = regrid(ua, zg, scheme='linear')

            egr = self.eady_growth_rate(fcor, ua, zg, brunt)

            cube_egr = ua.copy(egr * 86400)

            cube_egr.standard_name = None
            cube_egr.long_name = 'eady_growth_rate'
            cube_egr.var_name = 'egr'
            cube_egr.units = 'day-1'

            if self.time_statistic == 'annual_mean':
                cube_egr = annual_statistics(cube_egr)
                cube_egr = cube_egr.collapsed('time', iris.analysis.MEAN)
            elif self.time_statistic == 'seasonal_mean':
                cube_egr = seasonal_statistics(cube_egr)
                cube_egr = cube_egr.collapsed('time', iris.analysis.MEAN)
                self.seasonal_plots(cube_egr, alias)
            else:
                logger.info(
                    "Parameter time_statistic is not well set in the recipe."
                    "Must be 'annual_mean' or 'seasonal_mean'")
                sys.exit()

            self.save(cube_egr, alias, data)

    def potential_temperature(self, ta, plev):
        p0 = iris.coords.AuxCoord(self.ref_p,
                                  long_name='reference_pressure',
                                  units='hPa')
        p0.convert_units(plev.units)
        p = (p0.points / plev.points)**(2 / 7)
        theta = ta * iris.util.broadcast_to_shape(
            p, ta.shape, ta.coord_dims('air_pressure'))
        theta.long_name = 'potential_air_temperature'

        return theta

    def vertical_integration(self, x, y):
        # Perform a non-cyclic centered finite-difference to integrate
        # along pressure levels

        plevs = x.shape[1]

        dxdy_0 = ((x[:, 1, :, :].lazy_data() - x[:, 0, :, :].lazy_data()) /
                  (y[:, 1, :, :].lazy_data() - y[:, 0, :, :].lazy_data()))

        dxdy_centre = ((x[:, 2:plevs, :, :].lazy_data() -
                        x[:, 0:plevs - 2, :, :].lazy_data()) /
                       (y[:, 2:plevs, :, :].lazy_data() -
                        y[:, 0:plevs - 2, :, :].lazy_data()))

        dxdy_end = ((x[:, plevs - 1, :, :].lazy_data() -
                     x[:, plevs - 2, :, :].lazy_data()) /
                    (y[:, plevs - 1, :, :].lazy_data() -
                     y[:, plevs - 2, :, :].lazy_data()))

        bounds = [dxdy_end, dxdy_0]
        stacked_bounds = da.stack(bounds, axis=1)
        total = [dxdy_centre, stacked_bounds]

        # Concatenate arrays where the last slice is dxdy_0
        dxdy = da.concatenate(total, axis=1)

        # Move dxdy_0 to the beggining of the array
        dxdy = da.roll(dxdy, 1, axis=1)

        return dxdy

    def brunt_vaisala_frq(self, theta, zg):
        dthdz = self.vertical_integration(theta, zg)
        dthdz = da.where(dthdz > 0, dthdz, 0)
        buoy = (self.g / theta.lazy_data()) * dthdz
        brunt = da.sqrt(buoy)
        brunt = da.where(brunt != 0, brunt, self.fill_value)

        return brunt

    def coriolis(self, lats, ndim):
        fcor = 2.0 * self.omega * np.sin(np.radians(lats.points))
        fcor = fcor[np.newaxis, np.newaxis, :, np.newaxis]
        fcor = da.broadcast_to(fcor, ndim)

        return fcor

    def eady_growth_rate(self, fcor, ua, zg, brunt):
        dudz = self.vertical_integration(ua, zg)
        egr = self.con * abs(fcor) * abs(dudz) / brunt

        return egr

    def seasonal_plots(self, egr, alias):
        try:
            levels = self.cfg['plot_levels']
        except KeyError:
            logger.info("Parameter plot_levels is not set in the recipe."
                        "Plotting all pressure levels instead.")
            levels = egr.coord('air_pressure').points
        for level in levels:
            cube = extract_levels(egr, level, scheme='linear')
            crs_latlon = ccrs.PlateCarree()
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines(linewidth=1, color='black')
            # North Atlantic
            ax.set_extent((-90.0, 30.0, 20.0, 80.0), crs=crs_latlon)
            ax.set_yticks(np.linspace(25, 75, 6))
            qplt.contourf(cube, levels=np.arange(0, np.max(cube.data), 0.05))
            extension = self.cfg['output_file_type']
            diagnostic = self.cfg['script']
            plotname = '_'.join([alias, diagnostic,
                                 str(int(level))]) + f'.{extension}'
            plt.savefig(os.path.join(self.cfg[n.PLOT_DIR], plotname))
            plt.close()

    def save(self, egr, alias, data):
        script = self.cfg[n.SCRIPT]
        info = data[alias][0]
        keys = [
            str(info[key]) for key in ('project', 'dataset', 'exp', 'ensemble',
                                       'diagnostic', 'start_year', 'end_year')
            if key in info
        ]
        output_name = '_'.join(keys) + '.nc'
        output_file = os.path.join(self.cfg[n.WORK_DIR], output_name)
        iris.save(egr, output_file)

        caption = ("{script} between {start} and {end}"
                   "according to {dataset}").format(script=script.replace(
                       " ", '_'),
                                                    start=info['start_year'],
                                                    end=info['end_year'],
                                                    dataset=info['dataset'])
        ancestors = []
        for i in range(len(data[alias])):
            ancestors.append(data[alias][i]['filename'])
        record = {
            'caption': caption,
            'domains': ['global'],
            'authors': ['sanchez-gomez_emilia', 'moreno-chamarro_eduardo'],
            'references': ['acknow_project'],
            'ancestors': ancestors
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(output_file, record)


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        EadyGrowthRate(config).compute()


if __name__ == "__main__":
    main()
