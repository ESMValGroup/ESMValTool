import os
import logging
import string

import iris
import iris.cube
import iris.analysis
import iris.util

import numpy as np

from dask import array as da

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata


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

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            ta = iris.load_cube(var['ta'][0]['filename'])
            plev = ta.dim_coords[1]

            theta = self.potential_temperature(ta, plev)

            del ta

            zg = iris.load_cube(var['zg'][0]['filename'])

            brunt = self.brunt_vaisala_frq(theta, zg)

            lats = zg.dim_coords[2]
            fcor = self.coriolis(lats, zg.shape)

            ua = iris.load_cube(var['ua'][0]['filename'])

            egr = self.eady_growth_rate(fcor, ua, zg, brunt)

            cube_egr = ua.copy(egr * 86400)

            cube_egr = cube_egr.collapsed('time', iris.analysis.MEAN)

            self.save(cube_egr, alias, data)

    def potential_temperature(self, ta, plev):
        p0 = iris.coords.AuxCoord(self.ref_p, long_name='reference_pressure',
                                  units='hPa')
        p0.convert_units(plev.units)
        p = (p0.points/plev.points)**(2/7)
        theta = ta * iris.util.broadcast_to_shape(p, ta.shape, (1,))
        theta.long_name = 'potential_air_temperature'
        return theta

    def vertical_integration(self, x, y):
        # Perform a non-cyclic centered finite-difference to integrate
        # along pressure levels

        plevs = x.shape[1]

        dxdy_0 = ((x[:, 1, :, :].lazy_data() - x[:, 0, :, :].lazy_data()) /
                  (y[:, 1, :, :].lazy_data() - y[:, 0, :, :].lazy_data()))

        dxdy_centre = ((x[:, 2:plevs, :, :].lazy_data() -
                        x[:, 0:plevs-2, :, :].lazy_data()) /
                       (y[:, 2:plevs, :, :].lazy_data() -
                        y[:, 0:plevs-2, :, :].lazy_data()))

        dxdy_end = ((x[:, plevs-1, :, :].lazy_data() -
                     x[:, plevs-2, :, :].lazy_data()) /
                    (y[:, plevs-1, :, :].lazy_data() -
                     y[:, plevs-2, :, :].lazy_data()))

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
        fcor = 2.0*self.omega*np.sin(np.radians(lats.points))
        fcor = fcor[np.newaxis, np.newaxis, :, np.newaxis]
        fcor = da.broadcast_to(fcor, ndim)

        return fcor

    def eady_growth_rate(self, fcor, ua, zg, brunt):
        dudz = self.vertical_integration(ua, zg)
        egr = self.con * abs(fcor) * abs(dudz) / brunt

        return egr

    def save(self, egr, alias, data):
        egr.standard_name = None
        egr.long_name = 'eady_growth_rate'
        egr.var_name = 'egr'
        egr.units = 'day-1'

        script = self.cfg[n.SCRIPT]
        dataset = data[alias][0]['dataset']
        filename = '{dataset}_{script}.nc'.format(
            dataset=dataset,
            script=script
        )

        iris.save(egr, os.path.join(self.cfg[n.WORK_DIR], filename))


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        EadyGrowthRate(config).compute()


if __name__ == "__main__":
    main()
