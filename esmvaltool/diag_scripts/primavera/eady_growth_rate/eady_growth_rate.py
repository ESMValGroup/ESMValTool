"""Diagnostic for PRIMAVERA Eady Growth Rate."""
import logging
import os
import sys

import cartopy.crs as ccrs
import iris
import iris.analysis
import iris.cube
import iris.quickplot as qplt
import iris.util
import matplotlib.pyplot as plt
import numpy as np
from dask import array as da
from esmvalcore.preprocessor import (
    annual_statistics,
    extract_levels,
    regrid,
    seasonal_statistics,
)

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    names,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


class EadyGrowthRate:
    """Class used to compute the Eady Growth Rate."""

    def __init__(self, config):
        """
        Set diagnostic parameters and constants.

        Parameters
        ----------
            config : dict
                Dictionary containing configuration settings.
        """
        self.cfg = config
        self.fill_value = 1e20
        """Fill Value."""
        self.ref_p = 1000.0
        """Reference Pressure [Pa]."""
        self.gravity = 9.80665
        """Gravity [m/s]."""
        self.con = 0.3098
        """Constant."""
        self.omega = 7.292e-5
        """Rotation of the Earth [rad/s]."""
        self.time_statistic = self.cfg['time_statistic']
        """Time statistic to perform."""

    def compute(self):
        """Compute Eady Growth Rate and either it's annual or seasonal mean."""
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            temperature = iris.load_cube(var['ta'][0]['filename'])
            plev = temperature.coord('air_pressure')

            theta = self.potential_temperature(temperature, plev)

            del temperature

            geopotential = iris.load_cube(var['zg'][0]['filename'])

            brunt = self.brunt_vaisala_frq(theta, geopotential)

            lats = geopotential.coord('latitude')
            fcor = self.coriolis(lats, geopotential.shape)

            eastward_wind = iris.load_cube(var['ua'][0]['filename'])
            if eastward_wind.shape is not geopotential.shape:
                eastward_wind = regrid(eastward_wind,
                                       geopotential,
                                       scheme='linear')

            egr = self.eady_growth_rate(fcor, eastward_wind, geopotential,
                                        brunt)

            cube_egr = eastward_wind.copy(egr * 86400)

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

    def potential_temperature(self, temperature, plev):
        """Compute potential temperature.

        Parameters
        ----------
        temperature: iris.cube.Cube
            Cube of air temperature ta.
        plev: iris.coords.Coord
            Pressure level coordinates

        Returns
        -------
        theta: iris.cube.Cube
            Cube of potential temperature theta.
        """
        reference_pressure = iris.coords.AuxCoord(
            self.ref_p, long_name='reference_pressure', units='hPa')
        reference_pressure.convert_units(plev.units)
        pressure = (reference_pressure.points / plev.points)**(2 / 7)
        theta = temperature * iris.util.broadcast_to_shape(
            pressure, temperature.shape,
            temperature.coord_dims('air_pressure'))
        theta.long_name = 'potential_air_temperature'

        return theta

    @staticmethod
    def vertical_integration(var_x, var_y):
        """
        Vertical integration.

        Perform a non-cyclic centered finite-difference to integrate
        variable x with respect to variable y along pressure levels.

        Parameters
        ----------
        x: iris.cube.Cube
            Cube of variable x.
        y: iris.cube.Cube
            Cube of variable y.

        Returns
        -------
        dxdy: iris.cube.Cube
            Cube of variable integrated along pressure levels.
        """
        plevs = var_x.shape[1]

        dxdy_0 = (
            (var_x[:, 1, :, :].lazy_data() - var_x[:, 0, :, :].lazy_data()) /
            (var_y[:, 1, :, :].lazy_data() - var_y[:, 0, :, :].lazy_data()))

        dxdy_centre = ((var_x[:, 2:plevs, :, :].lazy_data() -
                        var_x[:, 0:plevs - 2, :, :].lazy_data()) /
                       (var_y[:, 2:plevs, :, :].lazy_data() -
                        var_y[:, 0:plevs - 2, :, :].lazy_data()))

        dxdy_end = ((var_x[:, plevs - 1, :, :].lazy_data() -
                     var_x[:, plevs - 2, :, :].lazy_data()) /
                    (var_y[:, plevs - 1, :, :].lazy_data() -
                     var_y[:, plevs - 2, :, :].lazy_data()))

        bounds = [dxdy_end, dxdy_0]
        stacked_bounds = da.stack(bounds, axis=1)
        total = [dxdy_centre, stacked_bounds]

        # Concatenate arrays where the last slice is dxdy_0
        dxdy = da.concatenate(total, axis=1)

        # Move dxdy_0 to the beggining of the array
        dxdy = da.roll(dxdy, 1, axis=1)

        return dxdy

    def brunt_vaisala_frq(self, theta, geopotential):
        """Compute Brunt-Väisälä frequency.

        Parameters
        ----------
        theta: iris.cube.Cube
            Cube of potential temperature.
        geopotential: iris.cube.Cube
            Cube of variable zg.

        Returns
        -------
        brunt: da.array
            Array containing Brunt-Väisälä frequency.
        """
        dthdz = self.vertical_integration(theta, geopotential)
        dthdz = da.where(dthdz > 0, dthdz, 0)
        buoy = (self.gravity / theta.lazy_data()) * dthdz
        brunt = da.sqrt(buoy)
        brunt = da.where(brunt != 0, brunt, self.fill_value)

        return brunt

    def coriolis(self, lats, ndim):
        """Compute Coriolis force.

        Parameters
        ----------
        lats: iris.coord.Coord
            Latitude coordinate.
        ndim: int
            Number of dimension.

        Returns
        -------
        fcor: da.array
            Array containing Coriolis force.
        """
        fcor = 2.0 * self.omega * np.sin(np.radians(lats.points))
        fcor = fcor[np.newaxis, np.newaxis, :, np.newaxis]
        fcor = da.broadcast_to(fcor, ndim)

        return fcor

    def eady_growth_rate(self, fcor, eastward_wind, geopotential, brunt):
        """Compute Eady Growth Rate.

        Parameters
        ----------
        fcor: da.array
            Array containing Coriolis force.
        eastward_wind: iris.cube.Cube
            Cube containing variable ua.
        geopotential: iris.cube.Cube
            Cube containing variable zg.
        brunt: da.array
            Array containing Brunt-Väisäla frequency

        Returns
        -------
        egr: da.array
            Array containing Eady Growth Rate.
        """
        dudz = self.vertical_integration(eastward_wind, geopotential)
        egr = self.con * abs(fcor) * abs(dudz) / brunt

        return egr

    def seasonal_plots(self, egr, alias):
        """
        Plot seasonal Eady Growth rate values.

        Parameters
        ----------
        egr: iris.cube.Cube
            Cube containing variable egr.
        alias: str
            Alias of the dataset.
        """
        try:
            levels = self.cfg['plot_levels']
        except KeyError:
            logger.info("Parameter plot_levels is not set in the recipe."
                        "Plotting all pressure levels instead.")
            levels = egr.coord('air_pressure').points
        for level in levels:
            cube = extract_levels(egr, level, scheme='linear')
            crs_latlon = ccrs.PlateCarree()
            axes = plt.axes(projection=ccrs.PlateCarree())
            axes.coastlines(linewidth=1, color='black')
            # North Atlantic
            axes.set_extent((-90.0, 30.0, 20.0, 80.0), crs=crs_latlon)
            axes.set_yticks(np.linspace(25, 75, 6))
            # Relevant range
            qplt.contourf(cube, levels=np.arange(0, 1.1, 0.05))
            extension = self.cfg['output_file_type']
            diagnostic = self.cfg['script']
            plotname = '_'.join([alias, diagnostic,
                                 str(int(level))]) + f'.{extension}'
            plt.savefig(os.path.join(self.cfg[names.PLOT_DIR], plotname))
            plt.close()

    def save(self, egr, alias, data):
        """Save results and write provenance."""
        script = self.cfg[names.SCRIPT]
        info = data[alias][0]
        keys = [
            str(info[key]) for key in ('project', 'dataset', 'exp', 'ensemble',
                                       'diagnostic', 'start_year', 'end_year')
            if key in info
        ]
        output_name = '_'.join(keys) + '.nc'
        output_file = os.path.join(self.cfg[names.WORK_DIR], output_name)
        iris.save(egr, output_file)

        script_name = script.replace(" ", '_')
        caption = (f"{script_name} between {info['start_year']} "
                   f"and {info['end_year']} according to {info['dataset']}")
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
    """Run Eady Growth Rate diagnostic."""
    with run_diagnostic() as config:
        EadyGrowthRate(config).compute()


if __name__ == "__main__":
    main()
