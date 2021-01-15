"""Diagnostic for PRIMAVERA Cyclone Tracker."""
import os
import logging
import shutil

import iris
import iris.analysis
import iris.cube
import iris.util

import cf_units

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    names,
    ProvenanceLogger,
    run_diagnostic
)

logger = logging.getLogger(os.path.basename(__file__))


class CycloneTracker:
    """Class used to prepare the input for the Cyclone Tracker."""

    def __init__(self, config):
        """Set diagnostic parameters and constants.

        Parameters
        ----------
            config : dict
                Dictionary containing configuration settings.
        """
        self.cfg = config
        self.atcffreq = None
        """Frequency of the output [x100 hr]."""
        self.westbd = self.cfg['westbd']
        """West boundary of the domain."""
        self.eastbd = self.cfg['eastbd']
        """East boundary of the domain."""
        self.northbd = self.cfg['northbd']
        """North boundary of the domain."""
        self.southbd = self.cfg['southbd']
        """South boundary of the domain."""
        self.tracktype = self.cfg['type']
        """Type of tracking to perform."""
        self.mslpthresh = self.cfg['mslpthresh']
        """Minimum MSLP gradient threshold [hPa/km]."""
        self.mslpthresh2 = self.cfg['mslpthresh2']
        """Maximum pressure threshold [hPa]."""
        self.v850thresh = self.cfg['v850thresh']
        """Minimum azimutally-average 850hPa windspeed threshold [m/s]."""
        self.contint = self.cfg['contint']
        """Interval minimum pressure of a new low and the closed isobar[Pa]."""
        self.wcore_depth = self.cfg['wcore_depth']
        """Contour interval [K]."""
        self.ikeflag = self.cfg['ikeflag']
        """Compute Integrated Kinetic Energy (y/n)."""
        self.verb = self.cfg['verb']
        """Verbosity (0/1/2/3)."""
        self.tracker_exe = self.cfg['tracker_exe']
        """Path to the tracker executable."""

    def compute(self):
        """
        Prepares the input and calls the tracker.
        """
        data = group_metadata(self.cfg['input_data'].values(), 'dataset')
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            psl = iris.load_cube(var['psl'][0]['filename'])
            uas = iris.load_cube(var['uas'][0]['filename'])
            vas = iris.load_cube(var['vas'][0]['filename'])
            if 'ua7h' not in var.keys():
                eastward_wind = iris.load_cube(var['ua'][0]['filename'])
                northward_wind = iris.load_cube(var['va'][0]['filename'])
                temperature = iris.load_cube(var['ta'][0]['filename'])
                try:
                    geopotential = iris.load_cube(var['zg'][0]['filename'])
                except KeyError:
                    geopotential = iris.load_cube(var['zg7h'][0]['filename'])
                    geopotential.var_name = 'zg'

            else:
                eastward_wind = iris.load_cube(var['ua7h'][0]['filename'])
                eastward_wind.var_name = 'ua'
                northward_wind = iris.load_cube(var['va7h'][0]['filename'])
                northward_wind.var_name = 'va'
                temperature = iris.load_cube(var['ta7h'][0]['filename'])
                temperature.var_name = 'ta'
                geopotential = iris.load_cube(var['zg7h'][0]['filename'])
                geopotential.var_name = 'zg'

                eastward_wind.remove_coord('time')
                eastward_wind.add_dim_coord(psl.coord('time'), 0)
                northward_wind.remove_coord('time')
                northward_wind.add_dim_coord(psl.coord('time'), 0)
                uas = iris.load_cube(var['uas'][0]['filename'])
                uas.remove_coord('time')
                uas.add_dim_coord(psl.coord('time'), 0)
                vas = iris.load_cube(var['vas'][0]['filename'])
                vas.remove_coord('time')
                vas.add_dim_coord(psl.coord('time'), 0)
                temperature.remove_coord('time')
                temperature.add_dim_coord(psl.coord('time'), 0)
                geopotential.remove_coord('time')
                geopotential.add_dim_coord(psl.coord('time'), 0)

            temperature = temperature.collapsed('air_pressure',
                                                iris.analysis.MEAN)
            temperature.coords('air_pressure')[0].points = [40100]
            self.atcffreq = '600'
            freq = psl.attributes['frequency']
            if 'day' in freq:
                self.atcffreq = '2400'

            try:
                years = set(psl.coord('year').points)
            except iris.exceptions.CoordinateNotFoundError:
                iris.coord_categorisation.add_year(psl, 'time')
                years = set(psl.coord('year').points)
            try:
                months = set(psl.coord('month_number').points)
            except iris.exceptions.CoordinateNotFoundError:
                iris.coord_categorisation.add_month_number(psl, 'time')
                months = set(psl.coord('month_number').points)

            total = [psl, eastward_wind, northward_wind,
                     uas, vas, temperature, geopotential]
            output = '{project}_' \
                     '{dataset}_' \
                     '{mip}_' \
                     '{exp}_' \
                     '{ensemble}_' \
                     'cyclone_information_' \
                     '{start}_' \
                     '{end}.txt'.format(project=data[alias][0]['project'],
                                        dataset=data[alias][0]['dataset'],
                                        exp=data[alias][0]['exp'],
                                        mip=data[alias][0]['mip'],
                                        ensemble=data[alias][0]['ensemble'],
                                        start=data[alias][0]['start_year'],
                                        end=data[alias][0]['end_year'])

            output_file = open(os.path.join(self.cfg[names.WORK_DIR], output),
                               'wb')

            try:
                start_day = psl.coord('day_of_month').points[0]
            except iris.exceptions.CoordinateNotFoundError:
                iris.coord_categorisation.add_day_of_month(psl, 'time')
                start_day = psl.coord('day_of_month').points[0]
            end_day = psl.coord('day_of_month').points[-1]
            self.run_custom_time(data[alias][0]['dataset'],
                                 total,
                                 years,
                                 months,
                                 start_day,
                                 end_day,
                                 output_file)
            output_file.close()
#            self.write_provenance(
#                alias, data, os.path.join(self.cfg[n.WORK_DIR], output)
#                )

    def run_custom_time(self, dataset, total, years,
                        months, start_day, end_day, output_file):
        """
        Prepares output file for a custom time period.

        Parameters
        ----------
        dataset: str
            Name of the dataset
        total: list
            List containing all the variables' cubes.
        years: set
            Value of the years present in the variables.
        months: set
            Value of the months present in the variables.
        start_day: int
            Value of the variables' start day.
        end_day: int
            Value of the variables' end day.
        output_file: file
            Output txt file.

        """
        total_period = []
        for i, _ in enumerate(total):
            total_period.append(total[i])
            calendar = total_period[i].coord('time').units.calendar
            total_period[i].coord('time').convert_units(cf_units.Unit(
                'hours since {0}-{1}-{2} 00:00:00'.format(
                    min(years), min(months), start_day), calendar=calendar))

        filename = '{dataset}_' \
                   '{start_year}{start_month}{start_day}_' \
                   '{end_year}{end_month}{end_day}'.format(
                       dataset=dataset,
                       start_year=list(years)[0],
                       start_month=list(months)[0],
                       start_day=start_day,
                       end_year=list(years)[-1],
                       end_month=list(months)[-1],
                       end_day=end_day)
        self.call_tracker(total, filename, output_file,
                          list(years)[0], list(months)[0])

    def call_tracker(self, variables, filename, output, year, month):
        """
        Calls the tracker.

        Parameters
        ----------
        variables: list.
            List containing the variables' cube.
        filename: str.
            NetCDF file containing the input for the tracker.
        output: str
            Output text file.
        year: int.
            First year present in the dataset.
        month: int.
            First month present in the dataset.

        """
        input_path = os.path.join(self.cfg[names.WORK_DIR], filename + '.nc')
        iris.save(variables, input_path)
        time = variables[0].coord('time').points
        path = os.path.join(self.cfg['run_dir'], filename)
        if not os.path.isdir(path):
            os.makedirs(path)
        os.system('ln -s {0} {1}/fort.11'.format(input_path, path))
        self.write_namelist(path, month, year)
        self.write_fort15(path, time)
        self.write_fort14(path)
        os.chdir(path)
        args = self.tracker_exe + ' < namelist'
        os.system(args)
        shutil.copyfileobj(open('fort.66', 'rb'), output)

    def write_fort15(self, path, time):
        """
        Writes the fort.15 file.

        Parameters
        ----------
        path: str
            Path to the fort.15 file
        time: list
            List containing the variables' timesteps.
        """
        fort15_file = open(os.path.join(path, 'fort.15'), 'w')
        for timestep, value in enumerate(time):
            fort15_file.write('{0:4} {1:5}\n'.format(timestep + 1,
                                                     int(value * 60)))
        fort15_file.close()

    def write_fort14(self, path):
        """
        Writes the fort.14 file.

        Parameters
        ----------
        path: str
            Path to the fort.14 file
        """
        fort14_file = open(os.path.join(path, 'fort.14'), 'w')
        fort14_file.close()

    def write_namelist(self, path, month, year):
        """
        Write namelist

        Parameters
        ----------
        path: str
            Path to the namelist.
        month: int
            Input month for the run.
        year: int
            Input year for the run.
        """
        namelist_file = open(os.path.join(path, 'namelist'), 'w')
        namelist_file.write('&datein \n')
        namelist_file.write('  inp%bcc={0}, \n'.format(str(year)[0:2]))
        namelist_file.write('  inp%byy={0}, \n'.format(str(year)[2:4]))
        namelist_file.write(f'  inp%bmm={month:02d}, \n')
        namelist_file.write('  inp%bdd=01, \n')
        namelist_file.write('  inp%bhh=00, \n')
        namelist_file.write('  inp%model=4, \n')
        namelist_file.write('  inp%lt_units=\'hours\', \n')
        namelist_file.write('  inp%file_seq=\'onebig\', \n')
        namelist_file.write('  inp%modtyp=\'global\', \n')
        namelist_file.write('  inp%nesttyp=\'fixed\', \n')
        namelist_file.write('  inp%filetype=\'netcdf\', \n')
        namelist_file.write('/\n')
        namelist_file.write('&atcfinfo \n')
        namelist_file.write('  atcfnum=81, \n')
        namelist_file.write('  atcfname=\'test\', \n')
        namelist_file.write('  atcfymdh={0}{1:02d}0100, \n'.format(
            year, month))
        namelist_file.write('  atcffreq={0}, \n'.format(self.atcffreq))
        namelist_file.write('/\n')
        namelist_file.write('&trackerinfo \n')
        namelist_file.write('  trkrinfo%westbd={0}, \n'.format(self.westbd))
        namelist_file.write('  trkrinfo%eastbd={0}, \n'.format(self.eastbd))
        namelist_file.write('  trkrinfo%northbd={0}, \n'.format(self.northbd))
        namelist_file.write('  trkrinfo%southbd={0}, \n'.format(self.southbd))
        namelist_file.write('  trkrinfo%type=\'{0}\', \n'.format(
            self.tracktype))
        namelist_file.write('  trkrinfo%mslpthresh={0}, \n'.format(
            self.mslpthresh))
        namelist_file.write('  trkrinfo%choose_t2=\'y\', \n')
        namelist_file.write('  trkrinfo%mslpthresh2={0}, \n'.format(
            self.mslpthresh2))
        namelist_file.write('  trkrinfo%v850thresh={0}, \n'.format(
            self.v850thresh))
        namelist_file.write('  trkrinfo%gridtype=\'global\', \n')
        namelist_file.write('  trkrinfo%contint={0}, \n'.format(self.contint))
        namelist_file.write('  trkrinfo%out_vit=\'y\', \n')
        namelist_file.write('/\n')
        namelist_file.write('&phaseinfo \n')
        namelist_file.write('  phaseflag=\'y\', \n')
        namelist_file.write('  phasescheme=\'both\', \n')
        namelist_file.write('  wcore_depth={0}, \n'.format(self.wcore_depth))
        namelist_file.write('/\n')
        namelist_file.write('&structinfo \n')
        namelist_file.write('  structflag=\'y\', \n')
        namelist_file.write('  ikeflag=\'{0}\', \n'.format(self.ikeflag))
        namelist_file.write('/\n')
        namelist_file.write('&fnameinfo \n')
        namelist_file.write('  gmodname=\'\', \n')
        namelist_file.write('  rundescr=\'\', \n')
        namelist_file.write('  atcfdescr=\'\', \n')
        namelist_file.write('/\n')
        namelist_file.write('&waitinfo \n')
        namelist_file.write('  use_waitfor=\'n\', \n')
        namelist_file.write('/\n')
        namelist_file.write('&verbose \n')
        namelist_file.write('  verb={0}, \n'.format(self.verb))
        namelist_file.write('/')
        namelist_file.close()

    def write_provenance(self, alias, data, output_file):
        """
        Write provenance

        Parameters
        ----------
        alias: str
            Alias of the dataset
        data: dict
            Dictionary containing the dataset information.
        output_file: str
            Path of the output netcdf file.

        """
        ancestors = []
        for i in range(len(data[alias])):
            ancestors.append(data[alias][i]['filename'])
        caption = ("Cyclone tracker output between "
                   "{start} and {end} according to {dataset}").format(
                       start=data[alias][0]['start_year'],
                       end=data[alias][0]['end_year'],
                       dataset=data[alias][0]['dataset']
                   )
        record = {
            'caption': caption,
            'domains': ['global'],
            'authors': ['caron_louis-philippe'],
            'references': ['primavera'],
            'ancestors': ancestors
            }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(output_file, record)


def main():
    """
    Run Cyclone Tracker diagnostic
    """
    with run_diagnostic() as config:
        CycloneTracker(config).compute()


if __name__ == "__main__":
    main()
