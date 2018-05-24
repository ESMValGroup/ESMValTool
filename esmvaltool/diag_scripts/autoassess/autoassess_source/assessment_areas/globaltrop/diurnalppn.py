'''
Module docstring
'''

import os

import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
import cf_units
import iris

from ..auto_assess_deprecated.supermeans import get_supermean

import custom_colourbars
import custom_plots
import diurnal_fft

def plot_comparison(runs):
    '''
    Routine to plot diurnal cycle of precipitation diagnostics.

    This is a type 2 multi-function routine that returns no metrics

    Arguments:
        runs - list of run dictionaries.  Each dictionary contains
               metadata for a single model run.  The first dictionary
               in this list is the control experiment.
               (see assessment_area.Area.py_program for description of
               the contents of this dictionary)

    Returns:
        doesn't return any objects - it only writes image files to the
        current working dir
    '''

    run_cntl = runs[0]
    run_expts = runs[1:]
    obsdir = run_cntl['clim_root']

    CNTL_SUPERMEAN_DATA_DIR = os.path.join(
        run_cntl['data_root'],
        run_cntl['runid'],
        run_cntl['_area'] + '_supermeans')

    for run_expt in run_expts:
        cntl = run_cntl.id
        expt = run_expt.id

        EXPT_SUPERMEAN_DATA_DIR = os.path.join(
            run_expt['data_root'],
            run_expt['runid'],
            run_expt['_area'] + '_supermeans')

        # For each season...
        for season in ['djf', 'jja']:

            # Load the files
            cubelist_cntl = get_supermean(
                    'precipitation_flux',
                    season,
                    CNTL_SUPERMEAN_DATA_DIR)  # m01s05i216
            cubelist_cntl = iris.cube.CubeList([cubelist_cntl])
            cubelist_expt = get_supermean('precipitation_flux', season,
                    EXPT_SUPERMEAN_DATA_DIR)  # m01s05i216
            cubelist_expt = iris.cube.CubeList([cubelist_expt])

            # Create plots for precip
            field = 'precip'
            _create_plot(field, season, cubelist_cntl, cubelist_expt, cntl,
                    expt, obsdir)

            # Load the files
            cubelist_cntl = get_supermean('toa_outgoing_longwave_flux', season,
                    CNTL_SUPERMEAN_DATA_DIR)  # m01s03i332
            cubelist_cntl = iris.cube.CubeList([cubelist_cntl])
            cubelist_expt = get_supermean('toa_outgoing_longwave_flux', season,
                    EXPT_SUPERMEAN_DATA_DIR)  # m01s03i332
            cubelist_expt = iris.cube.CubeList([cubelist_expt])

            # Create plots for OLR
            field = 'olr'
            _create_plot(field, season, cubelist_cntl, cubelist_expt, cntl,
                    expt, obsdir)


def _create_plot(field, season, cubelist_cntl, cubelist_test, cntlid, testid,
                 obsdir):
    '''
    '''

    # Construct plot labels
    label_cntl = 'Control : {}'.format(cntlid)
    label_test = 'Test : {}'.format(testid)

    # Setup climatology filename, plot label strings, and model data stash
    # codes; these all differ depending on which field we're plotting.
    if field == 'precip':
        file_clim = os.path.join(obsdir, 'TRMM',
                                 'TRMM_{}_010198_311206.pp'.format(season))
        field_label = 'precipitation / mm day-1'
        label_clim = 'Obs : TRMM'
        cf_name = 'precipitation_flux'  # m01s05i216
    elif field == 'olr':
        file_clim = os.path.join(obsdir, 'CERES',
                                 'CER_SRBAVG1_Terra-MODIS_Edition2D.' + \
                                 '2000-2005_{}.toa_up.pp'.format(season))
        field_label = 'OLR / W m-2'
        label_clim = 'Obs : CERES'
        cf_name = 'toa_outgoing_longwave_flux'  # m01s03i332

    # Create an Iris constraint object to extract latitudes between +/-50 deg
    lat_pm50 = iris.Constraint(latitude=
                               lambda cell: cell > -50.0 and cell <= 50.0)
    # Iris constraint to extract only fields with the desired stash code
    desired_field = iris.Constraint(cf_name)

    # Load the observed climatology file
    cube_clim = iris.load(file_clim, lat_pm50)

    # Extract the desired model data.
    cube_cntl = cubelist_cntl.extract_strict(lat_pm50 & desired_field)
    cube_test = cubelist_test.extract_strict(lat_pm50 & desired_field)

    # If we're looking at precipitation...
    if field == 'precip':

        # The precip climatology file has corrupted time data; set the time
        # coord for each field and merge them into a single cube
        for (i, cube) in enumerate(cube_clim):
            cube.coord('time').points = float(i) * 3.0
        cube_clim = cube_clim.merge()[0]

        # Convert from 3hr accumulations to mm day-1
        cube_clim = cube_clim * 8.0
        # Set basic meta-data that is missing in the file
        cube_clim.units = 'mm day-1'
        cube_clim.var_name = 'precipitation'

        # TODO unit conversion
        # Convert the units of the model data from kg m-2 s-2 to mm day-1
        cube_cntl = cube_cntl * 86400.0
        cube_test = cube_test * 86400.0
        cube_cntl.units = 'mm day-1'
        cube_test.units = 'mm day-1'

    # If we're looking at OLR...
    elif field == 'olr':

        # Extract the first and only cube in the list
        cube_clim = cube_clim[0]
        # The OLR climatology file has longitudes out by 180 degrees!  Fix it
        cube_clim.coord('longitude').points = \
            cube_clim.coord('longitude').points - 180.0
        # Add missing meta-data
        cube_clim.units = 'W m-2'
        cube_clim.var_name = 'OLR'

    # Create a blank cube with coords the same as the UM data, except that
    # the longitudes are shifted so that the prime meridian is in the middle
    coord_cube = iris.cube.Cube(np.zeros(np.shape(cube_cntl)),
                                dim_coords_and_dims=[
                                (cube_cntl.coord('latitude').copy(), 1),
                                (cube_cntl.coord('longitude').copy(), 2)])
    coord_cube.coord('longitude').points = \
        coord_cube.coord('longitude').points - 180.0

    # Regrid / interpolate all the data onto this coordinate set
    cube_clim = iris.analysis.interpolate.regrid(cube_clim, coord_cube)
    cube_cntl = iris.analysis.interpolate.regrid(cube_cntl, coord_cube)
    cube_test = iris.analysis.interpolate.regrid(cube_test, coord_cube)

    if field == 'olr':
        # The OLR climatology data is binned w.r.t. local time (whereas
        # everything else is in UTC).  Set longitudes to a constant zero to
        # prevent the conversion to local time in the phase calculation from
        # adding erroneous dependence on longitude.
        l_zero_lon = True
    else:
        l_zero_lon = False

    # Do the diurnal fft calculations for obs data, control run and test run.
    (amplitude_clim, phase_clim, local_diurnal_hour_clim) = \
        _calc_cube_diurnal_phase(cube_clim, field, l_zerolon=l_zero_lon)
    (amplitude_cntl, phase_cntl, local_diurnal_hour_cntl) = \
        _calc_cube_diurnal_phase(cube_cntl, field)
    (amplitude_test, phase_test, local_diurnal_hour_test) = \
        _calc_cube_diurnal_phase(cube_test, field)

    # Set base part of the filename string for saving plots
    figbase = '{0}_{1}_{2}_{3}'.format(testid, cntlid, field, season)

    # Plot local hour of phase of diurnal harmonic, fading the plot out
    # where the amplitude of the harmonic is small compared to the mean field
    _phase_plot(local_diurnal_hour_clim,
                field, label_clim, '{}_diurnal_phase_clim.png'.format(figbase))
    _phase_plot(local_diurnal_hour_cntl,
                field, label_cntl, '{}_diurnal_phase_cntl.png'.format(figbase))
    _phase_plot(local_diurnal_hour_test,
                field, label_test, '{}_diurnal_phase_test.png'.format(figbase))

    # Plot the difference in phase between the model runs and climatology
    _phase_diff_plot(local_diurnal_hour_clim, local_diurnal_hour_cntl,
                     field, label_cntl + '    minus    ' + label_clim,
                     '{}_diurnal_phase_cntl_minus_clim.png'.format(figbase))
    _phase_diff_plot(local_diurnal_hour_clim, local_diurnal_hour_test,
                     field, label_test + '    minus    ' + label_clim,
                     '{}_diurnal_phase_test_minus_clim.png'.format(figbase))
    _phase_diff_plot(local_diurnal_hour_cntl, local_diurnal_hour_test,
                     field, label_test + '    minus    ' + label_cntl,
                     '{}_diurnal_phase_test_minus_cntl.png'.format(figbase))

    # Plot the amplitude of the diurnal harmonic
    _amp_plot(amplitude_clim[1], field, field_label, label_clim,
              '{}_diurnal_amplitude_clim.png'.format(figbase))
    _amp_plot(amplitude_cntl[1], field, field_label, label_cntl,
              '{}_diurnal_amplitude_cntl.png'.format(figbase))
    _amp_plot(amplitude_test[1], field, field_label, label_test,
              '{}_diurnal_amplitude_test.png'.format(figbase))

    # Plot the difference in amplitude between the model runs and climatology
    _amp_diff_plot(amplitude_cntl[1]-amplitude_clim[1],
                   field, field_label,
                   label_cntl + '    minus    ' + label_clim,
                   '{}_diurnal_amplitude_cntl_minus_clim.png'.format(figbase))
    _amp_diff_plot(amplitude_test[1]-amplitude_clim[1],
                   field, field_label,
                   label_test + '    minus    ' + label_clim,
                   '{}_diurnal_amplitude_test_minus_clim.png'.format(figbase))
    _amp_diff_plot(amplitude_test[1]-amplitude_cntl[1],
                   field, field_label,
                   label_test + '    minus    ' + label_cntl,
                   '{}_diurnal_amplitude_test_minus_cntl.png'.format(figbase))


def _calc_cube_diurnal_phase(cube, field, l_zerolon=False):
    '''
    '''

    # Get the coordinate dimensions
    nt, ny, nx = cube.shape
    # Calculate the length of the frequency domain after real-input FFT
    nf = nt/2 + 1

    # Create a new Iris coord for frequency
    frequency = iris.coords.DimCoord(range(nf), var_name='frequency')

    # Create new Iris cubes for amplitude, phase, and local hour of peak of
    # diurnal harmonic
    amplitude = iris.cube.Cube(np.zeros([nf, ny, nx]),
                               dim_coords_and_dims=[
                               (frequency, 0),
                               (cube.coord('latitude'), 1),
                               (cube.coord('longitude'), 2)],
                               var_name='amplitude',
                               units=cube.units)
    phase = iris.cube.Cube(np.zeros([nf, ny, nx]),
                           dim_coords_and_dims=[
                           (frequency, 0),
                           (cube.coord('latitude'), 1),
                           (cube.coord('longitude'), 2)],
                           var_name='phase',
                           units=cf_units.Unit('radians'))
    local_diurnal_hour = iris.cube.Cube(np.zeros([ny, nx]),
                                        dim_coords_and_dims=[
                                        (cube.coord('latitude'), 0),
                                        (cube.coord('longitude'), 1)],
                                        var_name='local_diurnal_hour',
                                        units=cf_units.Unit('hours'))

    # Calculate amplitude and phase
    (amplitude.data, phase.data) = diurnal_fft.amplitude_and_phase(cube.data)

    # Use the longitude coordinate values to convert the phase of the diurnal
    # harmonic (index 1 in phase frequency space) to local time in hours.
    if l_zerolon:
        lon = np.zeros_like(cube.coord('longitude').points)
    else:
        lon = cube.coord('longitude').points
    local_diurnal_hour.data = (phase.data[1, :, :] + np.tile(lon, (ny, 1))
                               * np.pi/180.0) * (24.0/(2.0*np.pi))

    # Shift the local hour to be in the interval 0 < hour <= 24
    l_gr_24 = local_diurnal_hour.data > 24.0
    local_diurnal_hour.data[l_gr_24] = local_diurnal_hour.data[l_gr_24] - 24.0
    l_le_00 = local_diurnal_hour.data <= 0.0
    local_diurnal_hour.data[l_le_00] = local_diurnal_hour.data[l_le_00] + 24.0

    if field == 'olr':
        # For OLR, shift the phases by 12 hours; this is so that we calculate
        # the local time of minimum for OLR, rather than maximum.  Minimum
        # OLR corresponds to maximum deep convection.
        local_diurnal_hour.data = local_diurnal_hour.data + 12.0
        l_gr_24 = local_diurnal_hour.data > 24.0
        local_diurnal_hour.data[l_gr_24] = local_diurnal_hour.data[l_gr_24] \
            - 24.0

    return amplitude, phase, local_diurnal_hour


def _phase_plot(phase, field, plot_label, plot_file):
    '''
    '''
    if field == 'precip':
        pha_label = 'Local time of maximum in diurnal harmonic of ' +\
                    'precipitation / hour'
    elif field == 'olr':
        pha_label = 'Local time of minimum in diurnal harmonic of ' +\
                    'OLR / hour'
    pha_intervals = np.arange(24)+0.5
    colours = custom_colourbars.cyclic(24)
    plt.figure(figsize=(8, 3.2))
    axeid = plt.axes([0.02, 0.23, 0.96, 0.71], projection=ccrs.PlateCarree())
    qm = custom_plots.pcolor(phase, 'latitude', 'longitude', pha_intervals,
                             colours, ColourbarPos=[0.02, 0.15, 0.96, 0.06])
    plt.setp(qm, transform=ccrs.PlateCarree())
    axeid.set_extent([-180, 180, -50, 50], ccrs.PlateCarree())
    axeid.coastlines()
    axeid.gridlines()
    axeid.set_title(plot_label, loc='left')
    plt.setp(axeid, xticks=0.5+np.arange(0, 24+1, 3))
    plt.setp(axeid, xticklabels=np.asarray(np.arange(0, 24+1, 3), dtype=str))
    plt.xlim([0.5, 24.5])
#     plt.xlabel(pha_label)
    plt.text(0.5, -1.5, pha_label, horizontalalignment='center',
             verticalalignment='top', transform=plt.gca().transAxes)
    plt.draw()
    plt.savefig(plot_file)
    plt.close()


def _phase_diff_plot(phase1, phase2, field, plot_label, plot_file):
    '''
    '''
    if field == 'precip':
        pha_label = 'Difference in local time of maximum in diurnal ' +\
                    'harmonic of precipitation / hour'
    elif field == 'olr':
        pha_label = 'Difference in Local time of minimum in diurnal ' +\
                    'harmonic of OLR / hour'
    phase_difference = phase2 - phase1
    i = np.nonzero(phase_difference.data > 12.0)
    phase_difference.data[i] = phase_difference.data[i] - 24.0
    i = np.nonzero(phase_difference.data < -12.0)
    phase_difference.data[i] = phase_difference.data[i] + 24.0
    pha_intervals = np.concatenate(
        (np.arange(-12, 0), np.array([-0.25, 0.25]), np.arange(1, 12+1)))
    colours = custom_colourbars.diverging_two_sided_cyclic(12)
    plt.figure(figsize=(8, 3.2))
    axeid = plt.axes([0.02, 0.23, 0.96, 0.71], projection=ccrs.PlateCarree())
    qm = custom_plots.pcolor(phase_difference,
                             'latitude', 'longitude', pha_intervals,
                             colours, ColourbarPos=[0.02, 0.15, 0.96, 0.06])
    plt.setp(qm, transform=ccrs.PlateCarree())
    axeid.set_extent([-180, 180, -50, 50], ccrs.PlateCarree())
    axeid.coastlines()
    axeid.gridlines()
    axeid.set_title(plot_label, loc='left')
    plt.setp(axeid, xticks=[1, 4, 7, 10, 12, 13, 14, 15, 17, 20, 23, 26])
    plt.setp(axeid, xticklabels=['-12', '-9', '-6', '-3', '-1',
             r'-$\frac{1}{4}$', r'$\frac{1}{4}$', '1', '3', '6', '9', '12'])
    plt.xlim([1, 26])
#    plt.xlabel(pha_label)
    plt.text(0.5, -1.5, pha_label, horizontalalignment='center',
             verticalalignment='top', transform=plt.gca().transAxes)
    plt.savefig(plot_file)
    plt.close()


def _amp_plot(amplitude, field, field_label, plot_label, plot_file):
    '''
    '''
    if field == 'precip':
        amp_intervals = \
            np.array([0.1, 0.3, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0])
    elif field == 'olr':
        amp_intervals = \
            np.array([1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0])
    n_amp = len(amp_intervals)
    colours = custom_colourbars.diverging_one_sided(n_amp-1, style='blue')
    plt.figure(figsize=(8, 3.2))
    axeid = plt.axes([0.02, 0.23, 0.96, 0.71], projection=ccrs.PlateCarree())
    qm = custom_plots.pcolor(amplitude,
                             'latitude', 'longitude', amp_intervals,
                             colours, ColourbarPos=[0.02, 0.15, 0.96, 0.06])
    plt.setp(qm, transform=ccrs.PlateCarree())
    axeid.set_extent([-180, 180, -50, 50], ccrs.PlateCarree())
    axeid.coastlines()
    axeid.gridlines()
    axeid.set_title(plot_label, loc='left')
    plt.xlim([0.5, n_amp+0.5])
#    plt.xlabel('Amplitude of diurnal harmonic of '+field_label)
    plt.text(0.5, -1.5, 'Amplitude of diurnal harmonic of '+field_label,
             horizontalalignment='center', verticalalignment='top',
             transform=plt.gca().transAxes)
    plt.savefig(plot_file)
    plt.close()


def _amp_diff_plot(amplitude, field, field_label, plot_label, plot_file):
    '''
    '''
    if field == 'precip':
        amp_intervals = \
            np.array([0.1, 0.3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    elif field == 'olr':
        amp_intervals = \
            np.array([1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0])
    n_amp = len(amp_intervals)
    amp_intervals = np.concatenate((-amp_intervals[::-1], amp_intervals))
    colours = custom_colourbars.diverging_two_sided(n_amp-1)[::-1, :]
    plt.figure(figsize=(8, 3.2))
    axeid = plt.axes([0.02, 0.23, 0.96, 0.71], projection=ccrs.PlateCarree())
    qm = custom_plots.pcolor(amplitude,
                             'latitude', 'longitude', amp_intervals,
                             colours, ColourbarPos=[0.02, 0.15, 0.96, 0.06])
    plt.setp(qm, transform=ccrs.PlateCarree())
    axeid.set_extent([-180, 180, -50, 50], ccrs.PlateCarree())
    axeid.coastlines()
    axeid.gridlines()
    axeid.set_title(plot_label, loc='left')
    plt.xlim([0.5, 2*n_amp+0.5])
#    plt.xlabel('Amplitude of diurnal harmonic of '+field_label)
    plt.text(0.5, -1.5,
             'Difference in amplitude of diurnal harmonic of '+field_label,
             horizontalalignment='center', verticalalignment='top',
             transform=plt.gca().transAxes)
    plt.savefig(plot_file)
    plt.close()
