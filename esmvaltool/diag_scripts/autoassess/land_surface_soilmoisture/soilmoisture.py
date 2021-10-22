"""Run module for soil moisture metrics."""

import os
import logging
import numpy as np
import iris
from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger
from esmvaltool.diag_scripts.shared._supermeans import get_supermean

logger = logging.getLogger(__name__)


def get_provenance_record(caption, run):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'metrics',
        'authors': [
            'rumbold_heather',
            'sellar_alistair',
        ],
        'references': [
            'esacci-soilmoisture',
            'dorigo17rse',
            'gruber19essd',
        ],
        'ancestors': run,
    }

    return record


def land_sm_top(run):
    """
    Calculate median absolute errors for soil mosture against CCI data.

    Parameters
    ----------
    run: dict
        dictionary containing model run metadata
        (see auto_assess/model_run.py for description)

    Returns
    -------
    metrics: dict
        a dictionary of metrics names and values
    """
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')

    seasons = ['djf', 'mam', 'jja', 'son']

    # Constant: density of water
    rhow = 1000.

    # Work through each season
    metrics = dict()
    for season in seasons:
        fname = 'ecv_soil_moisture_{}.nc'.format(season)
        clim_file = os.path.join(run['climfiles_root'], fname)
        ecv_clim = iris.load_cube(clim_file)
        # correct invalid units
        if (ecv_clim.units == 'unknown' and
                'invalid_units' in ecv_clim.attributes):
            if ecv_clim.attributes['invalid_units'] == 'm^3m^-3':
                ecv_clim.units = 'm3 m-3'

        # m01s08i223
        # CMOR name: mrsos (soil moisture in top model layer kg/m2)
        mrsos = get_supermean('mass_content_of_water_in_soil_layer',
                              season,
                              supermean_data_dir)

        # Set soil moisture to missing data on ice points (i.e. no soil)
        np.ma.masked_where(mrsos.data == 0, mrsos.data, copy=False)

        # first soil layer depth
        dz1 = mrsos.coord('depth').bounds[0, 1] - \
            mrsos.coord('depth').bounds[0, 0]

        # Calculate the volumetric soil moisture in m3/m3
        # volumetric soil moisture = volume of water / volume of soil layer
        # = depth equivalent of water / thickness of soil layer
        # = (soil moisture content (kg m-2) / water density (kg m-3) )  /
        #      soil layer thickness (m)
        # = mosrs / (rhow * dz1)
        vol_sm1_run = mrsos / (rhow * dz1)
        vol_sm1_run.units = "m3 m-3"
        vol_sm1_run.long_name = "Top layer Soil Moisture"

        # update the coordinate system ECV data with a WGS84 coord system
        # unify coord systems for regridder
        vol_sm1_run.coord('longitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        vol_sm1_run.coord('latitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        ecv_clim.coord('longitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        ecv_clim.coord('latitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)

        # Interpolate to the grid of the climatology and form the difference
        vol_sm1_run = regrid(vol_sm1_run, ecv_clim, 'linear')

        # mask invalids
        vol_sm1_run.data = np.ma.masked_invalid(vol_sm1_run.data)
        ecv_clim.data = np.ma.masked_invalid(ecv_clim.data)

        # diff the cubes
        dff = vol_sm1_run - ecv_clim

        # save output and populate metric
        iris.save(dff, os.path.join(run['dump_output'],
                                    'soilmoist_diff_{}.nc'.format(season)))
        name = 'soilmoisture MedAbsErr {}'.format(season)
        dffs = dff.data
        dffs = np.ma.abs(dffs)
        metrics[name] = float(np.ma.median(dffs))

    # record provenance
    plot_file = "Autoassess soilmoisture metrics"
    caption = 'Autoassess soilmoisture MedAbsErr for {}'.format(str(seasons))
    provenance_record = get_provenance_record(caption, run)
    cfg = {}
    cfg['run_dir'] = run['out_dir']
    # avoid rewriting provenance when running the plot diag
    if not os.path.isfile(os.path.join(cfg['run_dir'],
                                       'diagnostic_provenance.yml')):
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(plot_file, provenance_record)

    return metrics
