'''
This file calculates monsoon metrics using observational data
'''

from . import monsoon_general as mgl
from . import south_asian_monsoon_metrics as sam
from . import east_asian_monsoon_metrics as eam


def eam_obs_area_metrics(run):
    'East Asian Monsoon metrics area averaged for observations'

    # What land sea mask do I use for land only metrics?

    land_jja = None  # mgl.extract_observations(run, 'land', 'jja')
    land_djf = None  # mgl.extract_observations(run, 'land', 'djf')
    temp_jja = mgl.extract_observations(run, 'temp', 'jja')
    temp_djf = mgl.extract_observations(run, 'temp', 'djf')
    ppn_jja = mgl.extract_observations(run, 'precip', 'jja')
    ppn_djf = mgl.extract_observations(run, 'precip', 'djf')
    mslp_jja = mgl.extract_observations(run, 'mslp', 'jja')
    mslp_djf = mgl.extract_observations(run, 'mslp', 'djf')
    z500_jja = mgl.extract_observations(run, 'geopotential', 'jja')
    z500_djf = mgl.extract_observations(run, 'geopotential', 'djf')
    v850_jja = mgl.extract_observations(run, 'Vwind850', 'jja')
    v850_djf = mgl.extract_observations(run, 'Vwind850', 'djf')
    u200_jja = mgl.extract_observations(run, 'Uwind200', 'jja')
    u200_djf = mgl.extract_observations(run, 'Uwind200', 'djf')

    metrics = eam.eam_area_metrics(temp_jja, temp_djf, ppn_jja, ppn_djf,
                                   mslp_jja, mslp_djf, z500_jja, z500_djf,
                                   v850_jja, v850_djf, u200_jja, u200_djf,
                                   land_jja=land_jja, land_djf=land_djf)

    return metrics


def eam_obs_monsoon_indices(run):
    '''
    Outputs the EAM Monsoon Index values for several different indices for
    model.
    '''

    v850_jja = mgl.extract_observations(run, 'Vwind850', 'jja')
    v850_djf = mgl.extract_observations(run, 'Vwind850', 'djf')
    mslp_jja = mgl.extract_observations(run, 'mslp', 'jja')
    mslp_djf = mgl.extract_observations(run, 'mslp', 'djf')
    u850_jja = mgl.extract_observations(run, 'Uwind850', 'jja')
    u300_djf = mgl.extract_observations(run, 'Uwind300', 'djf')

    metrics = eam.eam_monsoon_indices(v850_jja, v850_djf, mslp_jja, mslp_djf,
                                      u850_jja, u300_djf)

    return metrics


def eam_obs_other_metrics(run):
    '''
    Output the EAM metrics which are neither area averages or monsoon metrics
    for model
    '''

    mslp_jja = mgl.extract_observations(run, 'mslp', 'jja')
    mslp_djf = mgl.extract_observations(run, 'mslp', 'djf')
    u850_jja = mgl.extract_observations(run, 'Uwind850', 'jja')
    u850_djf = mgl.extract_observations(run, 'Uwind850', 'djf')
    u200_jja = mgl.extract_observations(run, 'Uwind200', 'jja')
    u200_djf = mgl.extract_observations(run, 'Uwind200', 'djf')
    v850_jja = mgl.extract_observations(run, 'Vwind850', 'jja')
    v850_djf = mgl.extract_observations(run, 'Vwind850', 'djf')

    metrics = eam.eam_other_metrics(mslp_jja, mslp_djf, u850_jja, u850_djf,
                                    u200_jja, u200_djf, v850_jja, v850_djf)

    return metrics


def eam_obs_rmse_metrics(run):
    '''
    Output the EAM metrics which are RMS Errors

    Note that RMSE will always be zero, and correlations will always be one.
    '''

    metrics = dict()

    # RMS(Temp:land)JJA 100-145E 20-50N
    metrics['RMS(Temp:land)JJA 100-145E 20-50N'] = 0.0

    # RMS(Temp:land)DJF 100-145E 20-50N
    metrics['RMS(Temp:land)DJF 100-145E 20-50N'] = 0.0

    # RMS(MSLP) JJA 100-145E 20-50N
    metrics['RMS(MSLP)JJA 100-145E 20-50N'] = 0.0

    # RMS(MSLP) DJF 100-145E 20-50N
    metrics['RMS(MSLP)DJF 100-145E 20-50N'] = 0.0

    # RMS(Prcp:land) JJA 100-145E 20-50N
    metrics['RMS(Prcp:land)JJA 100-145E 20-50N'] = 0.0

    # RMS(Prcp:land) DJF 100-145E 20-50N
    metrics['RMS(Prcp:land)DJF 100-145E 20-50N'] = 0.0

    # RMS(Prcp) JJA 100-145E 20-50N
    metrics['RMS(Prcp)JJA 100-145E 20-50N'] = 0.0

    # RMS(Prcp) DJF 100-145E 20-50N
    metrics['RMS(Prcp)DJF 100-145E 20-50N'] = 0.0

    # Pattern Correlation metrics
    # - These are currently outputting a two by two matrix.

    # Prcp JJA pattern Corr [100,20,145,50]
    metrics['Prcp JJA pattern Corr 100-145E 20-50N'] = 1.0

    # Prcp DJF pattern Corr [100,20,145,50]
    metrics['Prcp DJF pattern Corr 100-145E 20-50N'] = 1.0

    return metrics


def sam_obs_area_metrics(run):
    'South Asian Monsoon metrics area averaged for observations'

    # What land sea mask do I use for land only metrics?

    land_jjas = None  # mgl.extract_observations(run, 'land', 'jjas')
    ppn_jjas = mgl.extract_observations(run, 'precip', 'jjas')
    temp_jjas = mgl.extract_observations(run, 'temp', 'jjas')
    temp_djf = mgl.extract_observations(run, 'temp', 'djf')
    u850_jjas = mgl.extract_observations(run, 'Uwind850', 'jjas')
    u200_jjas = mgl.extract_observations(run, 'Uwind200', 'jjas')
    u10m_son = mgl.extract_observations(run, '10mUwind', 'son')

    metrics = sam.sam_area_metrics(temp_jjas, temp_djf, ppn_jjas,
                                   u850_jjas, u200_jjas, u10m_son,
                                   land_jjas=land_jjas)

    return metrics


def sam_obs_monsoon_indices(run):
    '''Output the SAM Seasonally and area averaged metrics for model'''

    u850_jjas = mgl.extract_observations(run, 'Uwind850', 'jjas')
    u200_jjas = mgl.extract_observations(run, 'Uwind200', 'jjas')
    v850_jjas = mgl.extract_observations(run, 'Vwind850', 'jjas')
    v200_jjas = mgl.extract_observations(run, 'Vwind200', 'jjas')

    metrics = sam.sam_monsoon_indices(u850_jjas, u200_jjas,
                                      v850_jjas, v200_jjas)

    return metrics


def sam_obs_jet_metrics(run):
    '''
    Output the SAM metrics which are specifically related to the size,
    position and shape of the jet for model.
    '''

    u850_as = mgl.extract_observations(run, 'Uwind850', 'as')

    metrics = sam.sam_jet_metrics(u850_as)

    return metrics


def sam_obs_other_metrics(run):
    '''
    Output the SAM metrics which are neither area averages or monsoon metrics

    Note that RMSE will always be zero, and correlations will always be one.
    '''

    metrics = dict()

    # RMSE metrics

    # Precip JJAS RMSE Region = [40,-20,160,50]
    metrics['Precip JJAS RMSE 40-160E 20S-50N'] = 0.0

    # UVWind JJAS RMSE Region = [40,-20,160,50] - Not sure what level????

    # Pattern Correlation metrics

    metrics['Precip JJAS Pattern Corr 40-160E 20S-50N'] = 1.0

    return metrics


if __name__ == '__main__':
    import csv
    from auto_assess.utils import combine_metrics
    from auto_assess.model_run import make_run_dict
    run = make_run_dict(control=True)
    metrics_dict = dict()
    metrics_dict['EAM_obs_area'] = eam_obs_area_metrics(run)
    metrics_dict['EAM_obs_monind'] = eam_obs_monsoon_indices(run)
    metrics_dict['EAM_obs_other'] = eam_obs_other_metrics(run)
    metrics_dict['EAM_obs_rmse'] = eam_obs_rmse_metrics(run)
    metrics_dict['SAM_obs_area'] = sam_obs_area_metrics(run)
    metrics_dict['SAM_obs_monind'] = sam_obs_monsoon_indices(run)
    metrics_dict['SAM_obs_jet'] = sam_obs_jet_metrics(run)
    metrics_dict['SAM_obs_other'] = sam_obs_other_metrics(run)
    metrics = combine_metrics(metrics_dict)

    for metric in sorted(metrics.keys()):
        print '{0} = {1}'.format(metric, metrics[metric])

    with open('monsoon_obs.csv', 'w') as csvfile:
        w = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        for metric in metrics.items():
            w.writerow(metric)
