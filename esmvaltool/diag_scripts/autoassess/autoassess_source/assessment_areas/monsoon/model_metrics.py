'''
This file calculates monsoon metrics using model data
'''

from . import monsoon_general as mgl
from . import south_asian_monsoon_metrics as sam
from . import east_asian_monsoon_metrics as eam


def eam_model_area_metrics(run):
    'East Asian Monsoon metrics area averaged for model'

    land_jja = mgl.extract_model_season(run, 'm01s03i395', 'jja')
    land_djf = mgl.extract_model_season(run, 'm01s03i395', 'djf')
    temp_jja = mgl.extract_model_season(run, 'm01s03i236', 'jja')
    temp_djf = mgl.extract_model_season(run, 'm01s03i236', 'djf')
    ppn_jja = mgl.extract_model_season(run, 'm01s05i216', 'jja')
    ppn_djf = mgl.extract_model_season(run, 'm01s05i216', 'djf')
    mslp_jja = mgl.extract_model_season(run, 'm01s16i222', 'jja')
    mslp_djf = mgl.extract_model_season(run, 'm01s16i222', 'djf')
    z500_jja = mgl.extract_model_season(run, 30207500, 'jja')
    z500_djf = mgl.extract_model_season(run, 30207500, 'djf')
    v850_jja = mgl.extract_model_season(run, 30202850, 'jja')
    v850_djf = mgl.extract_model_season(run, 30202850, 'djf')
    u200_jja = mgl.extract_model_season(run, 30201200, 'jja')
    u200_djf = mgl.extract_model_season(run, 30201200, 'djf')

    metrics = eam.eam_area_metrics(temp_jja, temp_djf, ppn_jja, ppn_djf,
                                   mslp_jja, mslp_djf, z500_jja, z500_djf,
                                   v850_jja, v850_djf, u200_jja, u200_djf,
                                   land_jja=land_jja, land_djf=land_djf)

    return metrics


def eam_model_monsoon_indices(run):
    '''
    Outputs the EAM Monsoon Index values for several different indices for
    model.
    '''

    v850_jja = mgl.extract_model_season(run, 30202850, 'jja')
    v850_djf = mgl.extract_model_season(run, 30202850, 'djf')
    mslp_jja = mgl.extract_model_season(run, 'm01s16i222', 'jja')
    mslp_djf = mgl.extract_model_season(run, 'm01s16i222', 'djf')
    u850_jja = mgl.extract_model_season(run, 30201850, 'jja')
    u300_djf = mgl.extract_model_season(run, 30201300, 'djf')

    metrics = eam.eam_monsoon_indices(v850_jja, v850_djf, mslp_jja, mslp_djf,
                                      u850_jja, u300_djf)

    return metrics


def eam_model_other_metrics(run):
    '''
    Output the EAM metrics which are neither area averages or monsoon metrics
    for model
    '''

    mslp_jja = mgl.extract_model_season(run, 'm01s16i222', 'jja')
    mslp_djf = mgl.extract_model_season(run, 'm01s16i222', 'djf')
    u850_jja = mgl.extract_model_season(run, 30201850, 'jja')
    u850_djf = mgl.extract_model_season(run, 30201850, 'djf')
    u200_jja = mgl.extract_model_season(run, 30201200, 'jja')
    u200_djf = mgl.extract_model_season(run, 30201200, 'djf')
    v850_jja = mgl.extract_model_season(run, 30202850, 'jja')
    v850_djf = mgl.extract_model_season(run, 30202850, 'djf')

    metrics = eam.eam_other_metrics(mslp_jja, mslp_djf, u850_jja, u850_djf,
                                    u200_jja, u200_djf, v850_jja, v850_djf)

    return metrics


def sam_model_area_metrics(run):
    '''Output the SAM Seasonally and area averaged metrics for model'''

    jjas = [6, 7, 8, 9]
    land_jjas = mgl.extract_model_month(run, 'm01s03i395', jjas)
    ppn_jjas = mgl.extract_model_month(run, 'm01s05i216', jjas)
    temp_jjas = mgl.extract_model_month(run, 'm01s03i236', jjas)
    temp_djf = mgl.extract_model_season(run, 'm01s03i236', 'djf')
    u850_jjas = mgl.extract_model_month(run, 30201850, jjas)
    u200_jjas = mgl.extract_model_month(run, 30201200, jjas)
    u10m_son = mgl.extract_model_season(run, 'm01s03i225', 'son')
    vmfl_jjas = mgl.extract_model_month(run, 'm01s30i428', jjas)

    metrics = sam.sam_area_metrics(temp_jjas, temp_djf, ppn_jjas,
                                   u850_jjas, u200_jjas, u10m_son,
                                   vmfl_jjas=vmfl_jjas, land_jjas=land_jjas)

    return metrics


def sam_model_monsoon_indices(run):
    '''Output the SAM Seasonally and area averaged metrics for model'''

    jjas = [6, 7, 8, 9]
    u850_jjas = mgl.extract_model_month(run, 30201850, jjas)
    u200_jjas = mgl.extract_model_month(run, 30201200, jjas)
    v850_jjas = mgl.extract_model_month(run, 30202850, jjas)
    v200_jjas = mgl.extract_model_month(run, 30202200, jjas)

    metrics = sam.sam_monsoon_indices(u850_jjas, u200_jjas,
                                      v850_jjas, v200_jjas)

    return metrics


def sam_model_jet_metrics(run):
    '''
    Output the SAM metrics which are specifically related to the size,
    position and shape of the jet for model.
    '''

    u850_as = mgl.extract_model_month(run, 30201850, [8, 9])

    metrics = sam.sam_jet_metrics(u850_as)

    return metrics
