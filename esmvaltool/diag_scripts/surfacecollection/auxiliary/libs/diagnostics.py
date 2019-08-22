#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:05:25 2019

@author: bmueller
"""
import iris
import logging
import os
from .utilities import set_metadata, checked_ref, data_correlation
from .utilities import calculate_anomalies, calculate_correlation
from .utilities import calculate_trend, corr_extract, unify_cubes
from .utilities import delete_aux_coords, calculate_climatology
from .utilities import weighted_spatial_mean, read_regions
from .utilities import get_masked_from_shape, timer
import pandas as pd
import numpy as np

logger = logging.getLogger(os.path.basename(__file__))

def time_series(data, **kwargs):
    """
    produces spatially aggregated mean time series
    ----------------------------------------------
    returns a combined time series cube
    """
    cubes  = []
        
    for cube in data.get_all():
        return_cube = weighted_spatial_mean(cube)
        cubes.append(return_cube)
        
    glob_cube = unify_cubes(cubes)
    
    return_dict = {"global": glob_cube}
    
    regions = kwargs.pop("subregions", None)
    if regions is not None:
        regions = read_regions(regions)
    
    region_names = np.array(regions.records())[:,kwargs.pop("names_column", 0)]
    
    for indx, shape in enumerate(regions.shapes()):
        loc_poly = shape.points
        masked_data = get_masked_from_shape(data, loc_poly)
        masked_cubes = []
        for cube in masked_data.get_all():
            return_cube = weighted_spatial_mean(cube)
            masked_cubes.append(return_cube)
        loc_cube = unify_cubes(masked_cubes)
        return_dict.update({region_names[indx]:loc_cube})
    
    return [return_dict]

def glob_temp_mean(data, **kwargs):
    """
    produces global temporal mean
    -----------------------------
    returns a list of mean cubes
    """
    cubes = []
    
    for cube in data.get_all():
        delete_aux_coords(cube)
        cubes.append(cube.collapsed("time", iris.analysis.MEAN))
    
    return cubes


def glob_temp_mean_absdiff(data, **kwargs):
    """
    produces global temporal mean absolute differences
    --------------------------------------------------
    returns a list of difference cubes
    """
    cubes = []
    
    ref = glob_temp_mean(data.ref_only())
    ref = checked_ref(ref, num=1)
    nonref = glob_temp_mean(data.nonref_only())

    for cube in nonref:
        diff = cube - ref
        diff.metadata = set_metadata(cube, ref,
                                     "global_temp_mean_absdiff")
        cubes.append(diff)
    
    return cubes


def glob_temp_mean_reldiff(data, **kwargs):
    """
    produces global temporal mean relative differences
    --------------------------------------------------
    returns a list of mean cubes
    """
    cubes = []
    
    ref = glob_temp_mean(data.ref_only())
    ref = checked_ref(ref, num=1)

    for cube in glob_temp_mean_absdiff(data):
        reldiff = cube / ref
        reldiff.metadata = set_metadata(cube, ref,
                                        "global_temp_mean_reldiff")
        cubes.append(reldiff)
    
    return cubes


def percentiles(data, **kwargs):
    """
    produces pixelwise temporal percentiles
    ---------------------------------------
    returns a dict of percentile cubes
    """
    
    if "percentiles" not in kwargs.keys():
        logger.error("option percentiles required for this diagnostic")
    else:
        percentiles = kwargs["percentiles"]
        if len(percentiles) == 0:
            percentiles = [0.5]
            logger.warning("no percentiles given (None), " +
                           "median (0.5) produced instead")
    
    ref = checked_ref(data.get_ref(), num=1)
    
    nonref = data.get_nonref()
    delete_aux_coords(ref)
    refperc = ref.collapsed("time", iris.analysis.PERCENTILE,
                            percent=[p*100 for p in percentiles])
    
    nonrefperc = []
    
    for cube in nonref:
        delete_aux_coords(cube)
        nonrefperc.append(cube.collapsed("time", iris.analysis.PERCENTILE,
                                         percent=[p*100 for p in percentiles]))
        
    res_list = []
    
    for p in percentiles:
        
        perc_check = lambda perc: perc==p*100
        perc_constraint = iris.Constraint(percentile_over_time=perc_check)
        
        refp = refperc.extract(perc_constraint)
        nonrefp = [nrp.extract(perc_constraint) for nrp in nonrefperc]
        
        corrs = [data_correlation(refp, nrp) for nrp in nonrefp]
        
        res_list.append({str(p):{"ref":refp,
                                 "nonref":nonrefp,
                                 "corr":corrs}})
    
    corr_data = corr_extract(res_list)
    
    corr_df = pd.DataFrame(data=corr_data,
                columns=[dset.metadata.attributes["source_file"].split(
                        os.sep)[-1] for dset in data.get_all()],
                index=percentiles)
    
    res_list.append(corr_df)
    
    return res_list

def trend(data, **kwargs):
    """
    produces pixelwise trends
    -------------------------
    returns a list of trend related cubes
    """
    
    pthreshold = kwargs.pop("pthreshold", None)
    
    cubes = []
    
    decadal = kwargs.pop("decadal", True)
    for cube in data.get_all():
        cubes.append(calculate_trend(cube, decadal = decadal))
        
    if pthreshold is not None:
        for c in cubes:
            thresholded = c["trend"].copy()
            thresholded.attributes.update({"pthreshold": pthreshold})
            iris.util.mask_cube(thresholded, c["p-value"].data > pthreshold)
            c.update({"threshold": thresholded})
    
    return cubes

def anomalytrend(data, **kwargs):
    """
    produces pixelwise trends for anomalies
    ---------------------------------------
    returns a list of trend related cubes
    """
    
    if "temporal_basis" not in kwargs.keys():
        logger.error("option temporal_basis required for this diagnostic")
    else:
        temporal_basis = kwargs["temporal_basis"]
        if len(temporal_basis) == 0:
            temporal_basis = "month"
            logger.warning("no temporal_basis given (None), " +
                           "monthly climatology produced instead")
            
    anomalies = calculate_anomalies(data, temporal_basis)
    
    kwargs.update({"decadal": True})
    
    cubes = trend(anomalies, **kwargs)
    
    return cubes

def climatologytrend(data, **kwargs):
    """
    produces pixelwise trends for climatologies
    -------------------------------------------
    returns a list of trend related cubes
    """
    
    if "temporal_basis" not in kwargs.keys():
        logger.error("option temporal_basis required for this diagnostic")
    else:
        temporal_basis = kwargs["temporal_basis"]
        if len(temporal_basis) == 0:
            temporal_basis = "month"
            logger.warning("no temporal_basis given (None), " +
                           "monthly climatology produced instead")
            
    climatology = calculate_climatology(data, temporal_basis)
    
    kwargs.update({"decadal": False})
    
    cubes = trend(climatology, **kwargs)
    
    return cubes

def correlation(data, **kwargs):
    """
    produces pixelwise correlations
    -------------------------------
    returns a list of correlation related cubes
    """
    
    pthreshold = kwargs.pop("pthreshold", None)
    
    ref = checked_ref(data.get_ref(), num=1)
    
    cubes = []
    
    for cube in data.get_nonref():
        cubes.append(calculate_correlation(cube, ref))
        
    if pthreshold is not None:
        for c in cubes:
            thresholded = c["correlation"].copy()
            thresholded.attributes.update({"pthreshold": pthreshold})
            iris.util.mask_cube(thresholded, c["p-value"].data > pthreshold)
            c.update({"threshold": thresholded})
    
    return cubes

def anomalycorrelation(data, **kwargs):
    """
    produces pixelwise correlations for anomalies
    ---------------------------------------------
    returns a list of correlation related cubes
    """
    
    if "temporal_basis" not in kwargs.keys():
        logger.error("option temporal_basis required for this diagnostic")
    else:
        temporal_basis = kwargs["temporal_basis"]
        if len(temporal_basis) == 0:
            temporal_basis = "month"
            logger.warning("no temporal_basis given (None), " +
                           "monthly climatology produced instead")
    
    anomalies = calculate_anomalies(data, temporal_basis)
    
    cubes = correlation(anomalies, **kwargs)
    
    return cubes

def climatologycorrelation(data, **kwargs):
    """
    produces pixelwise correlations for climatologies
    -------------------------------------------------
    returns a list of correlation related cubes
    """
    
    if "temporal_basis" not in kwargs.keys():
        logger.error("option temporal_basis required for this diagnostic")
    else:
        temporal_basis = kwargs["temporal_basis"]
        if len(temporal_basis) == 0:
            temporal_basis = "month"
            logger.warning("no temporal_basis given (None), " +
                           "monthly climatology produced instead")
    
    climatology = calculate_climatology(data, temporal_basis)
    
    cubes = correlation(climatology, **kwargs)
    
    return cubes