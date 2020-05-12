"""ATTRIBUTES FOR PROVENANCE TAGGING.

Module containing functions to create the metadata for the output files.

The module contains the following functions:
- get_prov_map: create a record of metadata for 2D outputs;
- get_prov_transp: create a record of metadata for 1D outputs (e.g. transports)
- meta_direntr: write metadata to a file containing one of the components of
                the material entropy production with the direct method;
- meta_indentr: write metadata to a file containing one of the components of
                the material entropy production with the indirect method;

@author: Valerio Lembo, University of Hamburg, 2019.
"""

import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared import ProvenanceLogger


def get_prov_map(attr, ancestor_files):
    """Create a provenance record for the 2D diagnostic outputs."""
    caption = (
        "Thermodynamic Diagnostic Tool - Monthly mean {} (lat, lon) fields"
        "for model {}.".format(attr[0], attr[1]))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': ['geo'],
        'authors': ['lembo_valerio'],
        'references': ['lembo16climdyn', 'lembo19gmdd', 'lucarini14revgeop'],
        'ancestors': ancestor_files,
    }
    return record


def get_prov_transp(attr, ancestor_file, plotname):
    """Create a provenance record for the 1d meridional transports."""
    caption = ("Thermodynamic Diagnostic Tool - Annual mean zonally averaged"
               " meridional {} transports"
               " for model {}.".format(attr[0], attr[1]))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': ['sect'],
        'plot_file': plotname,
        'authors': ['lembo_valerio'],
        'references': ['lembo16climdyn', 'lembo19gmdd', 'lucarini14revgeop'],
        'ancestors': ancestor_file,
    }
    return record


def meta_direntr(cfg, model, input_data, flist):
    """Write metadata to components of the direct entropy prod maps.

    Arguments:r
    - model: the name of the model;
    - inlist: the list of the input filenames;
    - flist: the list of the entropy filenames;

    @author: Valerio Lembo, University of Hamburg, 2019.
    """
    with ProvenanceLogger(cfg) as provlog:
        hfls_file = e.select_metadata(input_data,
                                      short_name='hfls',
                                      dataset=model)[0]['filename']
        hfss_file = e.select_metadata(input_data,
                                      short_name='hfss',
                                      dataset=model)[0]['filename']
        hus_file = e.select_metadata(input_data,
                                     short_name='hus',
                                     dataset=model)[0]['filename']
        pr_file = e.select_metadata(input_data, short_name='pr',
                                    dataset=model)[0]['filename']
        prsn_file = e.select_metadata(input_data,
                                      short_name='prsn',
                                      dataset=model)[0]['filename']
        ps_file = e.select_metadata(input_data, short_name='ps',
                                    dataset=model)[0]['filename']
        rlut_file = e.select_metadata(input_data,
                                      short_name='rlut',
                                      dataset=model)[0]['filename']
        ts_file = e.select_metadata(input_data, short_name='ts',
                                    dataset=model)[0]['filename']
        uas_file = e.select_metadata(input_data,
                                     short_name='uas',
                                     dataset=model)[0]['filename']
        vas_file = e.select_metadata(input_data,
                                     short_name='vas',
                                     dataset=model)[0]['filename']
        attr = ['sensible heat entropy production', model]
        ancestor = [
            hfss_file, hus_file, ps_file, rlut_file, uas_file, vas_file,
            ts_file
        ]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[0], record)
        attr = ['evaporation entropy production', model]
        ancestor = [hfls_file, ts_file]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[1], record)
        attr = ['rainfall precipitation entropy production', model]
        ancestor = [hus_file, pr_file, prsn_file, ps_file, ts_file]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[2], record)
        attr = ['snowfall precipitation entropy production', model]
        ancestor = [hus_file, prsn_file, ps_file, ts_file]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[3], record)
        attr = ['snow melt entropy production', model]
        ancestor = [prsn_file, ts_file]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[4], record)
        attr = ['potential energy entropy production', model]
        ancestor = [hus_file, pr_file, prsn_file, ps_file, rlut_file, ts_file]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[5], record)


def meta_indentr(cfg, model, input_data, flist):
    """Write metadata to components of the indirect entropy prod maps.

    Arguments:
    ---------
    - model: the name of the model;
    - inlist: the list of the input filenames;
    - flist: the list of the entropy filenames;

    @author: Valerio Lembo, University of Hamburg, 2019.
    """
    with ProvenanceLogger(cfg) as provlog:
        rlds_file = e.select_metadata(input_data,
                                      short_name='rlds',
                                      dataset=model)[0]['filename']
        rlus_file = e.select_metadata(input_data,
                                      short_name='rlus',
                                      dataset=model)[0]['filename']
        rlut_file = e.select_metadata(input_data,
                                      short_name='rlut',
                                      dataset=model)[0]['filename']
        rsds_file = e.select_metadata(input_data,
                                      short_name='rsds',
                                      dataset=model)[0]['filename']
        rsdt_file = e.select_metadata(input_data,
                                      short_name='rsdt',
                                      dataset=model)[0]['filename']
        rsus_file = e.select_metadata(input_data,
                                      short_name='rsus',
                                      dataset=model)[0]['filename']
        rsut_file = e.select_metadata(input_data,
                                      short_name='rsut',
                                      dataset=model)[0]['filename']
        ts_file = e.select_metadata(input_data, short_name='ts',
                                    dataset=model)[0]['filename']
        attr = ['horizontal entropy production', model]
        ancestor = [rlut_file, rsdt_file, rsut_file]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[0], record)
        attr = ['vertical entropy production', model]
        ancestor = [
            rlds_file, rlus_file, rlut_file, rsds_file, rsus_file, ts_file
        ]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[1], record)
