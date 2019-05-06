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
        'authors': ['lemb_va'],
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
        'authors': ['lemb_va'],
        'references': ['lembo16climdyn', 'lembo19gmdd', 'lucarini14revgeop'],
        'ancestors': ancestor_file,
    }
    return record


def meta_direntr(cfg, model, inlist, flist):
    """Write metadata to components of the direct entropy prod maps.

    Arguments:r
    - model: the name of the model;
    - inlist: the list of the input filenames;
    - flist: the list of the entropy filenames;

    @author: Valerio Lembo, University of Hamburg, 2019.
    """
    with ProvenanceLogger(cfg) as provlog:
        attr = ['sensible heat entropy production', model]
        ancestor = [
            inlist[1], inlist[2], inlist[5], inlist[15], inlist[17], inlist[19]
            ]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[0], record)
        attr = ['evaporation entropy production', model]
        ancestor = [inlist[0], inlist[2], inlist[5], inlist[15]]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[1], record)
        attr = ['rainfall precipitation entropy production', model]
        ancestor = [inlist[2], inlist[3], inlist[5], inlist[15]]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[2], record)
        attr = ['snowfall precipitation entropy production', model]
        ancestor = [inlist[2], inlist[4], inlist[5], inlist[15]]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[3], record)
        attr = ['snow melt entropy production', model]
        ancestor = [inlist[4], inlist[15]]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[4], record)
        attr = ['potential energy entropy production', model]
        ancestor = [
            inlist[2], inlist[3], inlist[4], inlist[5], inlist[8], inlist[15]
            ]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[5], record)


def meta_indentr(cfg, model, inlist, flist):
    """Write metadata to components of the indirect entropy prod maps.

    Arguments:
    - model: the name of the model;
    - inlist: the list of the input filenames;
    - flist: the list of the entropy filenames;

    @author: Valerio Lembo, University of Hamburg, 2019.
    """
    with ProvenanceLogger(cfg) as provlog:
        attr = ['horizontal entropy production', model]
        ancestor = [inlist[8], inlist[10], inlist[12]]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[0], record)
        attr = ['vertical entropy production', model]
        ancestor = [
            inlist[6], inlist[7], inlist[8], inlist[9], inlist[11], inlist[15]
            ]
        record = get_prov_map(attr, ancestor)
        provlog.log(flist[1], record)
