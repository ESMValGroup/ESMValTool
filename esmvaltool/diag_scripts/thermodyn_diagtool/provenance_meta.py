"""ATRRIBUTES FOR PROVENANCE TAGGING.

Module containing functions to tag plots from Thermodyn_diagtool.

The module provides plots for a single model of:
- climatological mean maps of TOA, atmospheric and surface energy budgets;
- annual mean time series of TOA, atmospheric and surface energy budgets anom.;
- climatological mean maps of latent energy and water mass budgets;
- annual mean time series of latent energy and water mass budget anom.;
- meridional section of meridional enthalpy transports;
- meridional section of meridional water mass transports;
- scatter plots of atmospheric vs. oceani peak magnitudes in the two hem.;
- climatological mean maps of every component of the entropy budget.

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
        'plot_type': '2D maps',
        'authors': 'lemb_va',
        'references': ['lembo16climdyn', 'lucarini14revgeop'],
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
        'plot_type': 'meridional sections',
        'plot_file': plotname,
        'authors': 'lemb_va',
        'references': ['lembo16climdyn', 'lucarini14revgeop'],
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
    provlog = ProvenanceLogger(cfg)
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
    """Write metadata to components of the inddirect entropy prod maps.

    Arguments:
    - model: the name of the model;
    - inlist: the list of the input filenames;
    - flist: the list of the entropy filenames;

    @author: Valerio Lembo, University of Hamburg, 2019.
    """
    provlog = ProvenanceLogger(cfg)
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
