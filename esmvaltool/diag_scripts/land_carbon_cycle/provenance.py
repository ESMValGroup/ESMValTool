"""Handle provenance record of land_carbon_cycle diagnostic."""

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    select_metadata,
)


def _get_project(cfg):
    """Extract project from cfg."""
    input_data = cfg['input_data'].values()
    projects = list(group_metadata(input_data, 'project').keys())
    projects = [p for p in projects if 'obs' not in p.lower()]
    if len(projects) == 1:
        return projects[0]
    return projects


def _get_ancestor_files(cfg, obs_name, projects=None):
    """Get ancestor files for provenance."""
    if projects is None:
        projects = _get_project(cfg)
    if isinstance(projects, str):
        projects = [projects]
    datasets = []
    for project in projects:
        datasets.extend(
            select_metadata(cfg['input_data'].values(), project=project))
    datasets.extend(
        select_metadata(cfg['input_data'].values(), dataset=obs_name))
    return [d['filename'] for d in datasets]


def _get_provenance_record(caption, statistics, plot_type, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'ancestors': ancestor_files,
        'authors': ['koirala_sujan'],
        'caption': caption,
        'domains': ['global'],
        'plot_type': plot_type,
        'realms': ['land'],
        'references': ['carvalhais14nature'],
        'statistics': statistics,
        'themes': ['bgchem', 'carbon', 'chem', 'ghg'],
    }
    return record
