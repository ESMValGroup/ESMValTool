"""Template recipe parser."""

import logging

import yaml
from esmvalcore.preprocessor._derive import get_required


logger = logging.getLogger(__name__)


def extract_requirements(recipe_path):
    """Extract all requirements of a single recipe."""
    requirements = {}
    with open(recipe_path, 'r') as recipe_file:
        recipe = yaml.safe_load(recipe_file)
    for (diag_name, diag_block) in recipe.get('diagnostics', {}).items():
        logger.debug("Found diagnostic '%s'", diag_name)
        var_dict = {}
        for (var_name, var_block) in diag_block.get('variables', {}).items():
            var_block.pop('additional_datasets', None)
            var_block.pop('preprocessor', None)
            var_block.pop('reference_dataset', None)
            var_block.pop('fx_files',None)
            var_block.pop('maps_range',None)
            var_block.pop('diff_range',None)
            var_block.pop('thresholds',None)
            var_block.setdefault('short_name', var_name)
            if 'derive' in var_block:
                derivation_input = get_required(var_block['short_name'],var_block['project'])
                for var in derivation_input:
                    derivation_block = dict(var_block)
                    derivation_block.pop('derive')
                    derivation_block.update(var)
                    var_dict[var.get('short_name',var_name)] = derivation_block
            else:
                var_dict[var_name] = var_block
        requirements[diag_name] = var_dict
    return requirements
