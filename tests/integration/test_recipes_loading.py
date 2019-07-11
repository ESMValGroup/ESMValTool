"""Test recipes are well formed."""

import glob
import os
import pytest

from esmvalcore._recipe import load_raw_recipe


def _get_recipes():
    recipes_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..',
        '..',
        'esmvaltool',
        'recipes',
    )
    recipes_path = os.path.abspath(recipes_path)
    print(recipes_path)
    recipes = glob.glob(os.path.join(recipes_path, '*.yml'))
    recipes += glob.glob(os.path.join(recipes_path, 'examples', '*.yml'))
    return recipes

@pytest.mark.parametrize('recipe_file', _get_recipes())
def test_diagnostic_run(recipe_file):
    """Check that recipe files are well formed."""
    load_raw_recipe(recipe_file)
