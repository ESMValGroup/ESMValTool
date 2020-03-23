"""Recipe tests."""
from pathlib import Path
import pytest
import yaml

from tests.integration.test_recipes_loading import RECIPES, IDS
import esmvaltool

REFERENCES_PATH = Path(esmvaltool.__file__).absolute().parent / 'references'


@pytest.mark.parametrize('recipe_file', RECIPES, ids=IDS)
def _test_reference_tags(recipe_file):
    recipe = yaml.safe_load(recipe_file.read_text())
    value = recipe['documentation']['references']
    _check_bibtex_exist(value)


def _check_bibtex_exist(tag):
    """Check bibtex file is added to REFERENCES_PATH."""
    bibtex_file = REFERENCES_PATH / f'{tag}.bibtex'
    assert not bibtex_file.is_file()
