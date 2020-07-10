"""Recipe tests."""
import textwrap
from pathlib import Path

import pytest
import yaml

import esmvaltool
from tests.integration.test_recipes_loading import IDS, RECIPES

REFERENCES_PATH = Path(esmvaltool.__file__).absolute().parent / 'references'


@pytest.mark.parametrize('recipe_file', RECIPES, ids=IDS)
def test_reference_tags(recipe_file):
    """Check bibtex file is added to REFERENCES_PATH."""
    recipe = yaml.safe_load(recipe_file.read_text())
    tags = recipe.get('documentation', {}).get('references', [])
    msg = textwrap.dedent("""
        The tag '{tag}' is mentioned in recipe '{recipe}'.
        However, its reference file '{tag}.bibtex' is not available in {path}.
        Please check instructions on how to add references at
        https://docs.esmvaltool.org/en/latest/community/diagnostic.html#adding-references
    """)
    for tag in tags:
        bibtex_file = REFERENCES_PATH / f'{tag}.bibtex'
        assert bibtex_file.is_file(), msg.format(tag=tag,
                                                 recipe=recipe_file,
                                                 path=REFERENCES_PATH)
