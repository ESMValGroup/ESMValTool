import copy
import logging
import os
from abc import abstractmethod

from prov.dot import prov_to_dot
from prov.model import ProvDocument

logger = logging.getLogger(__name__)


def get_recipe_provenance(documentation):
    """Create a provenance document describing a recipe run."""
    doc = ProvDocument()
    for nsp in ('recipe', 'author', 'attributes'):
        doc.add_namespace(nsp, 'http://www.esmvaltool.org/' + nsp)

    recipe = doc.entity(
        'recipe:recipe', {
            'recipe:description': documentation['description'],
            'recipe:projects': ', '.join(documentation['projects']),
            'recipe:references': ', '.join(documentation['references']),
        })

    for author in documentation['authors']:
        author = doc.agent(
            'author:' + author['name'],
            {'attributes:' + k: author[k]
             for k in author if k != 'name'})
        doc.wasAttributedTo(recipe, author)

    return doc


class TrackedFile(object):
    def __init__(self, filename, settings, ancestors=None):

        self._filename = filename
        self.settings = copy.deepcopy(settings)

        self.provenance = None
        self.entity = None
        self._ancestors = [] if ancestors is None else ancestors
        self._activity = None

    @property
    def filename(self):
        return self._filename

    def __hash__(self):
        return hash(self.filename)

    def initialize_provenance(self, task):
        """Initialize the provenance document."""
        if self.provenance is not None:
            raise ValueError(
                "Provenance of {} already initialized".format(self))
        self.provenance = ProvDocument()
        self._initialize_namespaces()
        self._initialize_activity(task)
        self._initialize_entity()
        self._initialize_ancestors()

    def _initialize_namespaces(self):
        """Inialize the namespaces."""
        for nsp in ('file', 'attribute', 'preprocessor', 'diagnostic', 'task'):
            self.provenance.add_namespace(nsp,
                                          'http://www.esmvaltool.org/' + nsp)

    def _initialize_activity(self, task):
        """Initialize the preprocessor task activity."""
        self._activity = self.provenance.activity('task:' + task.name)

    def _initialize_ancestors(self):
        """Register input Products/files for provenance tracking."""
        for ancestor in self._ancestors:
            if ancestor.provenance is None:
                raise ValueError("Uninitalized ancestor provenance.")
            self.provenance.update(ancestor.provenance)
            self.wasderivedfrom(ancestor)

    @abstractmethod
    def _initialize_entity(self):
        """Initialize the entity representing the file."""

    def wasderivedfrom(self, other):
        """Let the file know that it was derived from other."""
        entity = other.entity if hasattr(other, 'entity') else other
        if not self._activity:
            raise ValueError("Activity not initialized.")
        self.entity.wasDerivedFrom(entity, self._activity)

    def save_provenance(self):
        """Export provenance information."""
        filename = os.path.splitext(self.filename)[0]
        self.provenance.serialize(filename + '.xml', format='xml')
        figure = prov_to_dot(self.provenance)
        figure.write_png(filename + '.png')


def write_provenance(provenance, output_dir):
    """Write provenance information to output_dir."""
    filename = os.path.join(output_dir, 'provenance')
    logger.info("Writing provenance to %s.xml", filename)
    provenance.serialize(filename + '.xml', format='xml')

    graph = prov_to_dot(provenance)
    logger.info("Writing provenance to %s.png", filename)
    graph.write_png(filename + '.png')
    logger.info("Writing provenance to %s.pdf", filename)
    graph.write_pdf(filename + '.pdf')
