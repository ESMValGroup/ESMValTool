"""Provenance module."""
import copy
import logging
import os

from netCDF4 import Dataset
from PIL import Image
from PIL.PngImagePlugin import PngInfo
from prov.dot import prov_to_dot
from prov.model import ProvDocument

from ._config import replace_tags
from ._version import __version__

logger = logging.getLogger(__name__)

ESMVALTOOL_URI_PREFIX = 'https://www.esmvaltool.org/'


def update_without_duplicating(bundle, other):
    """Add new records from other provenance bundle."""
    for record in other.records:
        if record not in bundle.records:
            bundle.add_record(record)


def create_namespace(provenance, namespace):
    """Create an esmvaltool namespace."""
    provenance.add_namespace(namespace, uri=ESMVALTOOL_URI_PREFIX + namespace)


def get_esmvaltool_provenance():
    """Create an esmvaltool run activity."""
    provenance = ProvDocument()
    namespace = 'software'
    create_namespace(provenance, namespace)
    attributes = {}  # TODO: add dependencies with versions here
    activity = provenance.activity(
        namespace + ':esmvaltool==' + __version__, other_attributes=attributes)

    return activity


ESMVALTOOL_PROVENANCE = get_esmvaltool_provenance()


def attribute_to_authors(entity, author_tags):
    """Attribute entity to authors."""
    namespace = 'author'
    create_namespace(entity.bundle, namespace)

    authors = replace_tags('authors', author_tags)
    for author in authors:
        agent = entity.bundle.agent(
            namespace + ':' + author['name'],
            {'attribute:' + k: author[k]
             for k in author if k != 'name'})
        entity.wasAttributedTo(agent)


def attribute_to_projects(entity, project_tags):
    """Attribute entity to projecs."""
    namespace = 'project'
    create_namespace(entity.bundle, namespace)

    projects = replace_tags('projects', project_tags)
    for project in projects:
        agent = entity.bundle.agent(namespace + ':' + project)
        entity.wasAttributedTo(agent)


def get_recipe_provenance(documentation, filename):
    """Create a provenance entity describing a recipe."""
    provenance = ProvDocument()

    for namespace in ('recipe', 'attribute'):
        create_namespace(provenance, namespace)

    references = replace_tags('references', documentation.get(
        'references', []))

    entity = provenance.entity(
        'recipe:{}'.format(filename), {
            'attribute:description': documentation.get('description', ''),
            'attribute:references': ', '.join(references),
        })

    attribute_to_authors(entity, documentation.get('authors', []))
    attribute_to_projects(entity, documentation.get('projects', []))

    return entity


def get_task_provenance(task, recipe_entity):
    """Create a provenance activity describing a task."""
    provenance = ProvDocument()
    create_namespace(provenance, 'task')
    # TODO: add this to task output product instead
    #     create_namespace(provenance, 'attribute')
    #     attributes = {}
    #     if hasattr(task, 'settings'):
    #         for attr in task.settings:
    #             attributes['attribute:' + attr] = str(task.settings[attr])
    #     if hasattr(task, 'script'):
    #         attributes['attribute:script'] = task.script

    activity = provenance.activity('task:' + task.name)

    trigger = recipe_entity
    update_without_duplicating(provenance, recipe_entity.bundle)

    starter = ESMVALTOOL_PROVENANCE
    update_without_duplicating(provenance, starter.bundle)

    activity.wasStartedBy(trigger, starter)

    return activity


class TrackedFile(object):
    """File with provenance tracking."""

    def __init__(self, filename, attributes, ancestors=None):

        self._filename = filename
        self.attributes = copy.deepcopy(attributes)

        self.provenance = None
        self.entity = None
        self._ancestors = [] if ancestors is None else ancestors
        self._activity = None

    @property
    def filename(self):
        """Filename."""
        return self._filename

    def __hash__(self):
        """Return the hash value of the object."""
        return hash(self.filename)

    def initialize_provenance(self, activity):
        """Initialize the provenance document."""
        if self.provenance is not None:
            raise ValueError(
                "Provenance of {} already initialized".format(self))
        self.provenance = ProvDocument()
        self._initialize_namespaces()
        self._initialize_activity(activity)
        self._initialize_entity()
        self._initialize_ancestors(activity)

    def _initialize_namespaces(self):
        """Inialize the namespaces."""
        for namespace in ('file', 'attribute', 'preprocessor', 'task'):
            create_namespace(self.provenance, namespace)

    def _initialize_activity(self, activity):
        """Copy the preprocessor task activity."""
        self._activity = activity
        update_without_duplicating(self.provenance, activity.bundle)

    def _initialize_entity(self):
        """Initialize the entity representing the file."""
        attributes = {
            'attribute:' + k: str(v)
            for k, v in self.attributes.items()
        }
        self.entity = self.provenance.entity('file:' + self.filename,
                                             attributes)

    def _initialize_ancestors(self, activity):
        """Register ancestor files for provenance tracking."""
        for ancestor in self._ancestors:
            if ancestor.provenance is None:
                ancestor.initialize_provenance(activity)
            update_without_duplicating(self.provenance, ancestor.provenance)
            self.wasderivedfrom(ancestor)

    def wasderivedfrom(self, other):
        """Let the file know that it was derived from other."""
        if isinstance(other, TrackedFile):
            update_without_duplicating(self.provenance, other.provenance)
            entity = other.entity
        else:
            entity = other
        if not self._activity:
            raise ValueError("Activity not initialized.")
        self.entity.wasDerivedFrom(entity, self._activity)

    def _select_for_include(self):
        attributes = {
            'provenance': self.provenance.serialize(format='xml'),
        }
        if 'caption' in self.attributes:
            attributes['caption'] = self.attributes['caption']
        return attributes

    @staticmethod
    def _include_provenance_nc(filename, attributes):
        with Dataset(filename, 'a') as dataset:
            for key, value in attributes.items():
                setattr(dataset, key, value)

    @staticmethod
    def _include_provenance_png(filename, attributes):
        pnginfo = PngInfo()
        exif_tags = {
            'provenance': 'ImageHistory',
            'caption': 'ImageDescription',
        }
        for key, value in attributes.items():
            pnginfo.add_text(exif_tags.get(key, key), value, zip=True)
        with Image.open(filename) as image:
            image.save(filename, pnginfo=pnginfo)

    def _include_provenance(self):
        """Include provenance information as metadata."""
        attributes = self._select_for_include()

        # List of files to attach provenance to
        files = [self.filename]
        if 'plot_file' in self.attributes:
            files.append(self.attributes['plot_file'])

        # Attach provenance to supported file types
        for filename in files:
            ext = os.path.splitext(filename)[1].lstrip('.').lower()
            write = getattr(self, '_include_provenance_' + ext, None)
            if write:
                write(filename, attributes)

    def save_provenance(self):
        """Export provenance information."""
        self._include_provenance()
        filename = os.path.splitext(self.filename)[0] + '_provenance'
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
