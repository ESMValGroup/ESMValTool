import logging
import os

from prov.dot import prov_to_dot
from prov.model import ProvDocument

# from esmvaltool.preprocessor import MULTI_MODEL_FUNCTIONS, PreprocessingTask

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

    # def add_preprocessor_provenance(task, doc):
    #     for ancestor in task.ancestors:
    #         add_preprocessor_provenance(ancestor, doc)
    #
    #     if not isinstance(task, PreprocessingTask):
    #         return
    #
    #     task_activity = doc.activity('evt:' + task.name)
    #
    #     input_files = []
    #     if task._input_files:
    #         input_files.extend(task._input_files)
    #     for ancestor in task.ancestors:
    #         input_files.extend(f for f in ancestor.settings)
    #     grouped_input_files = _group_input(input_files, task.settings)
    #
    #     all_input = {}
    #     for output_file, input_files in grouped_input_files.items():
    #         output_entity = doc.entity('evt:' + output_file)
    #         input_collection = doc.collection('evt:input of ' + output_file)
    #         all_input[output_file] = input_collection
    #         for input_file in input_files:
    #             input_entity = doc.entity('evt:' + input_file)
    #             doc.hadMember(input_collection, input_entity)
    #
    #         # Add single dataset derivation relations
    #         settings = {
    #             'evt:' + k: str(v)
    #             for k, v in task.settings[output_file].items()
    #         }
    #         doc.wasDerivedFrom(
    #             output_entity,
    #             input_collection,
    #             task_activity,
    #             other_attributes=settings)
    #
    #     # Add multi dataset derivation relations
    #     steps = {j for i in task.settings.values() for j in i}
    #     for step in set(MULTI_MODEL_FUNCTIONS) - {'extract_metadata'}:
    #         if step in steps:
    #             multimodel_input_collection = doc.collection(
    #                 'evt:input of ' + task.name + ' ' + step)
    #             for output_file in all_input:
    #                 if step not in task.settings[output_file]:
    #                     continue
    #                 settings = task.settings[output_file][step]
    #                 if output_file in settings.get('exclude', {}).get(
    #                         '_filename', []):
    #                     continue
    #                 doc.hadMember(multimodel_input_collection,
    #                               all_input[output_file])
    #                 output_entities = []
    #                 if step == 'multi_model_statistics':
    #                     for filename in settings['filenames'].values():
    #                         output_entities.append(
    #                             doc.entity('evt:' + filename))
    #                 else:
    #                     output_entities.append(doc.entity(
    #                                                'evt:' + output_file))
    #                 settings = {'evt:' + step: str(settings)}
    #                 for output_entity in output_entities:
    #                     doc.wasDerivedFrom(
    #                         output_entity,
    #                         multimodel_input_collection,
    #                         task_activity,
    #                         other_attributes=settings)

    #     # Define the diagnostics
    #     for diagnostic in diagnostics:
    #         diagnostic_entity = doc.entity(
    #             'evt:diagnostic', {
    #                 'evt:description': diagnostic['description'],
    #                 'evt:references': 'Reference',
    #                 'evt:statistics': 'Statistics'
    #             })
    #         for script_name, script in diagnostics['scripts'].values():
    #             script_entity = doc.entity(
    #                 'evt:diagnostic_script', {
    #                     'evt:filename': script['script'],
    #                     'evt:settings': script['settings'],
    #                     'evt:ancestors': script['ancestors']
    #                 })
    #             script_run = doc.activity('evt:script_run')
    #
    #             doc.used(script_run, script_entity)
    #
    #     dataset1 = doc.entity('evt:dataset1')
    #     dataset2 = doc.entity('evt:dataset2')
    #     dataset3 = doc.entity('evt:dataset3')
    #     obs1 = doc.entity('evt:obs1')
    #     datasets = doc.entity('evt:datasets',
    #                           {'prov:type': 'prov:Collection'
    #                            })
    #     doc.hadMember(datasets, dataset1)
    #     doc.hadMember(datasets, dataset2)
    #     doc.hadMember(datasets, dataset3)
    #     doc.hadMember(datasets, obs1)
    #
    #     infile11 = doc.entity('evt:infile11.nc',
    #                           {'evt:trackingID': 'TrackingID'})
    #     infile12 = doc.entity('evt:infile12.nc',
    #                           {'evt:trackingID': 'TrackingID'})
    #     infile13 = doc.entity('evt:infile13.nc',
    #                           {'evt:trackingID': 'TrackingID'})
    #     infile2 = doc.entity('evt:infile2.nc',
    #                          {'evt:trackingID': 'TrackingID'})
    #     infile3 = doc.entity('evt:infile3.nc',
    #                          {'evt:trackingID': 'TrackingID'})
    #     infile4 = doc.entity('evt:infile4.nc',
    #                          {'evt:trackingID': 'TrackingID'})
    #     preprocfile = doc.entity('evt:preproc_file')
    #     software = doc.entity('evt:software', {
    #         'evt:ESMValTool': 'v2.0.0a',
    #         'evt:Python': '3.5',
    #         'evt:NCL': '6.2.1'
    #     })
    #     preproc_set = doc.entity(
    #         'evt:preproc_setting', {
    #             'evt:derivation': 'Derivation',
    #             'evt:timesel': 'Timesel',
    #             'evt:cmor': 'CMOR_fixes',
    #             'evt:levelint': 'Level interpolation',
    #             'evt:regridding': 'Regridding',
    #             'evt:masking': 'Masking',
    #             'evt:multimeanstat': 'Multimean statistics'
    #         })
    #     diag_set = doc.entity('evt:diag_setting')
    #
    #     diagrun = doc.activity('evt:diagrun')
    #     preprocrun = doc.activity('evt:preprocrun')
    #
    #     doc.wasGeneratedBy(preprocfile, preprocrun, datetime.now())
    #
    #     doc.used(preprocrun, datasets)
    #     doc.used(preprocrun, software)
    #     doc.used(preprocrun, preproc_set)
    #
    #     doc.used(diagrun, preprocfile)
    #     doc.used(diagrun, diag_set)
    #     doc.used(diagrun, software)
    #
    #     doc.used(dataset1, infile11)
    #     doc.used(dataset1, infile12)
    #     doc.used(dataset1, infile13)
    #     doc.used(dataset2, infile2)
    #     doc.used(dataset3, infile3)
    #     doc.used(obs1, infile4)

    return doc


def update_diagnostic_provenance(doc, output_dirs):
    """Update the provenance document with info provided by diagnostics."""
    outfile = doc.entity(
        'evt:outfile.png', {
            'evt:caption': 'Caption of the Image',
            'evt:variable': 'Variable',
            'evt:plottype': 'Plottype',
            'evt:domain': 'Domain',
            'evt:theme': 'Theme',
            'evt:realm': 'Realm'
        })


#     author_diag1 = doc.agent('evt:Author_diag1')
#     author_diag2 = doc.agent('evt:Author_diag2')
#     doc.wasAttributedTo(diagnostic_script, author_diag1)
#     doc.wasAttributedTo(diagnostic_script, author_diag2)
#     doc.wasDerivedFrom(outfile, recipe)
#     doc.wasGeneratedBy(outfile, diagrun, datetime.now())


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
