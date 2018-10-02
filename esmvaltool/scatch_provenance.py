import datetime as dt

# visualize the graph
from prov.dot import prov_to_dot
from prov.model import ProvDocument


def create_sketch():
    # Create a new provenance document
    d1 = ProvDocument()
    # Declaring namespaces for various prefixes
    d1.add_namespace('evt', 'http://www.esmvaltool.org/scheme')

    # Entity:
    outfile = d1.entity(
        'evt:outfile.png', {
            'evt:caption': 'Caption of the Image',
            'evt:variable': 'Variable',
            'evt:plottype': 'Plottype',
            'evt:domain': 'Domain',
            'evt:theme': 'Theme',
            'evt:realm': 'Realm'
        })
    diagnostic = d1.entity(
        'evt:diagnostic', {
            'evt:description': 'Diag description',
            'evt:references': 'Reference',
            'evt:statistics': 'Statistics'
        })
    recipe = d1.entity(
        'evt:recipe', {
            'evt:description': 'Nml description',
            'evt:project': 'Project',
            'evt:references': 'Reference'
        })
    dataset1 = d1.entity('evt:dataset1')
    dataset2 = d1.entity('evt:dataset2')
    dataset3 = d1.entity('evt:dataset3')
    obs1 = d1.entity('evt:obs1')
    datasets = d1.entity('evt:datasets',
                         {'prov:type': 'prov:Collection'
                          })  # c is a collection, with unknown content
    d1.hadMember(datasets, dataset1)
    d1.hadMember(datasets, dataset2)
    d1.hadMember(datasets, dataset3)
    d1.hadMember(datasets, obs1)

    infile11 = d1.entity('evt:infile11.nc', {'evt:trackingID': 'TrackingID'})
    infile12 = d1.entity('evt:infile12.nc', {'evt:trackingID': 'TrackingID'})
    infile13 = d1.entity('evt:infile13.nc', {'evt:trackingID': 'TrackingID'})
    infile2 = d1.entity('evt:infile2.nc', {'evt:trackingID': 'TrackingID'})
    infile3 = d1.entity('evt:infile3.nc', {'evt:trackingID': 'TrackingID'})
    infile4 = d1.entity('evt:infile4.nc', {'evt:trackingID': 'TrackingID'})
    preprocfile = d1.entity('evt:preproc_file')
    software = d1.entity('evt:software', {
        'evt:ESMValTool': 'v2.0.0a',
        'evt:Python': '3.5',
        'evt:NCL': '6.2.1'
    })
    preproc_set = d1.entity(
        'evt:preproc_setting', {
            'evt:derivation': 'Derivation',
            'evt:timesel': 'Timesel',
            'evt:cmor': 'CMOR_fixes',
            'evt:levelint': 'Level interpolation',
            'evt:regridding': 'Regridding',
            'evt:masking': 'Masking',
            'evt:multimeanstat': 'Multimean statistics'
        })
    diag_set = d1.entity('evt:diag_setting')

    # Agent:
    author_nml1 = d1.agent('evt:Author_nml1')
    author_nml2 = d1.agent('evt:Author_nml2')
    author_diag1 = d1.agent('evt:Author_diag1')
    author_diag2 = d1.agent('evt:Author_diag2')

    # Adding an activity
    diagrun = d1.activity('evt:diagrun')
    preprocrun = d1.activity('evt:preprocrun')

    d1.wasDerivedFrom(outfile, recipe)

    d1.wasGeneratedBy(outfile, diagrun, dt.datetime.now())
    d1.wasGeneratedBy(preprocfile, preprocrun, dt.datetime.now())

    d1.used(diagrun, diagnostic)
    d1.used(diagrun, preprocfile)
    d1.used(diagrun, diag_set)
    d1.used(diagrun, software)
    d1.used(preprocrun, datasets)
    d1.used(preprocrun, software)
    d1.used(preprocrun, preproc_set)
    d1.used(dataset1, infile11)
    d1.used(dataset1, infile12)
    d1.used(dataset1, infile13)
    d1.used(dataset2, infile2)
    d1.used(dataset3, infile3)
    d1.used(obs1, infile4)

    d1.wasAttributedTo(recipe, author_nml1)
    d1.wasAttributedTo(recipe, author_nml2)
    d1.wasAttributedTo(diagnostic, author_diag1)
    d1.wasAttributedTo(diagnostic, author_diag2)

    print(d1.get_provn())
    d1.serialize('article-prov.xml', format='xml')
    dot = prov_to_dot(d1)
    dot.write_png('article-prov.png')
    #
    # Or save to a PDF
    dot.write_pdf('article-prov.pdf')


if __name__ == '__main__':
    create_sketch()
