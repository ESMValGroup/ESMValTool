from prov.constants import PROV_ATTR_GENERATED_ENTITY, PROV_ATTR_USED_ENTITY
from prov.model import ProvDerivation


def get_file_record(prov, filename):
    records = prov.get_record('file:' + filename)
    assert records
    return records[0]


def check_provenance(product):
    prov = product.provenance

    entity = get_file_record(prov, product.filename)
    assert entity == product.entity

    check_product_wasderivedfrom(product)


def check_product_wasderivedfrom(product):
    """Check that product.filename was derived from product._ancestors."""
    print('checking provenance of file', product.filename)
    prov = product.provenance

    def get_identifier(filename):
        record = get_file_record(prov, filename)
        return {record.identifier}

    # Check that the input and output file records exist
    identifier = get_identifier(product.filename)

    relations = {r for r in prov.records if isinstance(r, ProvDerivation)}
    for ancestor in product._ancestors:
        input_identifier = get_identifier(ancestor.filename)
        for record in relations:
            if input_identifier == record.get_attribute(PROV_ATTR_USED_ENTITY):
                assert identifier == record.get_attribute(
                    PROV_ATTR_GENERATED_ENTITY)
                break
        else:
            assert False

    if not product._ancestors:
        assert 'tracking_id' in product.attributes
    else:
        for ancestor in product._ancestors:
            check_product_wasderivedfrom(ancestor)
