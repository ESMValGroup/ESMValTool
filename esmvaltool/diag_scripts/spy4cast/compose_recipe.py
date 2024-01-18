from esmvalcore.config import CFG
from esmvalcore.dataset import Dataset, datasets_to_recipe
import yaml
tos_template = Dataset(
    short_name='tos',
    mip='Omon',
    project='CMIP6',
    dataset='*',
    institute='*',
    ensemble='*',
    exp='historical',
    grid='*',
    timerange='*',
)
ps_template =  tos_template.copy(short_name='ps', mip='Amon')
datasets = list(tos_template.from_files()) + list(ps_template.from_files())
len(datasets)
for dataset in datasets:
    dataset.facets['diagnostic'] = 'mca'
print(yaml.safe_dump(datasets_to_recipe(datasets)))