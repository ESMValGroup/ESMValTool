from esmvalcore.config import CFG
from esmvalcore.dataset import Dataset, datasets_to_recipe
import yaml
dataset_template = Dataset(
    short_name='tos',
    mip='Omon',
    project='CMIP6',
    dataset='*',
    institute='*',
    ensemble='*',
    exp='historical',
    grid='gn',
    timerange='1850/2000',
)
datasets = list(dataset_template.from_files())
len(datasets)
for dataset in datasets:
    dataset.facets['diagnostic'] = 'mca'
print(yaml.safe_dump(datasets_to_recipe(datasets)))