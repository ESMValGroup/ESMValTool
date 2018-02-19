yfile="""
--- !CMIP6Dataset
project: PROJECT
name: NAME
product: PRODUCT
institute: INSTITUTE
model: MODEL
experiment: EXPERIMENT
mip: MIP
ensemble: ENSEMBLE
start_year: START_YEAR
end_year: END_YEAR
variable: VARIABLE
"""
import yaml
class CMIP6Dataset(yaml.YAMLObject):
    yaml_tag = u'!CMIP6Dataset'
    def __init__(self, project, name, product, institute, model, experiment,
            mip, ensemble, start_year, end_year, variable):
            self.project = project
            self.name = name
            self.product = product
            self.institute = institute
            self.model = model
            self.experiment = experiment
            self.mip = mip
            self.ensemble = ensemble
            self.start_year = start_year
            self.end_year = end_year
            self.variable = variable
    def __repr__(self):
        return "%s( project=%r, name=%r, product=%r, institute=%r, model=%r, experiment=%r, mip=%r, ensemble=%r, start_year=%r, end_year=%r, variable=%r )" % (
            self.__class__.__name__, self.project , self.name ,
            self.product , self.institute , self.model , self.experiment ,
            self.mip , self.ensemble , self.start_year , self.end_year ,
            self.variable )


h = yaml.load_all(yfile)
for d in h:
    print(d)
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../nml-utils/generateNML'))

from generateNML import get_namelist

print(get_namelist())
