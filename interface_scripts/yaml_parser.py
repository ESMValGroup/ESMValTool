import yaml

"""Quick parser for playing with the yaml syntax"""

test_nml = "../nml/namelist_template.yaml"

with open(test_nml, 'r') as stream:
    data_loaded = yaml.load(stream)

for ii in data_loaded["MODELS"]:
    print ii

print "---------------------------------------"

print data_loaded["PREPROCESS"]

print "---------------------------------------"

for diag in data_loaded["DIAGNOSTICS"]:
    for key, value in diag.items():
        for ii in value.items():
            print ii

