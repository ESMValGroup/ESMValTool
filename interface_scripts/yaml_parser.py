import yaml

"""Quick parser for playing with the yaml syntax"""

test_nml = "../nml/namelist_template.yaml"

with open(test_nml, 'r') as stream:
    data_loaded = yaml.load(stream)


print(data_loaded["GLOBAL"])

print(data_loaded["MODELS"])
