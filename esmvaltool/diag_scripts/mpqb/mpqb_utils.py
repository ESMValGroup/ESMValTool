import yaml
import os

# Get metrics dict, get dataset color(s) or get datasetnames
def get_mpqb_cfg(cfgtype, cfgkey):
    """cfgtype: colormap, datasetcolor, datasetname"""
    # First read cfg file
    cfg_filename = os.path.join(os.path.split(__file__)[0], 'mpqb_cfg.yml')
    with open(cfg_filename, 'r') as handle:
        mpqb_cfg = yaml.safe_load(handle)
    # Defaults to specified defaults in yml file
    if cfgtype == 'colormap':
        if cfgkey in mpqb_cfg['colormaps']:
            return mpqb_cfg['colormaps'][cfgkey]
        else:
            return mpqb_cfg['colormaps']['default']
    # Defaults to alias (provided as cfgkey)
    if cfgtype == 'datasetname':
        if cfgkey in mpqb_cfg['datasetnames']:
            return mpqb_cfg['datasetnames'][cfgkey]
        else:
            return cfgkey
    # Defaults to black.
    if cfgtype == 'datasetcolor':
        if cfgkey in mpqb_cfg['datasetcolors']:
            return mpqb_cfg['datasetcolors'][cfgkey]
        else:
            return 'k'
