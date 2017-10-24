from xml.etree import ElementTree as et
datafile = "namelist_DiurnalCycle_box_SFCflux.xml"
tree = et.parse(datafile)
for m in tree.findall('MODELS/model'):
    print(m.text)
for m in tree.findall('DIAGNOSTICS/diag/model'):
    print(m.text)
