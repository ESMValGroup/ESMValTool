# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 15:43:30 2016

@author: bmueller
"""
import xml.etree.ElementTree as XMLT
import xml.dom.minidom as XMLPP
import netCDF4 as nc
import os


class METAdata(object):

    def __init__(self, dtype="xml", modfile=None, data_dict={}, **kwargs):
        super(METAdata, self).__init__(**kwargs)
        """
        Default values
        """
        self.__avail__ = ["xml", "meta", "both"]
        # http://www.sno.phy.queensu.ca/~phil/exiftool/exiftool_pod.html;
        # relevant r/w formats
        #self.__meta_formats__ = ["jpg", "jpeg", "png", "eps", "mp4",
        #                         "tiff", "pdf", "ps"]
        self.__meta_formats__ = ["png"]
        self.__exif_maintag__ = "Exif.Image.ImageDescription"
        # self.__dtype__=dtype
        self.set_type(dtype)
        self.__data_dict__ = data_dict
        self.__modfile__ = modfile

        self.tags = None
        self.DataIDs = None

    def set_type(self, dtype):
        """
        set data type of metadata output
        currently: 'xml', 'meta', 'both'
        the routine also checks if the EXIV library is available
        if not, then only XML is written
        """
        if dtype in self.__avail__:
            self.__dtype__ = dtype
        else:
            assert False, "This type (" + dtype + ") is not defineded " + \
                "as available data type. Please consider writing [" + \
                ", ".join(self.avail) + "]!"
        # check if exif is available
        try:
            from PIL import Image, PngImagePlugin
        except:
            self.__dtype__ = 'xml'

    def get_avail(self):
        return self.__avail__

    def get_file(self):
        return self.__modfile__

    def get_dict(self):
        return self.__data_dict__

    def get_type(self):
        return self.__dtype__

    def set_dict(self, dictionary):
        if isinstance(dictionary, dict):
            self.__data_dict__ = dictionary
        else:
            assert False, "Input is not a dictionary!"

        if "DataIDs" in self.__data_dict__['ESMValTool'].keys():
            self.__get_DataIDs__()

    def set_file(self, modfile):
        if isinstance(modfile, (unicode, str)):
            self.__modfile__ = modfile
        else:
            print(type(modfile))
            assert False, "Input is not a string!"

    def write(self):
        if self.__dtype__ == "xml":
            self.__adjust_xml_file__()
            self.__write_xml__()
        elif self.__dtype__ == "meta":
            self.__write_meta__()
        elif self.__dtype__ == "both":
            self.__write_meta__()
            self.__adjust_xml_file__()
            self.__write_xml__()
        else:
            assert False, "Wrong data type! This should not happen!"

    def __write_meta__(self):

        from PIL import Image, PngImagePlugin

        if not self.__modfile__.split(".")[-1] in self.__meta_formats__:
            print("Warning! ." + self.__modfile__.split(".")[-1] +
                  " is not an acceptable " +
                  "meta data file format! Instead, XML-file for meta data " +
                  "will be produced!")
            self.__dtype__ = "xml"
            self.write()
            return

        if len(self.__data_dict__.keys()) == 1:

            meta = PngImagePlugin.PngInfo()

            image = Image.open(self.__modfile__)

            tags = self.__build_pretty_xml__(pretty=False)
            self.tags = self.__build_pretty_xml__(pretty=True)

            meta.add_itxt(self.__exif_maintag__, tags, zip=True)

            image.save(self.__modfile__, "png", pnginfo=meta)

        elif len(self.__data_dict__.keys()) == 0:
            print("No meta data to be written!")

        else:
            assert False, "The meta data dictionary" + \
                  " has too many entries on level 0!"

        return

    def __write_xml__(self):

        if len(self.__data_dict__.keys()) == 1:

            self.tags = self.__build_pretty_xml__()

            with open(self.__modfile__, "w") as XMLfile:
                XMLfile.write(self.tags)
                print("Meta data written to " + self.__modfile__ + "!")

        elif len(self.__data_dict__.keys()) == 0:
            print("No meta data to be written!")
            open(self.__modfile__, "w")

        else:
            assert False, "The meta data dictionary has too many entries" + \
                  " on level 0!"

        return

    def __adjust_xml_file__(self, modfile=None):

        if modfile is not None:
            self.__modfile__ = modfile

        file_path = self.__modfile__.split("/")
        file_elements = file_path[-1].split(".")

        if not file_elements[0] == "":
            file_elements = [""] + file_elements

        if not file_elements[-1] == "xml":
            file_elements = file_elements + ["xml"]

        file_path[-1] = ".".join(file_elements)

        self.set_file("/".join(file_path))
        return

    def __build_pretty_xml__(self, pretty=True):

        def prettify(elem):
            rough_string = XMLT.tostring(elem, 'utf-8')
            reparsed = XMLPP.parseString(rough_string)
            return reparsed.toprettyxml(indent="\t")

        def unprettify(elem):
            return XMLT.tostring(elem, 'utf-8')

        root = XMLT.Element(self.__data_dict__.keys()[0])
        self.__build_xml_tree__(root, self.__data_dict__[
                self.__data_dict__.keys()[0]])

        if pretty:
            root = prettify(root)
        else:
            root = unprettify(root)

        return root

    def __build_xml_tree__(self, curr_root, act_dict):

        for branch in act_dict.keys():
            if isinstance(act_dict[branch], dict):
                BRANCH = XMLT.SubElement(curr_root, branch)
                self.__build_xml_tree__(BRANCH, act_dict[branch])
            else:
                if isinstance(act_dict[branch], (int, str, float, unicode)):
                    XMLT.SubElement(curr_root, branch).text = act_dict[branch]
                elif isinstance(act_dict[branch], list):
                    XMLT.SubElement(curr_root, branch).text = "|".join(
                        act_dict[branch])
                elif act_dict[branch] is None:
                    XMLT.SubElement(curr_root, branch).text = ""
                else:
                    print(act_dict[branch])
                    assert False, "This data type is not implemented yet: " + \
                        str(type(act_dict[branch]))
        return

    def read(self, modfile=None):

        if modfile is None:
            modfile = self.__modfile__

        self.__dtype__ = modfile.split(".")[-1]

        if self.__dtype__ not in self.__meta_formats__:
            self.__dtype__ = "xml"
            if modfile.split("/")[-1].split(".")[0] != "":
                self.__adjust_xml_file__(modfile)
        else:
            self.__dtype__ = "meta"
            self.__modfile__ = modfile

        if self.__dtype__ == "xml":
            self.__read_xml__()
        elif self.__dtype__ == "meta":
            try:
                self.__read_meta__()
            except:
                self.__adjust_xml_file__()
                self.__read_xml__()
        else:
            assert False, "Wrong data type! This should not happen!"

        return self

    def __read_xml__(self):

        tree = XMLT.parse(self.__modfile__)
        root = tree.getroot()
        if len(root.getchildren()) == 0:
            assert False, "Object " + root.tag + " is misformed!"
        val = {}
        [val.update(self.__xml_to_dict__(branch))
         for branch in root.getchildren()]
        self.__data_dict__ = {root.tag: val}
        return self.__data_dict__

    def __read_meta__(self):

        from PIL import Image, PngImagePlugin
        from PIL.ExifTags import TAGS

        image = Image.open(self.__modfile__)
        info = image.info
        if self.__exif_maintag__ in info.keys():
            root = XMLT.fromstring(info[self.__exif_maintag__])
            if len(root.getchildren()) == 0:
                assert False, "Object " + root.tag + " is misformed!"
            val = {}
            [val.update(self.__xml_to_dict__(branch))
             for branch in root.getchildren()]
            self.__data_dict__ = {root.tag: val}
            return self.__data_dict__
        else:
            return None

    def __xml_to_dict__(self, branch):

        if len(branch.getchildren()) == 0:
            val = branch.text.split("|")
        else:
            val = {}
            [val.update(self.__xml_to_dict__(subbranch))
             for subbranch in branch.getchildren()]

        return {branch.tag: val}

    def __get_DataIDs__(self):

        def __ID_hunter__(file_list):

            new_list = []
            for f in file_list:
                if os.path.isfile(f):
                    try:
                        DataInfo = nc.Dataset(f)
                        ti = "tracking_id"
                        if ti in DataInfo.__dict__.keys():
                            new_list.append(getattr(DataInfo, ti))

                        elif "infile_0000" in DataInfo.__dict__.keys():
                            infiles = [t for t in DataInfo.__dict__.keys() if
                                       t.split("_")[0] == "infile"]
                            infiles = [getattr(DataInfo, i_f) for
                                       i_f in infiles]
                            new_list.extend(__ID_hunter__(infiles))

                        else:
                            assert False, \
                                "There is header information missing!" + \
                                "Either tracking_id or infile_**** needed."
                    except:
                        new_list.append("NO ID found in file " + f + "!")
                else:
                    new_list.append(f)

            return new_list

        self.__data_dict__['ESMValTool']['DataIDs'] = \
            __ID_hunter__(
                    self.__data_dict__['ESMValTool']['DataIDs'].split(","))

# USAGE

# MDE=METAdata()
#
# MDE.set_type("meta")
#
# MDE.set_file("TEST.png")
#
# MDE.set_dict({"ROOT1":{"TEST1":"Val1","TEST2":{"TEST3":"Val3","TEST4":"Val4"}}})
#
# MDE.write()
#
# MDX=METAdata()
#
# MDX.set_type("xml")
#
# MDX.set_file("test.sample.log")
#
# MDX.set_dict({"ROOT1":{"TEST1":"Val1","TEST2":{"TEST3":"Val3","TEST4":"Val4"}}})
#
# MDX.write()
#
# MDR=METAdata()
#
# MDR.set_file("test.sample.log")
#
# MDR.read()
#
# MDR.get_dict()
#
# MDR.read("TEST.png")
#
# MDR.get_dict()

# Standard Dict for ESMValTool:


# DICT={'ESMValTool':{
#    'built':'datetime',
#    'tags':['tag1','tag2','tag3'],
#    'caption':'CAPTIONTEXT',
#    'block':'#123'
#        }}
