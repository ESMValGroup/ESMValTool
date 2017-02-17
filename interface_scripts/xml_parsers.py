import pdb
import xml.sax
import nml_tags
import xml.dom.minidom
from auxiliary import info

verbosity = 0


class namelistHandler(xml.sax.handler.ContentHandler):
    """ @brief SAX handler class for parsing ESMValTool namelists

        SAX handler class for reading the ESMValTool XML-file and
        put its content into the following dictionary
        @code
        self.project_info = { 'MODELS':['first model string',
                                        'second model string', ...],
                              'GLOBAL':{var1:value1, var2:value2, ...},
                              'DIAGNOSTIC':['diagnostic1', 'diagnostic2', ...]
                              'AUXILIARIES':{'fx_files':{fx_id1:fx_file1,
                                                          fx_id2:fx_file2,
                                                          ...
                                                         }
                                            }

                            }
        @endcode
    """
    def __init__(self):
        """ @brief Initialize SAX namelist handler variables
        """
        ## A dictionary holding the data read from the XML-file
        self.project_info = {}
        self.str = ""
        self.attributes = None
        self.id = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        info("startElement: " + name, verbosity, required_verbosity=12)
        if name in vars(nml_tags):
            self.current_tag = vars(nml_tags)[name]()

        else:
            info("clearing str '" + self.str + "'", verbosity,
                 required_verbosity=12)
            self.attributes = {}
            self.str = ""
            if attr:
                for attr_key in attr.keys():
                    self.attributes[attr_key.lower()] = attr[attr_key]

    def characters(self, data):
        self.str = self.str + data

    def endElement(self, name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        info("endElement: " + name, verbosity, required_verbosity=12)

        if name == self.current_tag.__class__.__name__:
            self.project_info[name] = getattr(self.current_tag, 'closing_tag')(self.str, self.attributes)
        elif name == 'namelist':
            pass
        else:
            self.current_tag.add_nml_entry(name, self.str, self.attributes)
