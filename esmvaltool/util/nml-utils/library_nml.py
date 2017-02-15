import pdb
import re
import xml.sax


class NamelistError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class splitDiags(xml.sax.handler.ContentHandler):
    """ @brief SAX handler class for parsing ESMValTool namelists

        SAX handler class for reading the ESMValTool XML-file and
        rewriting each diagnostic to separate xml-files.
    """
    def __init__(self):
        """ @brief Initialize SAX namelist handler variables
        """
        ## A dictionary holding the data read from the XML-file
        self.namelist_sections = {}
        self.namelist_sections["header"] = ""
        self.namelist_sections["diags"] = []
        self.namelist_sections["footer"] = ""
        self.str = ""
        self.attributes = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        if name == "namelist":
            self.current_location = "header"
            self.str = self.str + '<namelist>'
        elif name == "diag":
            if self.current_location == "header":
                self.namelist_sections["header"] = self.str
            self.str = '<diag>'
            self.current_location = "diag"
        else:
            self.str = self.str + '<' + name
            if attr:
                for key in attr.keys():
                    self.str = self.str + " " + key + '="' + attr[key] + '"'
            self.str = self.str + '>'

    def characters(self, data):
        self.str = self.str + data

    def endElement(self, name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        if name == "diag":
            self.str = self.str + '</diag>\n'
            self.namelist_sections["diags"].append(self.str)
            self.str = ""
        elif name == "namelist":
            self.current_location = "footer"
            self.namelist_sections["footer"] = self.str + '</namelist>'
            self.str = ""

        else:
            self.str = self.str + '</' + name + '>'


class addModels(xml.sax.handler.ContentHandler):
    """ @brief SAX handler class for parsing ESMValTool namelists

        SAX handler class for reading the ESMValTool XML-file and
        adding models
    """
    def __init__(self):
        """ @brief Initialize SAX namelist handler variables
        """
        ## A dictionary holding the data read from the XML-file
        self.namelist_sections = {}
        self.namelist_sections["header"] = ""
        self.namelist_sections["MODELS"] = ""
        self.namelist_sections["footer"] = ""
        self.str = ""
        self.attributes = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        if name == "namelist":
            self.current_location = "header"
            self.str = self.str + '<namelist>'
        elif name == "MODELS":
            if self.current_location == "header":
                self.namelist_sections["header"] = self.str
            self.str = ""
            self.current_location = "MODELS"
        else:
            self.str = self.str + '<' + name
            if attr:
                for key in attr.keys():
                    self.str = self.str + " " + key + '="' + attr[key] + '"'
            self.str = self.str + '>'

    def characters(self, data):
        self.str = self.str + data

    def endElement(self, name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        if name == "MODELS":
            self.namelist_sections["MODELS"] = self.str
            self.str = ""
        elif name == "namelist":
            self.current_location = "footer"
            self.namelist_sections["footer"] = self.str + '</namelist>'
            self.str = ""

        else:
            self.str = self.str + '</' + name + '>'


class setGlobal(xml.sax.handler.ContentHandler):
    """ @brief SAX handler class for parsing ESMValTool namelists

        SAX handler class for setting global settings the ESMValTool
        XML-file and adding models to the <MODELS>-tag
    """
    def __init__(self, dict_repl):
        """ @brief Initialize SAX namelist handler variables
        """
        ## A dictionary holding the data read from the XML-file
        self.dict_repl = dict_repl
        self.str = ""
        self.attributes = None
        self.current_tag = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        self.str = self.str + '<' + name
        if attr:
            for key in attr.keys():
                self.str = self.str + " " + key + '="' + attr[key] + '"'
        self.str = self.str + '>'
        self.current_tag = name

    def characters(self, data):
        self.str = self.str + data

    def endElement(self, name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        if name in self.dict_repl.keys():
            self.str = re.sub('(.*>).*$', r'\1' + self.dict_repl[name], self.str)
        self.str = self.str + '</' + name + '>'


class addModel(xml.sax.handler.ContentHandler):
    """ @brief SAX handler class for parsing ESMValTool namelists

        SAX handler class for reading the ESMValTool XML-file and
        adding model to the <diag>-tag
    """
    def __init__(self, variable):
        """ @brief Initialize SAX namelist handler variables
        """
        ## A dictionary holding the data read from the XML-file
        self.variable = variable
        self.namelist_sections = {}
        self.namelist_sections["header"] = ""
        self.namelist_sections["footer"] = ""
        self.str = ""
        self.attributes = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        if name == "namelist":
            self.current_location = "header"
            self.str = self.str + '<namelist>'
        else:
            self.str = self.str + '<' + name
            if attr:
                for key in attr.keys():
                    self.str = self.str + " " + key + '="' + attr[key] + '"'
            self.str = self.str + '>'

    def characters(self, data):
        self.str = self.str + data

    def endElement(self, name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        if name == "variable":
            self.str = self.str + '</' + name + '>'
            element_value = re.search('.*>\s+(.*)</variable>', self.str, re.DOTALL).group(1).strip()
            if self.variable == element_value:
                if len(self.namelist_sections["header"]) < 1:
                    self.namelist_sections["header"] = self.str
                    self.str = ""

        elif name == "namelist":
            self.current_location = "footer"
            self.namelist_sections["footer"] = self.str + '</namelist>'
            self.str = ""

        else:
            self.str = self.str + '</' + name + '>'


class setGlobal(xml.sax.handler.ContentHandler):
    """ @brief SAX handler class for parsing ESMValTool namelists

        SAX handler class for setting global settings the ESMValTool
        XML-file and adding models
    """
    def __init__(self, dict_repl):
        """ @brief Initialize SAX namelist handler variables
        """
        ## A dictionary holding the data read from the XML-file
        self.dict_repl = dict_repl
        self.str = ""
        self.attributes = None
        self.current_tag = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        self.str = self.str + '<' + name
        if attr:
            for key in attr.keys():
                self.str = self.str + " " + key + '="' + attr[key] + '"'
        self.str = self.str + '>'
        self.current_tag = name

    def characters(self, data):
        self.str = self.str + data

    def endElement(self, name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        if name in self.dict_repl.keys():
            self.str = re.sub('(.*>).*$', r'\1' + self.dict_repl[name], self.str)
        self.str = self.str + '</' + name + '>'
