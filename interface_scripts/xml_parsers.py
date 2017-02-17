import pdb
import xml.sax
import nml_tags
import xml.dom.minidom
from auxiliary import info, error
import base64
import re
import os

from esgf_config import ESGFConfig
from esgf_config import ESGFTag

verbosity = 0

class configPrivateHandler(xml.sax.handler.ContentHandler):
	""" @brief SAX handler class for parsing the user defined configuration file

	"""
        def __init__(self):
            self.pathCollection = {}
	    self.ctag = []
   	    self.ccont = {}
            self.current_content = ""

        def startElement(self, name, attrs):
                self.current_content = ""
		self.ctag.append(name)
                self.ccont.update(attrs)

        def characters(self, content):
        	self.current_content += content.strip()

        def endElement(self, name):
		if len(self.ctag) == 4:
			if "/".join(self.ctag[-4:-1]) == "settings/pathCollection/usrpath":
	   	    		self.ccont.update({ name : self.current_content })
		if name == "usrpath":
			if 'id' in self.ccont.keys():
                                pathid = self.ccont.pop('id')
			else:
				pathid = hash(frozenset(self.ccont.items()))
				pathid = base64.b64encode(str(pathid))[:8]
			self.pathCollection[pathid] = self.ccont
			self.ccont = {}
		self.ctag.pop()

	def getPathByID(self,pathid):
		if pathid in self.pathCollection.keys():
			return self.pathCollection[pathid]["path"]
		else:
			return None


class ESGFConfigHandler(xml.sax.handler.ContentHandler):
    """ SAX handler class for parsing ESMValTool ESGF config files
    """
    def __init__(self):
        """ @brief Initialize SAX namelist handler variables
        """
        self.esgf_config_temp = {}
        self.string = ""
        self.attributes = None
        self.id = None
        self.elements = {}
        self.node = None

    def startElement(self, element_name, attr):
        """ default SAX startElement event handler
            :param name: default SAX startElement argument
            :param attr: default SAX startElement attribute
        """
        if element_name == "ESGF":
            self.current_tag = ESGFTag()

        elif element_name in ESGFConfig.valid_node_names:
            self.node = element_name

        elif element_name == ESGFConfig.user_cache_name:
            self.node = ESGFConfig.user_cache_name

        else:
            self.attributes = {}
            self.string = ""
            if attr:
                for attr_key in attr.keys():
                    self.attributes[attr_key.lower()] = attr[attr_key]

    def characters(self, data):
        self.string = self.string + data

    def endElement(self, element_name):
        """ @brief default SAX endElement event handler
            @param name default SAX endElement argument
        """
        info("endElement: " + element_name, verbosity, required_verbosity=12)

        if element_name == 'ESGF':
            pass

        elif element_name in ESGFConfig.valid_node_names\
            or element_name == "USER_CACHE":
            self.node = None

        else:
            self.current_tag.add_ESGF_entry(
                element_name,
                self.string.strip(),
                self.attributes,
                self.node
            )

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
	self.include = None
	self.conf = None

    def startElement(self, name, attr):
        """ @brief default SAX startElement event handler
            @param name default SAX startElement argument
            @param attr default SAX startElement attribute
        """
        info("startElement: " + name, verbosity, required_verbosity=12)
	if name == 'include':
		self.include = attr['href']
		self.conf = configPrivateHandler()
     		p = xml.sax.make_parser()
     		p.setContentHandler(self.conf)
     		p.parse(self.include)

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
	def replaceAttag(tbr):
	   pat = re.compile('@\\{..*\\}')
           if pat.search(tbr):
		patlist = pat.findall(tbr)
		patlist = [item[2:-1] for item in patlist]
		tmp = list(set(patlist) & set(self.conf.pathCollection.keys()))
		for t in tmp:
		   tbr = re.sub("@{"+t+"}", self.conf.getPathByID(t), tbr)
	   return tbr


        info("endElement: " + name, verbosity, required_verbosity=12)

	if self.include is not None:
           if name == self.current_tag.__class__.__name__:
               self.project_info[name] = getattr(self.current_tag, 'closing_tag')(replaceAttag(self.str), self.attributes)
           elif name == 'include':
               pass
           elif name == 'namelist':
               pass
           else:
               self.current_tag.add_nml_entry(name, replaceAttag(self.str), self.attributes)
	else:
           if name == self.current_tag.__class__.__name__:
               self.project_info[name] = getattr(self.current_tag, 'closing_tag')(self.str, self.attributes)
           elif name == 'include':
               pass
           elif name == 'namelist':
               pass
           else:
               self.current_tag.add_nml_entry(name, self.str, self.attributes)


    def endDocument(self):
	if self.include is not None:
		for k,v in self.conf.pathCollection.iteritems():
			if v['type'] not in  'output':
				if not os.path.exists(v['path']):
					error('Path {!s} defined in {!s} does not exist'.format(v['path'],
					      self.include))
				#	if k in self.project_info.keys():
				#		error('Path-ID {!s} in {!s} already used in different context. Please, choose  \
				#		      a different Path-ID'.format(k,self.include))
				#	else:
				#		self.project_info[k] = v['path']
				os.putenv("ESMValTool_" + k, str(v['path']))
