# Classes to hold config information for Earth System Grid Federation
#
# 2015-06-22  SR
# 2015-11-06  SR - Added user cache
# 2015-11-10  SR - Added ESGF search config
# 2015-11-26  SR - Switch to dummy X509 cert

import pdb
from re import findall as findall
from os import path as path

class ESGFConfigException(Exception):
    pass

# Note use of 'new-style' python classes in this file

class ESGFConfig(object):
    """
    Class to hold all singleton info items, plus all instances 
    of other types of ESGF config class 
    """
    # Valid ESGF node names (always in upper case)
    # All ESGF nodes included here for completeness, however
    # most won't be needed, as they don't have a file cache
    valid_node_names = [ 'ANL',
                         'BADC',
                         'BNU' 
                         'CMCC', 
                         'DKRZ',
                         'DKRZ_CMIP5',
                         'NOAA-ESRL',
                         'NOAA-GFDL',
                         'IPSL',
                         'NASA-GSFC',
                         'NASA-JPL',
                         'NCI',
                         'NERSC',
                         'ORNL',
                         'PCMDI' ]

    valid_nulls = ['NONE','N/A','NOT SPECIFIED']

    # Name of XML tag for user cache
    user_cache_name = "USER_CACHE"

    def __init__(self):
        self.local_node_name = None
        self.all_nodes = {} # This will include the user cache
        self.user_cache = None
        self.search_ESGF = False
        self.config_file_name = None # This is set manually
        self.user_openid = None
        self.search_service_url = None
        self.certif_service_url = None
        self.auth_realm = None
        self.X509_cert_file = None
        self.esgf_pyclient_dir = None
        # Other attributes can be added here as required

    def set_local_node(self, node_name):
        """
        Set current node (converts name to uppercase + checks valid) 
        :param node: Name of node, e.g. 'BADC' or 'DKRZ'
        """
        node_name_U = node_name.upper()
        if node_name_U in self.valid_nulls:
            # If node name in above list, set as None
            self.local_node_name = None
        elif not node_name_U in self.valid_node_names:
            # If node name invalid, raise error
            msg = "Node '%s' (converted to upper case " % node_name +\
                  "as '%s') is not valid. " % node_name_U +\
                  "Valid node names are %s. " % self.valid_nodes_names +\
                  "Valid null values for node name are %s" % self.valid_nulls 
            raise ESGFConfigException(msg)
        else:
            # If node name valid, set as current node
            self.local_node_name = node_name

    #def get_local_node_name(self): return self.local_node_name

    def get_local_node(self): 
        """
        Retrieve ESGF node config object for local node
        :params node: name of node
        :returns: ESGFNodeConfig instance
        """
        if self.local_node_name == None:
            return None
            """
            msg = 'No local node name specified in <ESGF> config section'
            raise ESGFConfigException(msg)
            """
        else:
            return self.get_node(self.local_node_name)

    #def get_search_ESGF(self): return self.search_ESGF

    def add_node(self, node_name):
        """
        Add new node (ignore if node exists already)
        :param node_name: Name of node
        :returns: Instance of ESGFNodeConfig belong to this node
        """
        if not node_name in self.all_nodes:
            self.all_nodes[node_name] = ESGFNodeConfig(node_name)
        return self.all_nodes[node_name]     

    def get_node(self, node_name):
        """
        Retrieve ESGF node config object for named node
        :params node_name: name of node
        :returns: ESGFNodeConfig instance
        """
        if node_name in self.all_nodes:
            return self.all_nodes[node_name]
        else:
            msg = "Node name '%s' not recognised."
            raise ESGFConfigException(msg)

    def add_user_cache(self):
        """
        Add user cache (ignore if exists already)
        :returns: Instance of ESGFNodeConfig representing user cache
        """
        if not self.user_cache:
            # Add user cache to list of nodes, this is so the ESGFTag 
            # class will process its <cache_root> and <cache_template> 
            # elements
            self.user_cache = self.add_node(self.user_cache_name)
        return self.user_cache

    def get_user_cache(self):
        """
        Retrieve ESGF node config object for user cache
        :returns: ESGFNodeConfig instance representing user cache
        """
        if self.user_cache:
            return self.user_cache
        else:
            return None

    def __str__(self):
        """
        Return string representation of ESGF config instance
        Note this calls __str__() method of ESGFNodeConfig class
        """
        # Print intro
        msg = '\nESGFConfig' +\
              '\n|' +\
              '\n|--local node: %s' % self.local_node_name +\
              '\n|' +\
              '\n|--search ESGF (if dataset not found locally): %s'\
              % self.search_ESGF +\
              '\n|' +\
              '\n|--node(s)' 
        # Print all nodes here except user cache
        for node_name in self.all_nodes:
            if not node_name == self.user_cache_name:
                msg += str(self.all_nodes[node_name])
        # Print user cache
        msg += '\n|' +\
               '\n|--user cache' +\
               str(self.user_cache)
        # Finish with horizontal rule
        msg += '\n' + "_" *54 
        return msg

class ESGFNodeConfig(object):
    """
    Class to hold information about a ESGF node,
    mainly the file cache
    """

    # Valid placeholder names in local node cache path template
    # These are all ESGF facets, except 'version'. Note 'project'
    # is NOT on this list as ESMValTool uses this to refer to the
    # project class. 'ESGF_project' should be used to refer to the
    # ESGF facet called 'project' 
    valid_placeholders = [ 'ESGF_project',
                           'product',
                           'institute',
                           'model',
                           'experiment',
                           'time_freq',
                           #'time_frequency',
                           'realm',
                           'mip',
                           #'cmor_table',
                           'ensemble',
                           'version',
                           'variable' ]

    # Version placeholder, used by get_version_path()
    version_ph = '[version]'

    def __init__(self, node_name):
        """
        :param node: Name of node that cache belongs to
        """
        self.node_name = node_name
        self.root = None
        self.path_templates = {}

    #def get_node_name(self): return self.node_name

    #def set_root(self, root): self.root = root

    #def get_root(self): return self.root

    def set_path_template(self, ptid, template, check_valid=True):
        """
        Set local node cache path template
        :param ptid: path template id (each node can have multiple templates)
        :param template: path template for all (or part) of local node cache
        :param check_valid: if True, checks validity of all placeholder names
        """
        # Raise error is path_template already exists for this node and id
        if ptid in self.path_templates:
            msg = "Cache path template with ptid '%s' " % ptid +\
                  "already exists for node '%s'." % self.node_name
            raise ESGFConfigException(msg)

        # If validity check selected, perform it now
        elif check_valid and not self._is_template_valid(template):
            msg = "node_cache_template '%s' " % template +\
                  'is invalid. The template should '\
                  'contain zero or more placeholders taken only '\
                  'from list %s, ' % self.valid_placeholders +\
                  'each surrounded by square brackets.'
            raise ESGFConfigException(msg)

        # If no check required, or check okay, set template
        else:    
            self.path_templates[ptid] = template

    def get_path_template(self, ptid):
        """
        Get path template
        :param ptid: path template id
        :returns: path template
        """
        return self.path_templates[ptid]

    """
    def get_all_path_templates(self):

        #Get all path templates
        #:returns: Dictionary containing path templates, referenced by id
 
        return self.path_templates
    """

    def get_dataset_path(self, ptid, **kwargs):
        """
        Returns the path to dataset in local node cache 
        :param ptid: id of path template
        :kwargs (optional): key/value pairs, where:
                key is name of placholder in template
                value to replaced placeholder  
        :returns: path to dataset (including cache root)
        """
        # Create a working copy of the template
        temp = self.path_templates[ptid]

        # Loop through kwargs, replacing placeholders one-by-one
        for key, value in kwargs.iteritems():
            temp = self._replace_placeholder(temp, key, value)

        # Prepend the root to the complete the path
        dataset_path = path.join(self.root,temp) 

        # Check if unfilled placeholders remain
        if self.has_placeholders(dataset_path):
            msg = "Dataset path '%s' contains " % dataset_path +\
                  'contains unfilled placeholders. ' +\
                  'Placeholders are text inside square brackets.'
            raise ESGFConfigException(mesg)

        # If not, return path
        else:
            return dataset_path

    def get_version_path(self, ptid, **kwargs):
        """
        Returns the path to version directory in local node cache 
        :param ptid: id of path template
        :kwargs (optional): key/value pairs, where:
                key is name of placholder in template
                value to replaced placeholder  
        :returns: path to version directory (including cache root)
                  or None if template has no version placeholder
        """
        # Create a working copy of the template
        temp = self.path_templates[ptid]

        # Check it contains a version placeholder
        if self.version_ph in temp:

            # Select substring to the left of the first 
            # '[version]' placeholder in the template
            temp = temp.split(self.version_ph)[0]

            # Loop through kwargs, replacing placeholders one-by-one
            for key, value in kwargs.iteritems():
                temp = self._replace_placeholder(temp, key, value)

            # Prepend the root to the complete the path
            version_path = path.join(self.root,temp) 

            # Check if unfilled placeholders remain
            if self.has_placeholders(version_path):
                msg = "Version path '%s' contains " % version_path +\
                      'contains unfilled placeholders. ' +\
                      'Placeholders are text inside square brackets.'
                raise ESGFConfigException(msg)

            # If all okay, return version path
            else: return version_path

        # If doesn't contain version placeholder, return None 
        else:
            msg = 'Path template %s contains no ' % temp +\
                  "version placeholder '%s'" % version
            raise ESGFConfigException(msg)

    @classmethod
    def _replace_placeholder(cls, template, key, value, strict=False):
        """
        Returns the template with one placeholder replaced
        :param template: A string containing the placeholder
        :param key: Name of placeholder to replace (not inc. '[' and ']')
        :param value: String with which placeholder is replaced
        :param strict: Raise exception if placeholder name not valid
        :returns: String holding template with placeholder replaced
        """
        if key not in cls.valid_placeholders:
            if strict:
                msg = "Placeholder '%s' is not valid. " % key\
                      + '\nValid names are %s' % cls.valid_placeholders
                raise ESGFConfigException(msg)
            else:
                # If strict = False, just return template unchanged
                return template
        else:
            placeholder = '[%s]' % key
            if placeholder in template:
                return template.replace(placeholder, value)
            else:
                if strict:
                    msg = 'Placeholder %s not present in template: %s'\
                          % (placeholder, template)
                    raise ESGFConfigException(msg)
                else:
                    # If strict = False, just return template unchanged
                    return template

    @classmethod
    def _is_template_valid(cls,template):
        """
        Checks if all placeholders in template are valid
        :param template: Template to check
        :returns: True all placeholders are valid, False otherwise
        """
        all_placeholder_names = findall(r'\[([^]]*)\]', template)
        for name in all_placeholder_names:
            if not name in cls.valid_placeholders:
                # Non-valid placeholder found
                return False
        # If we get to here, all placeholders are valid
        return True

    @staticmethod
    def has_placeholders(dataset_path):
        """
        Checks if all placeholders have been replaced with values
        :param template: template to be checked
        :returns: True if placeholders present, False otherwise
        A placeholder is zero of more characters in square brackets
        """
        all_placeholder_names = findall(r'\[([^]]*)\]', dataset_path)
        return len(all_placeholder_names)>0

    def __str__(self):
        """
        Prints a string representation of the instance,
        This is designed to be called from within ESGFConfig.__str__()
        """
        msg = '\n|  |' +\
              '\n|  |--%s' % self.node_name +\
              '\n|  |  |' +\
              '\n|  |  |--cache root: %s' % self.root +\
              '\n|  |  |' +\
              '\n|  |  |--cache template(s)' +\
              '\n|  |  |  |'
        for id in self.path_templates:
            # Get template (truncated to first 30 chars)
            template = self.path_templates[id][:30]
            msg += '\n|  |  |  |--%s: %s...' % (id, template)
        return msg


class ESGFTag(object):
    """
    Code to handle <ESGF> tag in XML ESGF config file
    and populate ESGFConfig() instance
    """
    def __init__(self):
        self.config = ESGFConfig()

    def add_ESGF_entry(self,
                       element_name,
                       element_string,
                       attributes,
                       node_name
                       ):

        # Handle elements that are not within a node element
        if element_name == 'local_node':
            self.config.set_local_node(element_string)

        elif element_name == 'search_ESGF':
            self.config.search_ESGF = _process_bool_element(
                element_name,
                element_string)

        # Handle online search configuration options
        #elif element_name == 'user_openid':
        #    self.config.user_openid = element_string

        elif element_name == 'search_service_url':
            self.config.search_service_url = element_string

        #elif element_name == 'certif_service_url':
        #    self.config.certif_service_url = element_string

        #elif element_name == 'auth_realm':
        #    self.config.auth_realm = element_string

        #elif element_name == 'X509_cert_file':
        #    self.config.X509_cert_file = element_string

        #elif element_name == 'esgf_pyclient_dir':
        #    self.config.esgf_pyclient_dir = element_string

        elif element_name == 'report_fullpath':
            self.config.report_fullpath = element_string

        # Otherwise, assume the element is node specfic
        else:

            # Check if node is local cache
            if node_name == self.config.user_cache_name:
                self.config.add_user_cache()

            # Otherwise, check node name is valid and add 
            # node to ESGF config. Note this does nothing 
            # if the node already exists, which is the case 
            # when processing a closing tag
            elif node_name in self.config.valid_node_names:
                self.config.add_node(node_name)

            # Raise exception if node name not recognised
            else:
                msg = "For element <%s> with value '%s', "\
                      % (element_name,element_string) +\
                      "node name '%s' is not recognised. " % node_name +\
                      'Is this element nested correctly within ' +\
                      'a valid node element such as <BADC> ... </BADC>?'
                raise ESGFConfigException(msg)         

            # Determine element name and process accordingly
            if element_name == 'cache_root':
                self.config.get_node(node_name).\
                    root = element_string

            if element_name == 'cache_template':
                if 'id' in attributes:
                    self.config.get_node(node_name).\
                        set_path_template(attributes['id'],element_string)
                else:
                    msg = "No 'id' attribute given for " +\
                          "<cache_template> of '%s' node " % node_name +\
                          "in ESGF config file."
                    raise ESGFConfigException(msg)


def _process_bool_element(element_name, element_string):
    """
    Processes element of ESGF config file and determines 
    if 'True' or 'False' specified, otherwise raises exception
    :param element_name: name of element (as given in opening tag)
    :param element_string: content of element
    :returns: True or False, accordingly
    """
    if element_string == 'True': return True
    elif element_string == 'False': return False
    else:
        msg = "Element <%s> in ESGF config file " % element_name +\
              "contains '%s'. " % element_string +\
              "Only 'True' or 'False' are acceptable."
        raise ESGFConfigException(msg)
