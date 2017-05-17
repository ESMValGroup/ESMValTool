# Online ESGF search functions for ESMValTool
#
# 2015-11-10  SR
# 2015-11-26  SR - Switch to dummy X509 cert

from os import path as os_path
from os import environ as os_environ
from sys import path as sys_path
from getpass import getpass
import urllib2

from pyesgf.search import SearchConnection as ESGFSearchConnection 


class ESGFSearchException(Exception):
    pass


class ESGFSearch:

    def __init__(self, esgf_config, info):
        """
        Initiates search class
        :param esgf_config: Instance of esgf_config.ESGFConfig class
        :param info: Message output function
        """
        self.config = esgf_config
        self.info = info

    def search(self, distrib=True, model_str='', **constraints):
        """
        Performs facet-based ESGF search using ESGF pyclient
        :param distrib: False = search local node, True = search all nodes
        :param model_str: ESMValTool model string (for text output only)
        :param constraints: Facet names and values, given as kwargs
        """
        # Apply facet name fixes
        for con in constraints:
            if con == "mip":
                constraints["cmor_table"] = constraints.pop("mip")
            elif con == "time_freq":
                constraints["time_frequency"] = constraints.pop("time_freq")

        connection = ESGFSearchConnection(\
            self.config.search_service_url,
            distrib=distrib)
        datasets = self._get_datasets(connection, **constraints)
        num_matches = len(datasets)
        self.info("num_matches = %s" % num_matches)

        #TODO: Handle num_matches, test each facet individually
        #      then in sequence

        # If a unambiguous dataset match is found, return download URLs 
        # of files in dataset
        if num_matches == 1:

            files = datasets[0].file_context().search(**constraints)
            download_urls = []
            for f in files:
                download_urls.append(f.download_url)

            return self._generate_download_instructions(
                download_urls,
                model_str)

        # If no matches, try to determine which facets are causing 
        # the problem 
        elif num_matches == 0:

            num_matches_by_facet = {}

            # Loop through the constraints, getting hit count for each facet
            for facet_name in constraints:
                facet_value = constraints[facet_name]
                single_constraint = {facet_name : facet_value}
                num_matches_by_facet[facet_name]\
                    = self._get_num_matches(connection, **single_constraint)

            # Return facet report
            return self._generate_zero_match_report(
                num_matches_by_facet,
                model_str,
                constraints)

        # If multiple matches, return error report
        else:
            return self._generate_multi_match_report(
                datasets,
                model_str,
                constraints)

    @staticmethod
    def _get_num_matches(connection, **constraints):
        """
        Get number of matching ESGF datasets for given constraint
        :param connection: esgf-pyclient connection object
        :param constraints: Facet names and values, given as kwargs
        """
        context = connection.new_context()
        context = context.constrain(**constraints)
        return context.hit_count

    @staticmethod
    def _get_datasets(connection, **constraints):
        """
        Get matching ESGF datasets for given constraint(s)
        :param connection: esgf-pyclient connection object
        :param constraints: Facet names and values, given as kwargs
        """
        context = connection.new_context()
        context = context.constrain(**constraints)
        return context.search()

    @staticmethod
    def _generate_download_instructions(download_urls,
                                        model_str):
        """
        Produce download instructions for user in ascii format
        """
        result = "To obtain this dataset, download the following files:"
        for url in download_urls:
            result += "\n%s" % url
        return result

    @staticmethod
    def _generate_zero_match_report(num_matches_by_facet,
                                    model_str,
                                    constraints):
        """
        Produce information about why zero matches where obtained
        for model facets (in ascii format)
        """
        at_least_one_facet_value_has_no_ESGF_match = False
        result = "No remote ESGF dataset was found to match the " +\
                 "following ESGF enabled model: \n\n" +\
                 "<model> %s </model>" % model_str +\
                 "\n\nFor your information, the number of ESGF datasets " +\
                 "matching each individual model specifier is as follows:\n\n"

        # Loop through each facet, reporting number of ESGF datasets
        # matching this facet alone
        for facet_name in num_matches_by_facet:
            facet_value = constraints[facet_name]
            num_matches = num_matches_by_facet[facet_name]
            result += "%s = '%s' matches %i datasets\n"\
                %(facet_name, facet_value, num_matches)
            if num_matches == 0:
                at_least_one_facet_value_has_no_ESGF_match = True

        # Add advice line if any facet has zero matches
        if at_least_one_facet_value_has_no_ESGF_match:
            result += "\nPay particular attention to any model " +\
                      "specifier above that matches zero " +\
                      "datasets, as these are sufficient in " +\
                      "themselves to prevent a match for the model.\n\n"

        return result

    @staticmethod
    def _generate_multi_match_report(datasets,
                                     model_str,
                                     constraints):
        """
        Produce information about ambiguous dataset match
        (in ascii format)

        TODO: Identify which dataset is original and which are replicas"
        """

        result = "Multiple remote ESGF datasets were found that match the " +\
                 "following ESGF enabled model: \n\n" +\
                 "<model> %s </model>" % model_str +\
                 "\n\nFor your information, the matching ESGF datasets are:\n"

        # Loop through each dataset, reporting the dataset id and number of files
        result += "-----------"
        for match_num, dataset in enumerate(datasets):
            result += "\nMatch {0}".format(match_num+1) +\
                      "\n{0} files".format(dataset.number_of_files) +\
                      "\nID = {0}".format(dataset.dataset_id)
            result += "\nURLS ="
            files = dataset.file_context().search(**constraints)
            for url_num, f in enumerate(files):
                result += "\n{0}:{1}".format(url_num+1,f.download_url)
            result += "\n-----------"

        return result
