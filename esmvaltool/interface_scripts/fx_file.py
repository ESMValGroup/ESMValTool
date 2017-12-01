# Class to hold information about fx-files
#
# 2014-10-23  SR
# 2015-05-04  ME - rewrote to handle arbitrary fx-file (hopefully)
#
# Currently holds just an unique ID value (string)
# and the full path of the file (including name and extension)
#
from __future__ import print_function


class FX_file_exception(Exception):
    pass


class FX_file(object):
    def __init__(self, identifier, fullpath):
        self.id = str(identifier)  # Note ID forced to string
        self.fullpath = fullpath

    def get_id(self):
        return self.id

    def get_fullpath(self):
        return self.fullpath

    def __str__(self):
        return "fx file '" + self.id + "' <" + self.fullpath + ">"


class AllFXfiles:
    def __init__(self):
        self.fx_files = {}

    def append(self, fx_file):
        """ @brief Add fx_file, or raise exception if fx_file ID not unique
            @param FX_file object
        """
        identifier = fx_file.get_id()
        if identifier in self.fx_files:
            msg = "a fx_file with ID '" + identifier +\
                  "' already exists in this namelist. " +\
                  "Each fx_file ID must be unique."
            raise FX_file_exception(msg)
        else:
            self.fx_files[identifier] = fx_file

    def __str__(self):
        """ @brief print all fx_files contained within AllFXfiles object

            Note: run this as .__str__(), not from within a print statement
        """
        for identifier in self.fx_files:
            print(str(self.fx_files[identifier]))

    def __iter__(self):
        for identifier in self.fx_files:
            yield identifier

    def __getitem__(self, identifier):
        return self.fx_files[identifier]
