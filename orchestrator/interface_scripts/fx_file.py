# Class to hold information about fx-files
#
# 2014-10-23  SR
# 2015-05-04  ME - rewrote to handle arbitrary fx-file (hopefully)
#
# Currently holds just an unique ID value (string)
# and the full path of the file (including name and extension)
#

import pdb


class FX_file_exception(Exception):
    pass


class FX_file(object):
    def __init__(self, Id, fullpath):
        self.Id = str(Id)  # Note ID forced to string
        self.fullpath = fullpath

    def get_Id(self):
        return self.Id

    def get_fullpath(self):
        return self.fullpath

    def __str__(self):
        return "fx file '" + self.Id + "' <" + self.fullpath + ">"


class AllFXfiles:
    def __init__(self):
        self.fx_files = {}

    def append(self, fx_file):
        """ @brief Add fx_file, or raise exception if fx_file ID not unique
            @param FX_file object
        """
        Id = fx_file.get_Id()
        if Id in self.fx_files:
            msg = "a fx_file with ID '" + Id +\
                  "' already exists in this namelist. " +\
                  "Each fx_file ID must be unique."
            raise FX_file_exception(msg)
        else:
            self.fx_files[Id] = fx_file

    def __str__(self):
        """ @brief print all fx_files contained within AllFXfiles object

            Note: run this as .__str__(), not from within a print statement
        """
        for Id in self.fx_files:
            print str(self.fx_files[Id])

    def __iter__(self):
        for Id in self.fx_files:
            yield Id

    def __getitem__(self, Id):
        return self.fx_files[Id]
