import pdb


class Folders(object):
    def __init__(self, folder_line):
        self.folder_line = folder_line

    def get_folder_line(self):
        return self.folder_line

    def split_entries(self):
        return self.folder_line.split()

    def __str__(self):
        folder_line = self.get_folder_line()
        return folder_line


class AllFolders:
    def __init__(self):
        self.folders = []

    def append(self, folder):
        self.folders.append(folder)

    def __str__(self):
        for folder in self.folders:
            print folder

    def __iter__(self):
        for folder in self.folders:
            yield folder
