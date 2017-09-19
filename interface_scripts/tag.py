import pdb


class Tags(object):
    def __init__(self, tag_line):
        self.tag_line = tag_line

    def get_tag_line(self):
        return self.tag_line

    def split_entries(self):
        return self.tag_line.split()

    def __str__(self):
        tag_line = self.get_tag_line()
        return tag_line


class AllTags:
    def __init__(self):
        self.tags = []

    def append(self, tag):
        self.tags.append(tag)

    def __str__(self):
        for tag in self.tags:
            print tag

    def __iter__(self):
        for tag in self.tags:
            yield tag
