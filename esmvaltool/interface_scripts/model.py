import pdb


class Model(object):
    def __init__(self, model_line, attributes, diag_specific_model):
        self.model_line = model_line
        self.diag_specific = diag_specific_model
        self.attributes = attributes

    def get_model_line(self):
        return self.model_line

    def get_model_attr_id(self):
        return [item.id for item in self.attributes]

    def is_diag_specific(self):
        return self.diag_specific

    def split_entries(self):
        return self.model_line.split()

    def __str__(self):
        model_line = self.get_model_line()
        return model_line


class AllModels:
    def __init__(self):
        self.models = []

    def append(self, model):
        self.models.append(model)

    def __str__(self):
        for model in self.models:
            print model.get_model_line()

    def __iter__(self):
        for model in self.models:
            yield model
