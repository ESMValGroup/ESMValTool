import importlib


class Fix(object):
    """
    Base class for fixes
    """
    def fix_file(self, filepath):
        return filepath

    def fix_metadata(self, cube):
        return cube

    def fix_data(self, cube):
        return cube

    def __eq__(self, other):
        return type(self) == type(other)

    def __ne__(self, other):
        return not (self == other)

    @staticmethod
    def get_fixes(project, model, variable):
        project = project.replace('-', '_')
        model = model.replace('-', '_')
        variable = variable.replace('-', '_')

        fixes = []
        try:
            fixes_module = importlib.import_module('orchestrator.interface_scripts.fixes.{0}.{1}'.format(project, model))
            for fix_name in ('allvars', variable):
                try:
                    fixes.append(getattr(fixes_module, fix_name)())
                except AttributeError:
                    pass
        except ImportError:
            pass
        return fixes
