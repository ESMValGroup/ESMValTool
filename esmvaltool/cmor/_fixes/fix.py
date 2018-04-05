"""Contains the base class for model fixes"""
import importlib
import os


class Fix(object):
    """
    Base class for model fixes.
    """

    def fix_file(self, filepath, output_dir):
        """
        Apply fixes to the files prior to creating the cube.

        Should be used only to fix errors that prevent loading or can
        not be fixed in the cube (i.e. those related with missing_value
        and _FillValue)

        Parameters
        ----------
        filepath: basestring
            file to fix
        output_dir: basestring
            path to the folder to store the fixe files, if required

        Returns
        -------
        basestring
            Path to the corrected file. It can be different from the original
            filepath if a fix has been applied, but if not it should be the
            original filepath

        """
        return filepath

    def fix_metadata(self, cube):
        """
        Apply fixes to the metadata of the cube.

        Changes applied here must not require data loading.

        These fixes should be applied before checking the metadata.

        Parameters
        ----------
        cube: iris.cube.Cube
            Cube to fix

        Returns
        -------
        iris.cube.Cube
            Fixed cube. It can be a difference instance.

        """
        return cube

    def fix_data(self, cube):
        """
        Apply fixes to the data of the cube.

        These fixes should be applied before checking the data.

        Parameters
        ----------
        cube: iris.cube.Cube
            Cube to fix

        Returns
        -------
        iris.cube.Cube
            Fixed cube. It can be a difference instance.

        """
        return cube

    def __eq__(self, other):
        return type(self) == type(other)

    def __ne__(self, other):
        return not (self == other)

    @staticmethod
    def get_fixes(project, model, variable):
        """
        Get the fixes that must be applied for a given dataset.

        It will look for them at the module
        esmvaltool.cmor._fixes.PROJECT in the file MODEL, and get
        the classes named allvars (which should be use for fixes that are
        present in all the variables of a model, i.e. bad name for the time
        coordinate) and VARIABLE (which should be use for fixes for the
        specific variable).

        Project, model and variable names will have '-' replaced by '_' before
        checking because it is not possible to use the character '-' in python
        names.

        Parameters
        ----------
        project: str
        model: str
        variable: str

        Returns
        -------
        list(Fix)
            Fixes to apply for the given data
        """
        project = project.replace('-', '_')
        model = model.replace('-', '_')
        variable = variable.replace('-', '_')

        fixes = []
        try:
            fixes_module = importlib.import_module(
                'esmvaltool.cmor._fixes.{0}.{1}'.format(
                    project, model))
            for fix_name in ('allvars', variable):
                try:
                    fixes.append(getattr(fixes_module, fix_name)())
                except AttributeError:
                    pass
        except ImportError:
            pass
        return fixes

    @staticmethod
    def get_fixed_filepath(output_dir, filepath):
        """
        Get the filepath for the fixed file

        Parameters
        ----------
        var_path: str
            Original path

        Returns
        -------
        str
            Path to the fixed file
        """
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        return os.path.join(output_dir, os.path.basename(filepath))
