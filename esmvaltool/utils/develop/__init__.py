"""Interface to utility commands for develop command group."""

import sys
from pathlib import Path

from . import compare
from . import recipe_filler


class DevelopCommand():
    """Development utilities."""

    def compare(self,
                reference_dir,
                current_dir,
                verbose=False):
        """Compare a recipe run to a reference run.

        Returns True if the runs were identical, False otherwise.

        Parameters
        ----------
        reference_dir : str
            Results directory from reference run
        current_dir : str
            Results directory from run to be tested
        verbose : bool
            Produce verbose output
        """
        same = compare.compare(Path(reference_dir), Path(current_dir), verbose)

        sys.exit(int(not same))


    def fill_recipe(self, recipe, output_recipe="recipe_autofilled.yml", config_file=""):
        """
        Fill in a partial recipe with additional datasets.
        
        Tool to obtain a set of additional datasets when given a partial recipe.
        The blank recipe should contain, to the very least, a list of diagnostics
        each with their variable(s). Example of minimum settings:
        
        diagnostics:
          diagnostic:
            variables:
              ta:
                mip: Amon
                start_year: 1850
                end_year: 1900
        
        Note that the tool will exit if any of these minimum settings are missing!
        
        Key features:
        
        - you can add as many variable parameters as are needed; if not added, the
          tool will use the "*" wildcard and find all available combinations;
        - you can restrict the number of datasets to be looked for with the `dataset:`
          key for each variable, pass a list of datasets as value, e.g.
          `dataset: [MPI-ESM1-2-LR, MPI-ESM-LR]`;
        - you can specify a pair of experiments eg `exp: [rcp26, rcp85]`
          for each variable; this will look for each available dataset per experiment
          and assemble an aggregated data stretch from each experiment; equivalent to
          esmvaltool's syntax of multiple experiments; this option needs an ensemble
          to be declared explicitly; it will return no entry if there are gaps in data
        - `start_year` and `end_year` are mandatory and are used to filter out the
          datasets that don't have data in the interval; if you want all possible years
          hence no filtering on years just use "*" for start and end years;
        - `config-user: rootpath: CMIPX` may be a list, rootpath lists are supported;
        
        Caveats:
        
        - the tool doesn't yet work for derived variables;
        - operation restricted to CMIP data.

        Parameters
        ----------

        recipe : str
            Path to partial recipe file
        output_recipe : str
            Path to output recipe
        config_file : str
            User configuration file
        """
        
        recipe_filler.run(recipe, output_recipe, config_file)
