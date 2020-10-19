"""Install Julia and R dependencies."""
import subprocess
import sys
from pathlib import Path


class Install:
    """Install extra dependencies.

    Diagnostics written in Julia or R need extra dependencies. Use this
    command to install them.

    Note that Julia or R must be pre-installed before running this command.
    """

    @staticmethod
    def _run(cmd, script):
        root = Path(__file__).parent
        try:
            subprocess.check_output(
                [cmd, str(root / script)],
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
        except subprocess.CalledProcessError as exc:
            print(exc.stdout)
            print("installation failed")
            sys.exit(1)
        else:
            print("Installation successful")

    def Julia(self):  # noqa: N802
        """Install dependencies needed to run Julia diagnostics."""
        print("installing Julia packages, please wait...")
        script = Path("Julia") / "setup.jl"
        self._run("julia", script)

    def R(self):  # noqa: N802
        """Install dependencies needed to run R diagnostics."""
        print("installing R packages, please wait...")
        print("Compiling may take up to 15 minutes or more.")
        script = Path("R") / "setup.R"
        self._run("Rscript", script)
