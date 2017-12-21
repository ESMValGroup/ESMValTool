"""ESMValtool task definition"""
import logging
import os
import pprint
import subprocess
import time

from .interface_scripts.data_interface import write_settings

logger = logging.getLogger(__name__)


def which(executable):
    """Find executable in PATH."""
    for path in os.environ["PATH"].split(os.pathsep):
        if os.access(os.path.join(path, executable), os.X_OK):
            return os.path.join(path, executable)

    return None


class AbstractTask(object):
    """Base class for defining task classes"""

    def __init__(self, settings, ancestors=None):
        """Initialize task."""
        self.settings = settings
        self.ancestors = [] if ancestors is None else ancestors
        self.output_data = None

    def run(self, input_data=None):
        """Run task."""
        if not self.output_data:
            if input_data is None:
                input_data = []
            for task in self.ancestors:
                input_data.extend(task.run())
            self.output_data = self._run(input_data)

        return self.output_data

    def _run(self, input_data):
        raise NotImplementedError(
            "Method should be implemented by child class")

    def str(self):
        """Return a nicely formatted description."""

        def _indent(txt):
            return '\n'.join('\t' + line for line in txt.split('\n'))

        txt = 'settings:\n{}\nancestors:\n{}'.format(
            pprint.pformat(self.settings, indent=2),
            '\n\n'.join(_indent(str(task)) for task in self.ancestors)
            if self.ancestors else 'None',
        )
        return txt


class DiagnosticError(Exception):
    """Error in diagnostic"""


class DiagnosticTask(AbstractTask):
    """Task for running a diagnostic"""

    def __init__(self, script, settings, ancestors=None):
        """Initialize"""
        super(DiagnosticTask, self).__init__(settings, ancestors)
        self.script = script
        self.cmd = self._initialize_cmd(script)
        self.log = os.path.join(settings['output_dir'], 'log.txt')

    @staticmethod
    def _initialize_cmd(script):
        """Create a an executable command from script."""
        diagnostics_root = os.path.join(
            os.path.dirname(__file__), 'diag_scripts')
        script_file = os.path.abspath(os.path.join(diagnostics_root, script))

        if not os.path.isfile(script_file):
            raise DiagnosticError(
                "Cannot execute script {} ({}): file does not exist.".format(
                    script, script_file))

        cmd = []
        if not os.access(script_file, os.X_OK):  # if not executable
            extension = os.path.splitext(script)[1].lower()[1:]
            executables = {
                'py': [which('python')],
                'ncl': [which('ncl'), '-n'],
                'r': [which('Rscript'), '--slave', '--quiet'],
            }
            if extension not in executables:
                raise DiagnosticError(
                    "Cannot execute script {} ({}): non-executable file "
                    "with unknown extension.".format(script, script_file))

            cmd = executables[extension]

        cmd.append(script_file)

        return cmd

    def _write_settings(self):
        """Write settings to file

        In yaml format or a custom format interface file if that is preferred.
        """
        ext = 'yml'
        if self.script.lower().endswith('.ncl'):
            ext = 'ncl'
        filename = os.path.join(self.settings['output_dir'],
                                'settings.{}'.format(ext))
        write_settings(self.settings, filename)
        return filename

    def _control_ncl_execution(self, process, lines):
        """Check if an error has occurred in an NCL script.

        Apparently NCL does not automatically exit with a non-zero exit code
        if an error occurs, so we take care of that here.
        """
        ignore_warnings = [
            warning.strip()
            for warning in self.settings.get('ignore_ncl_warnings', [])
        ]

        errors = ['error:', 'fatal:']
        if self.settings['exit_on_ncl_warning']:
            errors.append('warning:')

        msg = ("An error occurred during execution of NCL script {}, "
               "see the log in {}".format(self.script, self.log))

        warned = False
        for line in lines:
            if line.strip() in ignore_warnings:
                continue
            if 'warning:' in line:
                warned = True
            for error in errors:
                if error in line:
                    logger.error(msg)
                    try:
                        process.kill()
                    except OSError:  # ignore error if process already exited
                        pass
                    else:
                        logger.error("Killed process.")
                    raise DiagnosticError(msg)

        if warned:
            logger.warning(
                "There were warnings during the execution of NCL script %s, "
                "for details, see the log %s", self.script, self.log)

    def _run(self, input_data):
        """Run the diagnostic script."""
        if self.script is None:  # Run only preprocessor
            output_data = []
            return output_data

        self.settings['input_files'] = input_data
        settings_file = self._write_settings()

        cmd = list(self.cmd)
        cwd = None
        env = {str(k): str(v) for k, v in self.settings['env'].items()}

        is_ncl_script = self.script.lower().endswith('.ncl')
        if is_ncl_script:
            cwd = os.path.dirname(__file__)
        else:
            cmd.append(settings_file)

        logger.info("Running command %s", cmd)
        logger.debug("in environment\n%s", pprint.pformat(env))
        logger.debug("in current working directory: %s", cwd)
        logger.info("Writing (some) output to %s", self.settings['output_dir'])
        logger.info("Writing log to %s", self.log)

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=cwd,
            env=env,
            universal_newlines=True)

        returncode = None
        last_line = ['']

        with open(self.log, 'at') as log:
            while returncode is None:
                returncode = process.poll()
                txt = process.stdout.read()
                log.write(txt)

                # Check if an error occurred in an NCL script
                # Last line is treated separately to avoid missing
                # error messages spread out over multiple lines.
                lines = txt.split('\n')
                if is_ncl_script:
                    self._control_ncl_execution(process, last_line + lines)
                last_line = lines[-1:]

                # TODO: log resource usage

                # wait, but not long because the stdout buffer may fill up:
                # https://docs.python.org/3.6/library/subprocess.html#subprocess.Popen.stdout
                time.sleep(0.001)

        if returncode == 0:
            return self.settings['output_dir']

        raise DiagnosticError(
            "Diagnostic script {} failed with return code {}. See the log "
            "in {}.".format(self.script, returncode))

    def __str__(self):
        """Get human readable description."""
        txt = "{}:\nscript: {}\n{}".format(
            self.__class__.__name__,
            self.script,
            super(DiagnosticTask, self).str(),
        )
        return txt
