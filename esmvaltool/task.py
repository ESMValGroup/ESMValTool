"""ESMValtool task definition"""
import os
import subprocess
import time


class AbstractTask(object):
    def __init__(self, settings, ancestors=None):

        self.settings = settings
        self.ancestors = [] if ancestors is None else ancestors
        self.output_data = None

    def run(self, input_data=None):

        if input_data is None:
            input_data = []

        if not self.output_data:
            for task in self.ancestors:
                input_data.extend(task.run())
            self._run(input_data)

        return self.output_data

    def _run(self, input_data):
        raise NotImplementedError(
            "Method should be implemented by child class")


class DiagnosticError(Exception):
    """Error in diagnostic"""


class DiagnosticTask(AbstractTask):
    """Task for running a diagnostic"""
    
    def __init__(self, script, settings, ancestors=None):
        """Initialize"""
        super(DiagnosticTask, self).__init__(settings, ancestors)
        self.script = script
        self.executables = {
            'py': ['python'],
            'ncl': ['ncl'],
            'r': ['Rscript', '--slave', '--quiet'],
        }

    def _write_settings(self):
        """Write settings to file
        
        In yaml format or a custom format interface file if that is preferred.
        """
        # TODO; implement or copy write_data_to_interface from
        # interface_script.data_interface

    def _run(self, input_data):
        """Run the diagnostic script."""
        self.settings['input_files'] = input_data

        cmd = []        
        if os.path.isfile(self.script):
            extension = os.path.splitext(self.script).lower()[1:]
            cmd.append(self.executables[extension])
        cmd.append(self.script)

        settings_file = self._write_settings()
        cmd.append(settings_file)

        with open(self.stdout, 'a') as stdout, \
            open(self.stderr, 'a') as stderr:

            process = subprocess.Popen(
                cmd, stdout=stdout, stderr=stderr, env=self.settings['env'])
            returncode = None
            while returncode is None:
                # TODO: analyse stdout/stderr for ncl errors
                # TODO: log resource usage
                time.sleep(1)
                returncode = process.poll()

        if returncode in self.settings['allowed_returncodes']:
            return self.settings['output_files']

        raise DiagnosticError(
            "Diagnostic script {} failed with return code {}".format(
                self.script, returncode))
