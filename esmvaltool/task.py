"""ESMValtool task definition"""
import pprint
import subprocess
import time


class AbstractTask(object):
    def __init__(self, settings, ancestors=None):

        self.settings = settings
        self.ancestors = [] if ancestors is None else ancestors
        self.output_data = None

    def run(self, input_data=None):

        if not self.output_data:
            if input_data is None:
                input_data = []
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

    def _write_settings(self):
        """Write settings to file

        In yaml format or a custom format interface file if that is preferred.
        """
        # TODO; implement or copy write_data_to_interface from
        # interface_script.data_interface

    def _run(self, input_data):
        """Run the diagnostic script."""
        # Run only preprocessor
        if self.script is None:
            return

        self.settings['input_files'] = input_data

        cmd = list(self.script)
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

    def __str__(self):
        def indent(txt):
            return '\n'.join('\t' + line for line in txt.split('\n'))

        txt = 'DiagnosticTask:\nscript: {}\nsettings:\n{}\nancestors:\n{}'.format(
            self.script,
            pprint.pformat(self.settings, indent=2),
            '\n\n'.join(indent(str(task)) for task in self.ancestors)
            if self.ancestors else 'None',
        )
        return txt
