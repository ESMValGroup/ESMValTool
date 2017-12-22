"""ESMValtool task definition"""
import logging
import os
import pprint
import subprocess
import time
from multiprocessing import Pool, cpu_count

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

    def flatten(self):
        """Return a flattened set of all ancestor tasks and task itself."""
        tasks = set()
        for task in self.ancestors:
            tasks.update(task.flatten())
        tasks.add(self)
        return tasks

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
            return [self.settings['output_dir']]

        raise DiagnosticError(
            "Diagnostic script {} failed with return code {}. See the log "
            "in {}.".format(self.script, returncode, self.log))

    def __str__(self):
        """Get human readable description."""
        txt = "{}:\nscript: {}\n{}".format(
            self.__class__.__name__,
            self.script,
            super(DiagnosticTask, self).str(),
        )
        return txt


def get_flattened_tasks(tasks):
    """Return a set of all tasks and their ancestors in `tasks`."""
    return set(t for task in tasks for t in task.flatten())


def get_independent_tasks(tasks):
    """Return a set of independent tasks."""
    independent_tasks = set()
    all_tasks = get_flattened_tasks(tasks)
    for task in all_tasks:
        if not any(task in t.ancestors for t in all_tasks):
            independent_tasks.add(task)
    return independent_tasks


def run_tasks(tasks, parallel=True):
    """Run tasks."""
    if parallel:
        _run_tasks_parallel(tasks)
    else:
        _run_tasks_sequential(tasks)


def _run_tasks_sequential(tasks):
    """Run tasks sequentially"""
    n_tasks = len(get_flattened_tasks(tasks))
    logger.info("Running %s tasks sequentially", n_tasks)

    for task in get_independent_tasks(tasks):
        task.run()


def _run_tasks_parallel(tasks):
    """Run tasks in parallel"""
    scheduled = list(get_flattened_tasks(tasks))
    running = []
    results = []

    n_scheduled, n_running = len(scheduled), len(running)
    n_tasks = n_scheduled

    pool = Pool()

    logger.info("Running %s tasks using at most %s processes", n_tasks,
                cpu_count())

    def done(task):
        """Assume a task is done if it not scheduled or running."""
        return not (task in scheduled or task in running)

    while scheduled or running:
        # Submit new tasks to pool
        just_scheduled = []
        for task in scheduled:
            if not task.ancestors or all(done(t) for t in task.ancestors):
                result = pool.apply_async(_run_task, [task])
                results.append(result)
                running.append(task)
                just_scheduled.append(task)
        for task in just_scheduled:
            scheduled.remove(task)

        # Handle completed tasks
        for task, result in zip(running, results):
            if result.ready():
                task.output_data = result.get()
                running.remove(task)
                results.remove(result)

        # Wait if there are still tasks running
        if running:
            time.sleep(0.1)

        # Log progress message
        if len(scheduled) != n_scheduled or len(running) != n_running:
            n_scheduled, n_running = len(scheduled), len(running)
            n_done = n_tasks - n_scheduled - n_running
            logger.info("Progress: %s tasks running or queued, %s tasks "
                        "waiting for ancestors, %s/%s done", n_running,
                        n_scheduled, n_done, n_tasks)

    pool.close()
    pool.join()


def _run_task(task):
    """Run task and return the result."""
    return task.run()
