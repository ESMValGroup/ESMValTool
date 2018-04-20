"""ESMValtool task definition"""
import contextlib
import datetime
import logging
import os
import pprint
import subprocess
import threading
import time
from multiprocessing import Pool, cpu_count

import psutil
import yaml

logger = logging.getLogger(__name__)

MODEL_KEYS = {
    'mip',
}


def which(executable):
    """Find executable in PATH."""
    for path in os.environ["PATH"].split(os.pathsep):
        if os.access(os.path.join(path, executable), os.X_OK):
            return os.path.join(path, executable)

    return None


def _get_resource_usage(process, start_time, children=True):
    """Get resource usage."""
    # yield header first
    entries = [
        'Date and time (UTC)',
        'Real time (s)',
        'CPU time (s)',
        'CPU (%)',
        'Memory (GB)',
        'Memory (%)',
        'Disk read (GB)',
        'Disk write (GB)',
    ]
    fmt = '{}\t' * len(entries[:-1]) + '{}\n'
    yield fmt.format(*entries)

    # Compute resource usage
    gigabyte = float(2**30)
    precision = [1, 1, None, 1, None, 3, 3]
    cache = {}
    while process.is_running():
        try:
            if children:
                # Include child processes
                processes = process.children(recursive=True)
                processes.append(process)
            else:
                processes = [process]

            # Update resource usage
            for proc in cache:
                # Set cpu percent and memory usage to 0 for old processes
                if proc not in processes:
                    cache[proc][1] = 0
                    cache[proc][2] = 0
                    cache[proc][3] = 0
            for proc in processes:
                # Update current processes
                cache[proc] = [
                    proc.cpu_times().user + proc.cpu_times().system,
                    proc.cpu_percent(),
                    proc.memory_info().rss / gigabyte,
                    proc.memory_percent(),
                    proc.io_counters().read_bytes / gigabyte,
                    proc.io_counters().write_bytes / gigabyte,
                ]
        except (OSError, psutil.AccessDenied, psutil.NoSuchProcess):
            # Try again if an error occurs because some process died
            continue

        # Create and yield log entry
        entries = [sum(entry) for entry in zip(*cache.values())]
        entries.insert(0, time.time() - start_time)
        entries = [round(entry, p) for entry, p in zip(entries, precision)]
        entries.insert(0, datetime.datetime.utcnow())
        yield fmt.format(*entries)


@contextlib.contextmanager
def resource_usage_logger(pid, filename, interval=1, children=True):
    """Log resource usage."""
    halt = threading.Event()

    def _log_resource_usage():
        """Write resource usage to file."""
        process = psutil.Process(pid)
        start_time = time.time()
        with open(filename, 'w') as file:
            for msg in _get_resource_usage(process, start_time, children):
                file.write(msg)
                time.sleep(interval)
                if halt.is_set():
                    return

    thread = threading.Thread(target=_log_resource_usage)
    thread.start()
    try:
        yield
    finally:
        halt.set()
        thread.join()


def write_ncl_settings(settings, filename, mode='wt'):
    """Write settings to NCL file."""
    logger.debug("Writing NCL configuration file %s", filename)

    def _format(value):
        """Format string or list as NCL"""
        if value is None or isinstance(value, str):
            txt = '"{}"'.format(value)
        elif isinstance(value, (list, tuple)):
            # TODO: convert None to fill value?
            # If an array contains a str, make all items str
            if any(isinstance(v, str) or v is None for v in value):
                value = [(str(v)) for v in value]
            txt = '(/{}/)'.format(', '.join(_format(v) for v in value))
        else:
            txt = str(value)
        return txt

    def _format_dict(name, dictionary):
        """Format dict as NCL"""
        lines = ['{} = True'.format(name)]
        for key, value in sorted(dictionary.items()):
            lines.append('{}@{} = {}'.format(name, key, _format(value)))
        txt = '\n'.join(lines)
        return txt

    def _header(name):
        """Delete any existing NCL variable known as `name`."""
        return ('if (isvar("{name}")) then\n'
                '    delete({name})\n'
                'end if\n'.format(name=name))

    lines = []
    for key, value in sorted(settings.items()):
        txt = _header(name=key)
        if isinstance(value, dict):
            txt += _format_dict(name=key, dictionary=value)
        else:
            txt += '{} = {}'.format(key, _format(value))
        lines.append(txt)
    with open(filename, mode) as file:
        file.write('\n\n'.join(lines))
        file.write('\n')


class AbstractTask(object):
    """Base class for defining task classes"""

    def __init__(self, settings, output_dir, ancestors=None):
        """Initialize task."""
        self.settings = settings
        self.ancestors = [] if ancestors is None else ancestors
        self.output_dir = output_dir
        self.output_files = None

    def flatten(self):
        """Return a flattened set of all ancestor tasks and task itself."""
        tasks = set()
        for task in self.ancestors:
            tasks.update(task.flatten())
        tasks.add(self)
        return tasks

    def run(self, input_files=None):
        """Run task."""
        if not self.output_files:
            if input_files is None:
                input_files = []
            for task in self.ancestors:
                input_files.extend(task.run())
            self.output_files = self._run(input_files)

        return self.output_files

    def _run(self, input_files):
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

    def __init__(self, script, settings, output_dir, ancestors=None):
        """Initialize"""
        super(DiagnosticTask, self).__init__(
            settings=settings, output_dir=output_dir, ancestors=ancestors)
        self.script = script
        self.cmd = self._initialize_cmd(script)
        self.log = os.path.join(settings['run_dir'], 'log.txt')
        self.resource_log = os.path.join(settings['run_dir'],
                                         'resource_usage.txt')

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
                'ncl': [which('ncl'), '-n', '-p'],
                'r': [which('Rscript'), '--slave', '--quiet'],
            }
            if extension not in executables:
                raise DiagnosticError(
                    "Cannot execute script {} ({}): non-executable file "
                    "with unknown extension.".format(script, script_file))

            cmd = executables[extension]

        cmd.append(script_file)

        return cmd

    def write_settings(self):
        """Write settings to file"""
        run_dir = self.settings['run_dir']
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)

        filename = os.path.join(run_dir, 'settings.yml')

        with open(filename, 'w') as file:
            yaml.safe_dump(self.settings, file)

        # If running an NCL script:
        if self.script.lower().endswith('.ncl'):
            # Also write an NCL file and return the name of that instead.
            return self._write_ncl_settings()

        return filename

    def _write_ncl_settings(self):
        """Write settings to NCL file"""
        filename = os.path.join(self.settings['run_dir'], 'settings.ncl')

        config_user_keys = {
            'run_dir',
            'plot_dir',
            'work_dir',
            'max_data_filesize',
            'output_file_type',
            'log_level',
            'write_plots',
            'write_netcdf',
        }
        settings = {'diag_script_info': {}, 'config_user_info': {}}
        for key, value in self.settings.items():
            if key in config_user_keys:
                settings['config_user_info'][key] = value
            elif not isinstance(value, dict):
                settings['diag_script_info'][key] = value
            else:
                settings[key] = value

        write_ncl_settings(settings, filename)

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
                logger.warning("NCL: %s", line)
                warned = True
            for error in errors:
                if error in line:
                    logger.error(msg)
                    logger.error("NCL: %s", line)
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

    def _run(self, input_files):
        """Run the diagnostic script."""
        if self.script is None:  # Run only preprocessor
            output_files = []
            return output_files

        is_ncl_script = self.script.lower().endswith('.ncl')
        if is_ncl_script:
            input_files = [
                f for f in input_files
                if f.endswith('.ncl') or os.path.isdir(f)
            ]
        else:
            input_files = [
                f for f in input_files
                if f.endswith('.yml') or os.path.isdir(f)
            ]

        self.settings['input_files'] = input_files

        cmd = list(self.cmd)
        cwd = None
        env = self.settings.pop('env', None)
        if env:
            env = {str(k): str(v) for k, v in env.items()}

        settings_file = self.write_settings()

        if is_ncl_script:
            cwd = os.path.dirname(__file__)
            env = dict(os.environ)
            env['settings'] = settings_file
        else:
            cmd.append(settings_file)

        logger.info("Running command %s", cmd)
        logger.debug("in environment\n%s", pprint.pformat(env))
        logger.debug("in current working directory: %s", cwd)
        logger.info("Writing output to %s", self.output_dir)
        logger.info("Writing plots to %s", self.settings['plot_dir'])
        logger.info("Writing log to %s", self.log)

        rerun_msg = '' if cwd is None else 'cd {}; '.format(cwd)
        if env:
            rerun_msg += ' '.join('{}="{}"'.format(k, env[k]) for k in env
                                  if k not in os.environ)
        rerun_msg += ' ' + ' '.join(cmd)
        logger.info("To re-run this diagnostic script, run:\n%s", rerun_msg)

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=cwd,
            env=env)

        returncode = None
        last_line = ['']

        with resource_usage_logger(process.pid, self.resource_log),\
                open(self.log, 'at') as log:
            while returncode is None:
                returncode = process.poll()
                txt = process.stdout.read()
                txt = txt.decode(encoding='utf-8', errors='ignore')
                log.write(txt)

                # Check if an error occurred in an NCL script
                # Last line is treated separately to avoid missing
                # error messages spread out over multiple lines.
                lines = txt.split('\n')
                if is_ncl_script:
                    self._control_ncl_execution(process, last_line + lines)
                last_line = lines[-1:]

                # wait, but not long because the stdout buffer may fill up:
                # https://docs.python.org/3.6/library/subprocess.html#subprocess.Popen.stdout
                time.sleep(0.001)

        if returncode == 0:
            return [self.output_dir]

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


def run_tasks(tasks, max_parallel_tasks=None):
    """Run tasks."""
    if max_parallel_tasks == 1:
        _run_tasks_sequential(tasks)
    else:
        _run_tasks_parallel(tasks, max_parallel_tasks)


def _run_tasks_sequential(tasks):
    """Run tasks sequentially"""
    n_tasks = len(get_flattened_tasks(tasks))
    logger.info("Running %s tasks sequentially", n_tasks)

    for task in get_independent_tasks(tasks):
        task.run()


def _run_tasks_parallel(tasks, max_parallel_tasks=None):
    """Run tasks in parallel"""
    scheduled = get_flattened_tasks(tasks)
    running = []
    results = []

    n_scheduled, n_running = len(scheduled), len(running)
    n_tasks = n_scheduled

    pool = Pool(processes=max_parallel_tasks)

    logger.info("Running %s tasks using at most %s processes", n_tasks,
                max_parallel_tasks or cpu_count())

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
                task.output_files = result.get()
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
