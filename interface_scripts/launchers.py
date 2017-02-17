from auxiliary import nclExecuteError, nclExecuteWarning, error, info
import os
import pdb
import re
import subprocess
import sys
import string
import StringIO
import contextlib

@contextlib.contextmanager
def stdoutIO(std_outerr=None):
    """
    Routine to capture stdout/err from python routines
    """
    old = sys.stdout
    old_err = sys.stderr
    if std_outerr is None:
        std_outerr = StringIO.StringIO()
    sys.stdout = sys.stderr = std_outerr
    yield std_outerr
    sys.stdout = old
    sys.stderr = old_err


class launchers(object):
    def __init__(self):
        self.persistent_env_variables = ['ESMValTool_data_root']
        self.filename = os.path.join(os.path.dirname(__file__),                   
                                        '../interface_data/curr_trace_indent.txt') 
    def convert_arguments(self):
        """
        convert launcher arguments to a dictionary

        If the launcher had been assigned an attribute 'arguments'
        then this is assumed to be a string which can be converted to a
        dictionary using the dict() command.

        A simple example would be [('a': 'b'), ('run_in_shell' : True)]
        which would give a dictionary like
        {'a' : 'b', 'run_in_shell' : True}

        In case that no argument is avalaible, the launch_args is set to None
        """
        self.launch_args = None
        if hasattr(self, 'arguments'):
            # the launcher arguments are a string which is interpreted as the
            # input to the dict() command. As a result one obtains a dictionary
            exec 'self.launch_args = dict(' + self.arguments + ')'

    def write_stdouterr(self, std_outerr, verbosity, exit_on_warning):
        """ @brief Filder and write subprocess output to screen
            @param std_outerr Array with stdout/stderr mixed
            @param verbosity Integer controling the verbosity of the output
            @param exit_on_warning Boolean defining whether the wrapper should
                                   exit on warnings
        """
        # Supress some warnings if they've been written to the ignore file
        # (this should be done by the diagnostic script)
        warnings_ignore_file = "interface_data/warnings_to_ignore.txt"
        if os.path.isfile(warnings_ignore_file):
            fin = open(warnings_ignore_file)
            diag_specific_warnings_to_ignore = fin.readlines()
        else:
            diag_specific_warnings_to_ignore = []

        # Suppress version info output and empty lines
        info = [line for line in range(len(std_outerr)) if ((line > self.filter_max_line) and len(std_outerr[line]) > 0)]

        fatal = [line for line in range(len(std_outerr))
                 if re.search(self.fatal_string, std_outerr[line], re.IGNORECASE)]

        error = [line for line in range(len(std_outerr))
                 if re.search(self.error_string, std_outerr[line], re.IGNORECASE)]

        warning = [line for line in range(len(std_outerr))
                   if re.search(self.warning_string, std_outerr[line], re.IGNORECASE)]

        # Remove warnings to be skipped
        warnings_to_skip = []
        if len(diag_specific_warnings_to_ignore) == 0:
            warnings_to_keep = warning
        else:
            warnings_to_keep = []
            for line in warning:
                if std_outerr[line].strip() in [w.strip() for w in diag_specific_warnings_to_ignore]:
                    warnings_to_skip.append(line)
                else:
                    warnings_to_keep.append(line)

        if len(fatal) > 0:
            for line in [currLine for currLine in std_outerr]:
                sys.stderr.write(self.lang + " ERROR MESSAGE: " + line + '\n')
            raise nclExecuteError(self.lang + " ERROR (see full NCL output above)")

        if len(error) > 0:
            for line in [currLine for currLine in std_outerr]:
                sys.stderr.write(self.lang + " ERROR MESSAGE: " + line + '\n')
# A-laue_ax+
            # In contrast to "fatal" errors, we treat "normal" errors as warnings
            # and do not abort (except if "exit_on_warning" in the namelist is set
            # to "True".)
            if exit_on_warning:
# A-laue_ax-
                raise nclExecuteError(self.lang + " ERROR (see full NCL output above)")

        all_output_as_warning = False
        if len(warnings_to_keep) > 0:
            all_output_as_warning = True
            for line in [currLine for currLine in std_outerr]:
                sys.stderr.write(self.lang + " WARNING MESSAGE: " + line + '\n')
            if exit_on_warning:
                raise nclExecuteWarning(self.lang + " WARNING (see full NCL output above)")

        if len(info) > 0 and verbosity <= 10 and not all_output_as_warning:
            for line_no in info:
                remove_quotation_marks = re.sub('"', '', std_outerr[line_no])
                if re.search('.*info: (.*)', remove_quotation_marks) is not None:
                    output_string = re.search(".*info: (.*)", remove_quotation_marks).group(1)
                    sys.stdout.write(self.lang + " info: " + output_string + '\n')
                else:
                    sys.stdout.write(std_outerr[line_no] + '\n')

        if verbosity > 10:
            for line in [currLine for currLine in std_outerr]:
                sys.stdout.write(line + '\n')

        if len(warnings_to_skip) > 0:
            sys.stdout.write(self.lang + " info: The following warnings were ignored (specified to" + '\n')
            sys.stdout.write(self.lang + " info: be ignored in the diagnostic script)" + '\n')
            for line in warnings_to_skip:
                sys.stdout.write(self.lang + " info: " + std_outerr[line] + '\n')


class ncl_launcher(launchers):
    def __init__(self):
        launchers.__init__(self)
        self.lang = "NCL"
        self.filter_max_line = 4
        self.error_string = 'error:'
        self.fatal_string = 'fatal:'
        self.warning_string = 'warning:'

    def execute(self, ncl_executable, project_info, verbosity, exit_on_warning):
        """ @brief Wrapper to execute NCL scripts
            @param ncl_command Full path to the NCL script to execute
            @param project_info Current namelist in dictionary format
            @param verbosity Set the verbosity level (0 minimum verbosity)
            @param exit_on_warning Boolean defining whether the wrapper should
            crash on any NCL warnings

            This wrapper will take an NCL script, execute it then scan the
            stdout for the keywords 'fatal' and 'warning'. If they occur an
            exception is raised. The wrapper will also delete all
            'ESMValTool_'-prefixed environment variables after execution.
        """
        # Reset NCL trace back indent (available with verbosity=2)
	f_ncl_indent = open(self.filename, "w")
        f_ncl_indent.write("0")
        f_ncl_indent.close()

        if not os.path.exists(ncl_executable):
            raise nclExecuteError(self.lang + " ERROR (file to execute, \""
                                            + ncl_executable
                                            + "\", is missing)")

        run_application = subprocess.Popen("ncl " + ncl_executable, shell=True,
                                           stdin=open(os.devnull),
                                           stdout=subprocess.PIPE)
        #run_application.wait()
        std_outerr = run_application.communicate()[0].split('\n')
        self.write_stdouterr(std_outerr, verbosity, exit_on_warning)

        for key in [env for env in os.environ if re.search('^ESMValTool_*', env)]:
            if key not in self.persistent_env_variables:
                del(os.environ[key])


class r_launcher(launchers):
    def __init__(self, **kwargs):
        super(r_launcher, self).__init__(**kwargs)
        self.filter_max_line = 0
        self.lang = "R  "
        self.error_string = 'Error '
        self.fatal_string = 'fatal:'
        self.warning_string = 'Warning '

    def execute(self, r_script, project_info, verbosity, exit_on_warning):
        """ @brief Wrapper to execute R scripts
            @param r_script Full path to the R script to execute
            @param project_info Current namelist in dictionary format
            @param verbosity Set the verbosity level (0 minimum verbosity)
            @param exit_on_warning Boolean defining whether the wrapper should
            crash on any R warnings

            This wrapper will take an R script and executes it.
            The wrapper will also delete all 'ESMValTool_'-prefixed
            environment variables after execution.
        """
        # Reset NCL trace back indent (available with verbosity=2)
	f_r_indent = open(self.filename, "w")
        f_r_indent.write("0" + '\n')
        f_r_indent.close()

        if 'r_pre_launch' in project_info['GLOBAL']:
            r_pre_launch = project_info['GLOBAL']['r_pre_launch']
        else:
            r_pre_launch = ""

        if not os.path.exists(r_script):
            raise nclExecuteError(self.lang + " ERROR (file to execute, \""
                                            + r_script
                                            + "\", is missing)")

        # Default launch command for R scripts
        r_launch = ' Rscript --slave --quiet '

        # Namelist override for R script launch command
        self.convert_arguments()
        if self.launch_args is not None:
            if 'r_launch' in self.launch_args.keys():
                r_launch = self.launch_args['r_launch']

        r_run = r_pre_launch + r_launch + r_script

        run_application = subprocess.Popen(r_run, shell=True,
                                           stdin=open(os.devnull),
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT)
        run_application.wait()
        std_outerr = run_application.communicate()[0].split('\n')
        self.write_stdouterr(std_outerr, verbosity, exit_on_warning)

        for key in [env for env in os.environ if re.search('^ESMValTool_*', env)]:
            del(os.environ[key])


class py_launcher(launchers):
    """
    general python launcher
    """
    def __init__(self, execute_as_shell=False, **kwargs):
        """
        Parameters
        ----------
        execute_as_shell : bool
            execute launcher in shell using e.g. subprocess
            if False, then the diagnostic script is called directly
        """
        super(py_launcher, self).__init__(**kwargs)
        self.lang = "PY "
        self.filter_max_line = 4
        self.error_string = 'error:'
        self.fatal_string = 'fatal:'
        self.warning_string = 'warning:'
        self.execute_as_shell = execute_as_shell

    def execute(self, python_executable, project_info, verbosity, exit_on_warning):
        """ @brief Wrapper to execute PYTHON scripts
            @param python_executable: Full path to the python script to execute
            @param project_info Current namelist in dictionary format
            @param verbosity Set the verbosity level (0 minimum verbosity)
            @param exit_on_warning Boolean defining whether the wrapper should
            crash on any python warnings

            This wrapper will take a PYTHON script, execute it then scan the
            stdout for the keywords 'fatal' and 'warning'. If they occur an
            exception is raised. The wrapper will also delete all
            'ESMValTool_'-prefixed environment variables after execution.

            An alternative way to launch python scripts would be to call them
            directly. This would imply however an import of the diagnostic
            python script at this point here. To be more flexible, an execution
            via the command line is preferred. See here for details:

            * http://stackoverflow.com/questions/4230725/how-to-execute-a-python-script-file-with-an-argument-from-inside-another-python
            * http://docs.python.org/2/tutorial/modules.html
        """

        # check if additional arguments were given
        self.convert_arguments()
        if self.launch_args is not None:
            if 'execute_as_shell' in self.launch_args.keys():
                self.execute_as_shell = self.launch_args['execute_as_shell']

        if self.execute_as_shell:
            self._execute_shell(python_executable, project_info, verbosity, exit_on_warning)
        else:  # Default option
            self._execute_script(python_executable, project_info, verbosity, exit_on_warning)

    def _execute_script(self, python_executable, project_info, verbosity, exit_on_warning):
        """
        excute python script directly

        Parameters
        ----------
        python_executable : str
            name of executable script
        project_info : dict
            dictionary with relevant project info from namelist
        """
        if not os.path.exists(python_executable):
            raise ValueError('Python executable not existing: %s' % python_executable)

        # try to import the script. Note the script needs to be in
        # the pythonpath. This is normally ensured already due to
        # the sys.path.append statements in main.py
        cmd = 'import ' + os.path.splitext(os.path.basename(python_executable))[0] + ' as usr_script'

        try:
            exec cmd
        except:
            print cmd
            raise ValueError('The script %s can not be imported!' % python_executable)

        # import was successfull. Now call the script with project_info
        # as argument
        # Capturing stdout/err for warnings, however, this currently doesn't capture Tracebacks..
#        with stdoutIO() as s:
#            usr_script.main(project_info)
#        self.write_stdouterr(string.split(s.getvalue(), '\n'), verbosity, exit_on_warning)

        # This catpures Traceback but won't let us analyse the stdout/err for text warnings
        usr_script.main(project_info)

    def _execute_shell(self, python_executable, project_info, verbosity, exit_on_warning):
        """
        execute python script in shell as subprocess
        """
        run_application = subprocess.Popen("python " + python_executable, shell=True,
                                           stdin=open(os.devnull),
                                           stdout=subprocess.PIPE)
        run_application.wait()
        std_outerr = run_application.communicate()[0].split('\n')
        self.write_stdouterr(std_outerr, verbosity, exit_on_warning)

        for key in [env for env in os.environ if re.search('^ESMValTool_*', env)]:
            del(os.environ[key])

class shell_launcher(launchers):                                                           
      """ @brief general unix shell launcher                                             
            @param shell: name of shell (bash or csh)                                      
      """                                                                                  
      def __init__(self,shell):                                                            
              if shell in ['csh','bash']:                                                  
                      self.lang = shell                                                    
              else:                                                                        
                      raise error('Unknown shell: {0}'.format(shell))                      
              super(shell_launcher, self).__init__()                                       
                                                                                           
      def execute(self, executable, project_info, verbosity, exit_on_warning):             
              try:                                                                         
                      with open(self.filename, "w") as f:                                  
                              f.write("0")                                                 
                              f.close()                                                    
              except IOError:                                                              
                      raise error("IOError while open/write: '{0}'".format(filename))      
                                                                                           
              if not os.path.exists(executable):                                           
                      raise IOError("file to execute is missing: {0}".format(executable))  
                                                                                           
              cmd = self.lang + " {0}".format(executable)                                  
              run_application = subprocess.Popen(cmd, shell=True,                          
                                           stdin=open(os.devnull),                         
                                           stdout=subprocess.PIPE,                         
                                           stderr=subprocess.PIPE)                         
              run_application.wait()                                                       
              output = run_application.communicate()                                       
              std_out = filter(None,output[0].split('\n'))                                 
              std_err = filter(None,output[1].split('\n'))                                 
              for i in [self.lang.upper() + ' INFO : ' + item for item in std_out]:        
                      print i                                                              
              for i in [self.lang.upper() + ' ERROR: ' + item for item in std_err]:        
                                                                                           
                      print i                                                              
                                                                                           
              for key in [env for env in os.environ if re.search('^ESMValTool_*', env)]:   
                      del(os.environ[key])                                                 
                                                                                           
class csh_launcher(shell_launcher):                                                        
      """ @brief csh-shell script launcher                                               
      """                                                                                  
      def __init__(self):                                                                  
              super(csh_launcher, self).__init__('csh')                                    
                                                                                           
class bash_launcher(shell_launcher):                                                       
      """ @brief bash-shell script launcher                                              
      """                                                                                  
      def __init__(self):                                                                  
              super(bash_launcher, self).__init__('bash')                                  
