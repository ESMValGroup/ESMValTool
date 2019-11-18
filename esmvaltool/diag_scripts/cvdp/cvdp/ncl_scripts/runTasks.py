import subprocess
import sys
import time
import os

#--------------------------- BEGIN USER MODIFICATIONS ---------------------------
EXEC_STR = "ncl -Q"
POLL_INTERVAL = 15.          # seconds
MAX_CONCURRENT =  int(os.environ.get('MAX_TASKS'))    #  previous settings:  = 4 or  = sys.maxint
#--------------------------- END USER MODIFICATIONS -----------------------------

def launchTask(script):
#    print "Launching: ", script
    task = subprocess.Popen(EXEC_STR + " " + script, shell=True, executable="/bin/bash") 
    return task
 
# ------------------------- main -----------------------------------------------
 
# get command-line args, strip out 1st element, which is the name of this script...
scripts = sys.argv[1:]
#print scripts       # debugging -- remove or comment out 

# fire off up-to MAX_CONCURRENT subprocesses...
tasks = list()
for i,script in enumerate(scripts):
    if i >= MAX_CONCURRENT:
        break
    tasks.append( launchTask(script) )

#print scripts
scripts = scripts[len(tasks):]  # remove those scripts we've just launched...
#print scripts

#for task in tasks:     # debugging -- remove or comment out 
#    print(task.pid)

while len(tasks) > 0:
    finishedList = []
    for task in tasks:
         retCode = task.poll()
         if retCode != None:
#             print "Task status ", task.pid, ": ", task.poll()
             finishedList.append(task)

             # more scripts to be run?
             if len(scripts) > 0:
                 tasks.append( launchTask(scripts[0]) )
                 del scripts[0]

    for task in finishedList:
        tasks.remove(task)

    time.sleep(POLL_INTERVAL)
#    print "."      # Feedback to show the script is doing something; not necessary

print "runTasks.py: Done with CVDP calculation scripts"


