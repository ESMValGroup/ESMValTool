'''
Module to control ftp''ing to collaboration server
'''

import os
import stat
import tempfile


def runftp(basename):
    '''ftp files to the collaboration twiki web server'''

    # Open a temporary file that will hold your ftp commands
    home = os.path.expanduser('~')
    hometmp = os.path.join(home, 'tmp')
    if not os.path.exists(hometmp):
        os.mkdir(hometmp)
    ftp_file = hometmp+'/'+basename+'.ftp'

    with open(ftp_file, 'w') as f:
        # Fill that file with ftp commands
        f.write('mkdir valid\n')
        f.write('cd valid\n')
        f.write('lcd /project/gmed/valnote/'+basename+'\n')
        # f.write('lcd /project/cma/valid/'+basename+'\n')
        f.write('mkdir '+basename+'\n')
        f.write('cd '+basename+'\n')
        f.write('mput *.html\n')
        f.write('mput *.css\n')
        f.write('mput *.txt\n')
        f.write('mput *.js\n')
        f.write('mkdir images\n')
        f.write('cd images\n')
        f.write('lcd images\n')
        f.write('mput *.png\n')
        f.write('mkdir ../thumbs\n')
        f.write('cd ../thumbs\n')
        f.write('lcd ../thumbs\n')
        f.write('mput *.png\n')
        f.write('bye\n')

    # Issue the ftp command (via els037)
    os.system('ssh els056 "cat '+ftp_file+' | ftp -i exxvmcollab"')
    os.remove(ftp_file)


def check_netrc(basename):
    '''Check .netrc file to make sure you have access to server'''
    extweb = 'x'
    netrcfile = os.path.join(os.path.expanduser('~'), '.netrc')
    if os.path.isfile(netrcfile):
        print
        print "FTPing files to collaborative twiki web server"
        filestat = os.stat(netrcfile)[0]

        # Check that only you have read and write permission
        if not ((filestat & stat.S_IRWXG) or (filestat & stat.S_IRWXO)):

            # Check to see if this is going on the development web or the
            #  user's personal web
            with open(netrcfile, 'r') as f:
                foundline = False
                for line in f:

                    # If you have found exxvmcollab then the next line
                    #  must contain the user login details
                    if foundline:

                        # If all the information is on seperate lines
                        #  then process the login line here
                        loginline = line.split()
                        login = loginline[1].split('-')
                        extweb = login[len(login)-1]
                        foundline = False

                    # If you have found the entry to do with the external
                    #  web server then turn on processing
                    if line.find('exxvmcollab') > 0:
                        foundline = True

                        # If all the information is all on one line then
                        #  process it here
                        if line.find('login') > 0:
                            loginline = line.split()
                            login = loginline[3].split('-')
                            extweb = login[len(login)-1]
                            foundline = False

            if extweb == 'x':
                print "~/.netrc file does not contain 'machine exxvmcollab' " \
                    + "and associated login details."

            else:
                # Do the transfer
                runftp(basename)

        else:
            print "Permissions not set up on ~/.netrc file. You should only " \
                + "have read and write permission for yourself and no-one " \
                + "else. Run 'chmod 600 ~/.netrc'"
    else:
        print "You do not have a .netrc file in your home directory. Create " \
            + "one, make sure only you have read and write permission " \
            + "(chmod 600 ~/.netrc) and put in the lines:"
        print "     machine exxvmcollab"
        print "             login DanCopsey-development"
        print "             password mypassword"

    return extweb
