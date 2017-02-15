"""
sample python script to be executed by ESMValTool via the command line
This will produce some sample output in the current directory
"""


print('')
print('WELCOME TO THE EXTERNAL PYTHON SCRIPT')
print('HELLO WORLD')

o=open('python_test_out.txt', 'w')
o.write('WELCOME TO THE EXTERNAL PYTHON SCRIPT\n')
o.write('HELLO WORLD')
o.close()

