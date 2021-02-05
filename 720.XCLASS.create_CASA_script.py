###################################################################################################
# run XCLASS within CASA
###################################################################################################

def create_CASA_script(XCL_files, savepath=os.path.join(XCLASSdir,'all_XCLASS.py')):

    with open(savepath, 'w') as f:
        f.write(f"""
import os
import time,datetime

if (os.environ['CONDA_DEFAULT_ENV'] != 'python2' ):
    raise EnvironmentError(u"\u001b[31mYOU NEED TO RUN CASA WITHIN A PYTHON2 ENVIRONMENT!\u001b[0m")
try:
    print(casa_util.version_string())
except:
    raise EnvironmentError(u"\u001b[31mYOU NEED TO RUN XCLASS IN CASA!\u001b[0m")

start = time.time()
for f in {XCL_files}:
    execfile(f)

stop = time.time()
exec_time = np.round(stop-start, 1)
casalog.post("\\nFinished XCLASS runs \\nExecution took "+str(datetime.timedelta(seconds=exec_time))+"hours.\\n")
""")

    print("Run this command\nexecfile('"+savepath+"')\nin CASA in a python2 environment to perform all fits.")


###################################################################################################
#
###################################################################################################
