###################################################################################################
# create XCLASS run script
###################################################################################################


def create_XCLASS_script(num, SSC, specie, alg_file, obs_file, mol_file, mode=None):

    resultpath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'results'))
    XCLASS = f"""import os
AlgorithmXMLFile = '{alg_file}'
experimentalData = '{obs_file}'
MolfitsFileName  = '{mol_file}'

newmolfit, modeldata, JobDir = myXCLASSFit()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')
"""

    # save to disk
    if mode==None:
        savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'XCLASS.py'))
    elif mode=='manual':
        savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'XCLASS_manual.py'))
    elif mode=='final':
        savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'XCLASS_final.py'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(XCLASS)
    return savepath


###################################################################################################
#
###################################################################################################
