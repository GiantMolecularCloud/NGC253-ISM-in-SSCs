###################################################################################################
# create model spectrum
###################################################################################################

def create_model_spectrum(num, SSC, detected_species, mode='final'):
    """
    Load molfit files and stitch them up to get a molfit file for all species.
    """

    if mode=='final':
        modestr = '_final'
    elif mode=='manual':
        modestr = '_manual'
    else:
        modestr = ''

    # load molfit files
    molfit_files = []
    for specie in detected_species:
        molfit_files.append(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'results','molecules'+modestr+'__LM__call_1.out.molfit')))
    molfit_files.sort()

    combined_molfit = ''
    for f in molfit_files:

        # get specie
        sp = [s for s in f.split('/') if '=' in s][0]
        specie = [s for s in unique_species if escape_fname(s)==sp][0]

        # handle specie names correctly
        if '#1' in specie:
            spx = specie
        else:
            spx = specie+';'

        # load molfit file
        molfit = open(f, 'r').readlines()
        for idx,ll in enumerate(molfit):
            ll = ll.replace('\n','')
            ll = [l for l in ll.split(' ') if not l=='']
            molfit[idx] = ll

        # get which lines contain the components of this specie
        for idx,ll in enumerate(molfit):
            if ll[0]==spx:
                idx_start = idx+1
            if ll[0]!='n' and idx>1:
                idx_stop  = idx
                break
            if idx==len(molfit)-1:
                idx_stop  = idx+1

        try:
            idx_start, idx_stop
        except:
            print(SSC, specie)

        for idx in np.arange(idx_start-1, idx_stop):
            combined_molfit += '  '.join(molfit[idx]).replace('y','n')
            combined_molfit += '\n'

    savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','all_molecules.run_'+str(num)+'.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(combined_molfit)
    return savepath


def model_spectrum_XCLASS(num, SSC, model_molfit_file):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    savepath   = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','XCLASS_model_spectrum.run_'+str(num)+'.py'))
    resultpath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_'+str(num),'combined_model'))
    modelpath  = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_'+str(num),'combined_model.spectrum.pickle'))
    energypath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_'+str(num),'combined_model.energy.pickle'))
    taupath    = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_'+str(num),'combined_model.optical.pickle'))

    with open(savepath, 'w') as f:
        f.write(f"""
import pickle

FreqMin          = 342500
FreqMax          = 358500
FreqStep         = 2.5
TelescopeSize    = 0.15
Inter_Flag       = True
t_back_flag      = True
tBack            = 0.0
tslope           = 0.0
nH_flag          = False
N_H              = 1e+24
beta_dust        = 0.1
kappa_1300       = 0.01
MolfitsFileName  = '{model_molfit_file}'
iso_flag         = False
IsoTableFileName = ' '
RestFreq         = 0.0
vLSR             = 0.0

modeldata, log, TransEnergies, IntOptical, JobDir = myXCLASS()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')

pickle.dump(modeldata, open('{modelpath}', 'w'))
pickle.dump(TransEnergies, open('{energypath}', 'w'))
pickle.dump(IntOptical, open('{taupath}', 'w'))
""")
    return savepath


###################################################################################################
#
###################################################################################################
