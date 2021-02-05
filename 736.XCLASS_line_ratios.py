#######################
# GAS IN SSCS: XCLASS #
#######################


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))

Xfiles      = fnunpickle(os.path.join(Xfinaldir, 'Xfiles.pickle'))
temperature = fnunpickle(os.path.join(Xfinaldir, 'temperature.pickle'))
nums        = fnunpickle(os.path.join(Xfinaldir, 'nums.pickle'))
final_data = fnunpickle(os.path.join(Xfinaldir, 'data.pickle'))

execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))


###################################################################################################
# model independently
###################################################################################################

def create_line_molfit(num, SSC):

    # fit molfit files & model molfit files
    fit   = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'fit','results','molecules__LM__call_1.out.molfit'))
    model = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','molecules.molfit'))

    # read molfit file
    with open(fit, 'r') as file :
        filedata = file.read()

    # replace the fit switch
    filedata = filedata.replace(' y ', ' n ')

    # write the model molfit file
    os.system('mkdir -p '+os.path.dirname(model))
    with open(model, 'w') as file:
        file.write(filedata)

    Xfiles[SSC['num']][str(num)]['model'] = {'molfit': model}


def model_spectrum_XCLASS(num, SSC):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    modelpath  = Xfiles[SSC['num']][str(num)]['model']['molfit']
    savepath   = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','XCLASS_model.py'))
    resultpath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','results'))
    specpath   = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','spectrum.pickle'))
    energypath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','energy.pickle'))
    taupath    = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','opacity.pickle'))

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
MolfitsFileName  = '{modelpath}'
iso_flag         = False
IsoTableFileName = ' '
RestFreq         = 0.0
vLSR             = 0.0

modeldata, log, TransEnergies, IntOptical, JobDir = myXCLASS()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')

pickle.dump(modeldata, open('{specpath}', 'w'))
pickle.dump(TransEnergies, open('{energypath}', 'w'))
pickle.dump(IntOptical, open('{taupath}', 'w'))
""")
    Xfiles[SSC['num']][str(num)]['model']['file']     = savepath
    Xfiles[SSC['num']][str(num)]['model']['spectrum'] = specpath
    Xfiles[SSC['num']][str(num)]['model']['energy']   = energypath
    Xfiles[SSC['num']][str(num)]['model']['tau']      = taupath


for num in tqdm(nums):
    for SSC in SSCs:
        create_model_molfit(num, SSC)
        model_spectrum_XCLASS(num, SSC)


###################################################################################################
# model independently
###################################################################################################

def CASA_command(SSC):
    XCL_files = [Xfiles[SSC['num']][str(num)]['model']['file'] for num in nums]
    Xcommand = f"""
import time,datetime

start = time.time()
for f in {XCL_files}:
\texecfile(f)

stop = time.time()
exec_time = np.round(stop-start, 1)
casalog.post("\\nFinished XCLASS runs \\nExecution took "+str(datetime.timedelta(seconds=exec_time))+"hours.\\n")
"""
    return Xcommand

commands = [CASA_command(SSC) for SSC in SSCs]
pool = Pool(len(SSCs))
pool.map(run_in_casa, commands)
pool.close()
pool.join()

fnpickle(Xfiles, os.path.join(Xfinaldir, 'Xfiles.pickle'))


###################################################################################################
#
###################################################################################################
