#######################
# GAS IN SSCS: XCLASS #
#######################


###################################################################################################
# build data structure
###################################################################################################

# save number of components
for SSC in SSCs:
    for specie, data in data_XCLASS_refit[SSC['num']].items():
        ncomp = len(data['velocity']['median'])
        data_XCLASS_refit[SSC['num']][specie]['components'] = ncomp


# create structure to save XCLASS information
intensity_files = {}
for SSC in SSCs:
    try:
        intensity_files[SSC['num']]
    except:
        intensity_files[SSC['num']] = {}

    for specie in data_XCLASS_refit[SSC['num']].keys():
        try:
            intensity_files[SSC['num']][specie]
        except:
            intensity_files[SSC['num']][specie] = {}

        for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
            try:
                intensity_files[SSC['num']][specie][component]
            except:
                intensity_files[SSC['num']][specie][component] = {}

            for num in nums:
                try:
                    intensity_files[SSC['num']][specie][component][num]
                except:
                    intensity_files[SSC['num']][specie][component][num] = {}


###################################################################################################
# create model molfit file
###################################################################################################

def create_model_molfit(num, SSC, specie, component):

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset"""

    if '#' in specie:
        spx = specie
    else:
        spx = specie+';'

    # fix ncomp to 1 because only a single component is modeled in each run
    molfit += f"""
{spx}   1"""

    T = data_XCLASS_refit[SSC['num']][specie]['temperature']['all'][component][num]
    v = data_XCLASS_refit[SSC['num']][specie]['velocity']['all'][component][num]
    w = data_XCLASS_refit[SSC['num']][specie]['linewidth']['all'][component][num]
    N = data_XCLASS_refit[SSC['num']][specie]['column density']['all'][component][num]

    molfit += f"""
n  0.01  10  1   n  {'{:.2E}'.format(T)}  {'{:.2E}'.format(T)}  {'{:.2E}'.format(T)}   n  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}   n  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}   n  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}   c"""

    # save to disk
    savepath = escape_fname(os.path.join(intensitydir,'SSC_'+str(SSC['no']),specie,'component_'+str(component),'run_'+str(num),'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)

    # keep savepath
    intensity_files[SSC['num']][specie][component][num]['molfit'] = savepath


for num in tqdm(nums):
    for SSC in SSCs:
        for specie,data in data_XCLASS_refit[SSC['num']].items():
            for component in np.arange(data['components']):
                create_model_molfit(num, SSC, specie, component)


###################################################################################################
# create model script
###################################################################################################

def model_spectrum_XCLASS(num, SSC, specie, component):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    modelpath  = intensity_files[SSC['num']][specie][component][num]['molfit']
    savepath   = escape_fname(os.path.join(intensitydir,'SSC_'+str(SSC['no']),specie,'component_'+str(component),'run_'+str(num),'XCLASS_model.py'))
    resultpath = escape_fname(os.path.join(intensitydir,'SSC_'+str(SSC['no']),specie,'component_'+str(component),'run_'+str(num),'results'))
    specpath   = escape_fname(os.path.join(intensitydir,'SSC_'+str(SSC['no']),specie,'component_'+str(component),'run_'+str(num),'spectrum.pickle'))
    energypath = escape_fname(os.path.join(intensitydir,'SSC_'+str(SSC['no']),specie,'component_'+str(component),'run_'+str(num),'energy.pickle'))
    taupath    = escape_fname(os.path.join(intensitydir,'SSC_'+str(SSC['no']),specie,'component_'+str(component),'run_'+str(num),'opacity.pickle'))
    intensity_files[SSC['num']][specie][component][num]['scriptfile'] = savepath
    intensity_files[SSC['num']][specie][component][num]['spectrum']   = specpath
    intensity_files[SSC['num']][specie][component][num]['energy']     = energypath
    intensity_files[SSC['num']][specie][component][num]['tau']        = taupath

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


for num in tqdm(nums):
    for SSC in SSCs:
        for specie,data in data_XCLASS_refit[SSC['num']].items():
            for component in np.arange(data['components']):
                model_spectrum_XCLASS(num, SSC, specie, component)


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSC,specie,component):
    XCL_files = [intensity_files[SSC['num']][specie][component][num]['scriptfile'] for num in nums]
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

commands = []
for SSC in SSCs:
    for specie in data_XCLASS_refit[SSC['num']].keys():
        for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
            commands.append(CASA_command(SSC,specie,component))

pool = Pool(len(SSCs))
pool.map(run_in_casa, commands)
pool.close()
pool.join()

fnpickle(intensity_files, join(intensitydir, 'intensity_files.pickle'))


###################################################################################################
# load model spectra
###################################################################################################

# XCLASS does not allow to model only certain transitions but always models all transitions of a
# given species. To overcome this: fit the model spectrum with Gaussians


def specie_lines(specie):
    """
    Select all lines of a given species.
    """
    blines = []
    for line in lines:
        if line['XCLASS']==specie:
            blines.append(line)
    return blines


# create structure to save data
###################################################################################################

intensity_spectra = {}
for SSC in SSCs:
    try:
        intensity_spectra[SSC['num']]
    except:
        intensity_spectra[SSC['num']] = {}

    for specie in data_XCLASS_refit[SSC['num']].keys():
        try:
            intensity_spectra[SSC['num']][specie]
        except:
            intensity_spectra[SSC['num']][specie] = {}

        for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
            try:
                intensity_spectra[SSC['num']][specie][component]
            except:
                intensity_spectra[SSC['num']][specie][component] = {}

            for num in nums:
                try:
                    intensity_spectra[SSC['num']][specie][component][num]
                except:
                    intensity_spectra[SSC['num']][specie][component][num] = {}


# read in spectra
###################################################################################################

for num in tqdm(nums):
    for SSC in SSCs:
        for specie in data_XCLASS_refit[SSC['num']].keys():
            for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
                with open(intensity_files[SSC['num']][specie][component][num]['spectrum'], 'rb') as f:
                    m = pickle.load(f, encoding="latin1")
                    f = (m[:,0]*u.MHz).to(u.GHz)
                    i = (m[:,1]*u.K)
                    intensity_spectra[SSC['num']][specie][component][num]['model spectrum'] = [f,i]
                # print_overwrite("loaded num "+str(num)+"/101, SSC "+str(SSC['num'])+"/14, "+specie, end='\r')


# cut out bands
###################################################################################################

spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
for SSC in tqdm(SSCs):
    for band in ['LSB','USB']:
        obsfreq = spectra[SSC['num']][band]['frequency'].to(u.GHz)
        shiftfreq = [(-vsys-SSC['velshift']).to(u.GHz, equivalencies=u.doppler_optical(f)).value for f in obsfreq]*u.GHz
        fmin = np.nanmin(shiftfreq)
        fmax = np.nanmax(shiftfreq)
        for specie in data_XCLASS_refit[SSC['num']].keys():
            for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
                for num in nums:
                    frequency = intensity_spectra[SSC['num']][specie][component][num]['model spectrum'][0]
                    intensity = intensity_spectra[SSC['num']][specie][component][num]['model spectrum'][1]
                    bandmask = np.logical_and(frequency>fmin, frequency<fmax)
                    intensity_spectra[SSC['num']][specie][component][num][band] = [frequency[bandmask], intensity[bandmask]]

fnpickle(intensity_spectra, join(intensitydir, 'intensity_spectra.pickle'))


# get percentiles model spectrum
###################################################################################################

model_spectra = {}
for SSC in SSCs:
    try:
        model_spectra[SSC['num']]
    except:
        model_spectra[SSC['num']] = {}

    for specie in data_XCLASS_refit[SSC['num']].keys():
        try:
            model_spectra[SSC['num']][specie]
        except:
            model_spectra[SSC['num']][specie] = {}

        for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
            try:
                model_spectra[SSC['num']][specie][component]
            except:
                model_spectra[SSC['num']][specie][component] = {}

            for band in ['LSB','USB']:
                try:
                    model_spectra[SSC['num']][specie][component][band]
                except:
                    model_spectra[SSC['num']][specie][component][band] = {}

for SSC in tqdm(SSCs):
    for band in ['LSB','USB']:
        for specie in data_XCLASS_refit[SSC['num']].keys():
            for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):

                frequency = intensity_spectra[SSC['num']][specie][component][num][band][0].value

                intensities = []
                for num in nums:
                    intensities.append( intensity_spectra[SSC['num']][specie][component][num][band][1].value )
                p16, median, p84 = np.percentile(intensities, (16,50,84), axis=0)

                model_spectra[SSC['num']][specie][component][band]['frequency'] = frequency
                model_spectra[SSC['num']][specie][component][band]['p16']       = p16
                model_spectra[SSC['num']][specie][component][band]['median']    = median
                model_spectra[SSC['num']][specie][component][band]['p84']       = p84

fnpickle(model_spectra, join(intensitydir, 'model_spectra.pickle'))


###################################################################################################
#
###################################################################################################
