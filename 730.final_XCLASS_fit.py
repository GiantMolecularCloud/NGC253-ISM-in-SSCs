#######################
# GAS IN SSCS: XCLASS #
#######################

# Prepare the data in python3 in a way that keeps the XCLASS code simple and manageable.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
fitable_species = fnunpickle(os.path.join(XCLASSdir, 'fitable_species.pickle'))

nums = np.arange(0,101)
temperature = 130*u.K

Xfiles = {SSC['num']: {str(num): {} for num in nums} for SSC in SSCs}

execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))


###################################################################################################
# get spectra
###################################################################################################

def export_spectrum(num, SSC, rms):
    try:
        Xfiles[SSC['num']][str(num)]
    except:
        Xfiles[SSC['num']][str(num)] = {}
    Xfiles[SSC['num']][str(num)]['dat'] = {}
    for band in ['LSB','USB']:

        # get data
        spectrum  = spectra[str(SSC['no'])][band]
        frequency = spectrum['frequency'].to(u.MHz)
        intensity = spectrum['spectrum'].to(u.K)

        # shift spectrum to rest frequency
        velshift  = SSC['velshift']
        frequency = [(-vsys-velshift).to(u.GHz, equivalencies=u.doppler_optical(f)).value for f in frequency]*u.GHz

        # remove NaNs
        frequency, intensity = crossmatch(frequency.to(u.MHz).value, intensity.to(u.K).value)

        # get limits
        fmin = frequency[0]
        fmax = frequency[-1]

        # add noise
        if not num==0:
            randstate = np.random.RandomState(num)
            noise =  np.random.normal(loc=0., scale=rms.to(u.K).value, size=len(frequency))
            intensity += noise

        # export to ASCII file
        savepath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),band+'.dat'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        np.savetxt(savepath, np.transpose([frequency, intensity]), fmt=('%10.4f', '%10.4f'))

        # keep savepath and fmin/fmax
        Xfiles[SSC['num']][str(num)]['dat'][band] = {'file': savepath, 'fmin': fmin, 'fmax': fmax}


for num in tqdm(nums):
    for SSC in SSCs:
        export_spectrum(num, SSC, rms=0.46*u.K)


###################################################################################################
# get observation file
###################################################################################################

def create_observations_file(num, SSC, fstep=2.5*u.MHz):

    observation = f"""<?xml version="1.0" encoding="UTF-8"?>
<ExpFiles>
<NumberExpFiles>2</NumberExpFiles>"""

    for band,dat in Xfiles[SSC['num']][str(num)]['dat'].items():
        observation += f"""
<file>
    <FileNamesExpFiles>{dat['file']}</FileNamesExpFiles>
    <ImportFilter>xclassASCII</ImportFilter>
    <NumberExpRanges>1</NumberExpRanges>
    <FrequencyRange>
        <MinExpRange>{dat['fmin']}</MinExpRange>
        <MaxExpRange>{dat['fmax']}</MaxExpRange>
        <StepFrequency>{fstep.to(u.MHz).value}</StepFrequency>
        <t_back_flag>False</t_back_flag>
        <BackgroundTemperature>0.0</BackgroundTemperature>
        <TemperatureSlope>0.0</TemperatureSlope>
        <HydrogenColumnDensity>0.e+0</HydrogenColumnDensity>
        <DustBeta>0.0</DustBeta>
        <Kappa>0.0</Kappa>
    </FrequencyRange>
    <GlobalvLSR>0.0</GlobalvLSR>
    <TelescopeSize>0.15</TelescopeSize>
    <Inter_Flag>True</Inter_Flag>
    <ErrorY>no</ErrorY>
    <NumberHeaderLines>0</NumberHeaderLines>
    <SeparatorColumns> </SeparatorColumns>
</file>"""

    observation += f"""
<iso_flag>False</iso_flag>
</ExpFiles>"""

    # save to disk
    savepath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'observation.xml'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(observation)

    # keep savepath
    Xfiles[SSC['num']][str(num)]['obs'] = {'file': savepath}


for num in tqdm(nums):
    for SSC in SSCs:
        create_observations_file(num, SSC)


###################################################################################################
# assemble final molfit file
###################################################################################################

# parse model molfit files
# build a complete molfit file for joint fitting of the full band

# load independent fits
detected_species = fnunpickle(os.path.join(XCLASSdir, 'detected_species.pickle'))
line_data_N = fnunpickle(os.path.join(XCLASSdir, 'line_column_density_data.pickle'))

def get_molfit_params(SSC):

    def get_molfit(SSC,spx):
        molfit = open(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),spx,'molecules_final.molfit')), 'r').readlines()[2:]
        molfit = molfit_list(molfit,spx)
        return molfit

    def molfit_list(molfit, spx):
        for idx,ll in enumerate(molfit):
            ll = ll.replace('\n','')
            ll = [l for l in ll.split(' ') if not l=='']
            molfit[idx] = ll
        return molfit

    def get_specie_index(molfit,spx):
        if not '#1' in spx:
            spx += ';'
        for idx,ll in enumerate(molfit):
            if ll[0]==spx:
                idx_start = idx+1
            if ll[0]!='n' and idx>1:
                idx_stop  = idx
                break
            if idx==len(molfit)-1:
                idx_stop  = idx+1
        return idx_start,idx_stop

    def get_values(molfit,spx,comp):
        t_llim = float(molfit[comp][5])
        t_ulim = float(molfit[comp][6])
        t_val  = float(molfit[comp][7])
        N_llim = float(molfit[comp][9])
        N_ulim = float(molfit[comp][10])
        N_val  = float(molfit[comp][11])
        w_llim = float(molfit[comp][13])
        w_ulim = float(molfit[comp][14])
        w_val  = float(molfit[comp][15])
        v_llim = float(molfit[comp][17])
        v_ulim = float(molfit[comp][18])
        v_val  = float(molfit[comp][19])
        params = {'temperature':    {'min': t_llim, 'med': t_val, 'max': t_ulim},
                  'velocity':       {'min': v_llim, 'med': v_val, 'max': v_ulim},
                  'linewidth':      {'min': w_llim, 'med': w_val, 'max': w_ulim},
                  'column density': {'min': N_llim, 'med': N_val, 'max': N_ulim}}
        return params


    Xfiles[SSC['num']]['molfit'] = {}
    for spx in detected_species[SSC['num']]:
        molfit = get_molfit(SSC,spx)
        idx_start,idx_stop = get_specie_index(molfit,spx)
        n_comps = idx_stop-idx_start
        Xfiles[SSC['num']]['molfit'][spx] = {'n': n_comps, 'components': []}

        for comp in np.arange(idx_start,idx_stop):
            params = get_values(molfit,spx,comp)
            Xfiles[SSC['num']]['molfit'][spx]['components'].append(params)


def create_molfit_file(SSC):

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset"""

    for specie in detected_species[SSC['num']]:
        if '#1' in specie:
            spx = specie
        else:
            spx = specie+';'

        molfit += f"""
{spx}   {Xfiles[SSC['num']]['molfit'][specie]['n']}"""

        for component in Xfiles[SSC['num']]['molfit'][specie]['components']:
            v = component['velocity']
            w = component['linewidth']
            N = component['column density']

            molfit += f"""
n  0.01  10  1   n  50  250  {temperature.value}   y  {'{:.2E}'.format(N['min'])}  {'{:.2E}'.format(N['max'])}  {'{:.2E}'.format(N['med'])}   y  {'{:.1f}'.format(w['min'])}  {'{:.1f}'.format(w['max'])}  {'{:.1f}'.format(w['med'])}   y  {'{:.1f}'.format(v['min'])}  {'{:.1f}'.format(v['max'])}  {'{:.1f}'.format(v['med'])}   c"""

    # save to disk
    savepath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)

    # keep savepath
    Xfiles[SSC['num']]['molfit']['file'] = savepath


for SSC in tqdm(SSCs):
    get_molfit_params(SSC)

for SSC in tqdm(SSCs):
    create_molfit_file(num, SSC)


###################################################################################################
# create XCLASS script
###################################################################################################

def create_XCLASS_script(num, SSC):

    alg_file   = os.path.join(basescriptdir,'NGC253','paper_20a','algorithm.xml')
    obs_file   = Xfiles[SSC['num']][str(num)]['obs']['file']
    mol_file   = Xfiles[SSC['num']]['molfit']['file']
    resultpath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'fit','results'))
    XCLASS = f"""import os
AlgorithmXMLFile = '{alg_file}'
experimentalData = '{obs_file}'
MolfitsFileName  = '{mol_file}'

newmolfit, modeldata, JobDir = myXCLASSFit()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')
"""

    # save to disk
    savepath = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'fit','XCLASS_fit.py'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(XCLASS)

    Xfiles[SSC['num']][str(num)]['fit'] = {'file': savepath}


for num in tqdm(nums):
    for SSC in SSCs:
        create_XCLASS_script(num, SSC)


###################################################################################################
# save settings
###################################################################################################

fnpickle(Xfiles, os.path.join(Xfinaldir, 'Xfiles.pickle'))
fnpickle(temperature, os.path.join(Xfinaldir, 'temperature.pickle'))
fnpickle(nums, os.path.join(Xfinaldir, 'nums.pickle'))


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSC):
    XCL_files = [Xfiles[SSC['num']][str(num)]['fit']['file'] for num in nums]
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


###################################################################################################
#
###################################################################################################
