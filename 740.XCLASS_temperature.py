#######################
# GAS IN SSCS: XCLASS #
#######################

# Temperature fitting
# Keep other lines fixed and fit only for temperature sensitive lines.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))

data_N   = fnunpickle(os.path.join(Xfinaldir, 'data.pickle'))
data_I   = fnunpickle(os.path.join(mandir,'line_intensity_data.pickle'))
ratios_N = fnunpickle(os.path.join(Xfinaldir, 'ratios_column_density.pickle'))
ratios_I = fnunpickle(os.path.join(mandir, 'ratios_intensity.pickle'))


nums = np.arange(0,101)
tempfiles = {SSC['num']: {str(num): {} for num in nums} for SSC in SSCs}

temperature_species = ['H2CS;v=0;#1','SO2;v=0']
fnpickle(temperature_species, os.path.join(tempdir, 'temperature_species.pickle'))


###################################################################################################
# get spectra
###################################################################################################

# copy over spectra from XCLASS_final run

for num in tqdm(nums):
    for SSC in SSCs:
        tempfiles[SSC['num']][str(num)] = {'dat': {}}
        for band in ['LSB','USB']:
            final_spec = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),band+'.dat'))
            temp_spec  = escape_fname(os.path.join(tempdir,'SSC_'+str(SSC['no']),'run_'+str(num),band+'.dat'))
            os.system('mkdir -p '+os.path.dirname(temp_spec))
            os.system('cp '+final_spec+' '+temp_spec)
            tempfiles[SSC['num']][str(num)]['dat'][band] = temp_spec


###################################################################################################
# get observation file
###################################################################################################

# copy over observation file from XCLASS_final run

for num in tqdm(nums):
    for SSC in SSCs:
        final_obs = escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),'observation.xml'))
        temp_obs  = escape_fname(os.path.join(tempdir,'SSC_'+str(SSC['no']),'run_'+str(num),'observation.xml'))
        os.system('mkdir -p '+os.path.dirname(temp_obs))
        os.system('cp '+final_obs+' '+temp_obs)
        tempfiles[SSC['num']][str(num)]['obs'] = temp_obs


###################################################################################################
# assemble final molfit file
###################################################################################################

# Keep other lines fixed and only fit temperature sensitive lines for temperature and column density.
# Allow a small variation in linewidth to compensate line shape variations due to temperature effects.

def create_molfit_file(SSC):

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset"""

    for specie,data in data_N[SSC['num']].items():
        if '#1' in specie:
            spx = specie
        else:
            spx = specie+';'

        ncomp = len(data['velocity']['median'])
        molfit += f"""
{spx}   {ncomp}"""

        for c in np.arange(ncomp):
            v = data['velocity']['median'][c]
            w = data['linewidth']['median'][c]
            N = data['column density']['median'][c]

            if specie in temperature_species:
                molfit += f"""
n  0.01  10  1   y  25  1000  130   y  {'{:.2E}'.format(N/100)}  {'{:.2E}'.format(N*100)}  {'{:.2E}'.format(N)}   y  {'{:.1f}'.format(w/0.8)}  {'{:.1f}'.format(w*1.2)}  {'{:.1f}'.format(w)}   n  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}   c"""
            else:
                molfit += f"""
n  0.01  10  1   n  50  250  130   n  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}   n  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}   n  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}   c"""

    # save to disk
    savepath = escape_fname(os.path.join(tempdir,'SSC_'+str(SSC['no']),'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)

    # keep savepath
    tempfiles[SSC['num']]['molfit'] = savepath


for SSC in tqdm(SSCs):
    create_molfit_file(SSC)


###################################################################################################
# create XCLASS script
###################################################################################################

def create_XCLASS_script(num, SSC):

    alg_file   = os.path.join(basescriptdir,'NGC253','paper_20a','algorithm.xml')
    obs_file   = tempfiles[SSC['num']][str(num)]['obs']
    mol_file   = tempfiles[SSC['num']]['molfit']
    resultpath = escape_fname(os.path.join(tempdir,'SSC_'+str(SSC['no']),'run_'+str(num),'fit','results'))
    XCLASS = f"""import os
AlgorithmXMLFile = '{alg_file}'
experimentalData = '{obs_file}'
MolfitsFileName  = '{mol_file}'

newmolfit, modeldata, JobDir = myXCLASSFit()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')
"""

    # save to disk
    savepath = escape_fname(os.path.join(tempdir,'SSC_'+str(SSC['no']),'run_'+str(num),'fit','XCLASS_fit.py'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(XCLASS)

    tempfiles[SSC['num']][str(num)]['fit'] = savepath


for num in tqdm(nums):
    for SSC in SSCs:
        create_XCLASS_script(num, SSC)


###################################################################################################
# save settings
###################################################################################################

fnpickle(tempfiles, os.path.join(tempdir, 'tempfiles.pickle'))
fnpickle(nums, os.path.join(tempdir, 'nums.pickle'))


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSC):
    XCL_files = [tempfiles[SSC['num']][str(num)]['fit'] for num in nums]
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
