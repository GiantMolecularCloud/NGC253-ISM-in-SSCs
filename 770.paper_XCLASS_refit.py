#######################
# GAS IN SSCS: XCLASS #
#######################

# address referee report: derive intensities from XCLASS fits

# refit a couple of CO lines because the results were inconsistent
# refit H13CN in some sources because the values were too low


###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '700.info.py'))
execfile(join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(join(mandir, 'spectra.pickle'))
detected_species = fnunpickle(join(XCLASSdir, 'detected_species.pickle'))

data_XCLASS = fnunpickle(join(resultsdir, 'data_XCLASS.pickle'))
data_Gauss  = fnunpickle(join(resultsdir, 'data_Gauss.pickle'))

# select one of the following setups and run script:

# # CO refit
# nums = np.arange(0,101)
# refitfiles   = {s: {str(num): {} for num in nums} for s in ['3','6','7','9','14']}
# refitspecies = 'CO;v=0'
# refitband    = 'LSB'

# # H13CN refit
# nums = np.arange(0,101)
# refitfiles   = {s: {str(num): {} for num in nums} for s in ['5','8','9','11','14']}
# refitspecies = 'HC-13-N;v=0'
# refitband    = 'LSB'

# # HCN v2=1 refit
# nums = np.arange(0,101)
# refitfiles   = {s: {str(num): {} for num in nums} for s in ['11']}
# refitspecies = 'HCN;v2=1'
# refitband    = 'USB'

# H13CN once again
nums = np.arange(0,101)
refitfiles   = {s: {str(num): {} for num in nums} for s in ['1','2','3','4','6','10','12','13']}
refitspecies = 'HC-13-N;v=0'
refitband    = 'LSB'


###################################################################################################
# get spectra
###################################################################################################

# copy over spectra from XCLASS_final run

for num in tqdm(nums):
    for SSCnum in refitfiles.keys():
        refitfiles[SSCnum][str(num)] = {'dat': {}}
        for band in ['LSB','USB']:
            final_spec = escape_fname(join(Xfinaldir,'SSC_'+SSCnum,'run_'+str(num),band+'.dat'))
            refit_spec  = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),band+'.dat'))
            os.system('mkdir -p '+os.path.dirname(refit_spec))
            os.system('cp '+final_spec+' '+refit_spec)
            refitfiles[SSCnum][str(num)]['dat'][band] = refit_spec


###################################################################################################
# get observation file
###################################################################################################

# copy over observation file from XCLASS_final run

for num in tqdm(nums):
    for SSCnum in refitfiles.keys():
        final_obs = escape_fname(join(Xfinaldir,'SSC_'+SSCnum,'run_'+str(num),'observation.xml'))
        refit_obs = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'observation.xml'))
        os.system('mkdir -p '+os.path.dirname(refit_obs))
        os.system('cp '+final_obs+' '+refit_obs)
        refitfiles[SSCnum][str(num)]['obs'] = refit_obs


###################################################################################################
# assemble final molfit file
###################################################################################################

# Keep other lines fixed and only refit the affected lines.
# Allow a small variation in linewidth to compensate line shape variations due to temperature effects.

def create_molfit_file(SSCnum):

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset"""

    for specie,data in data_XCLASS[SSCnum].items():

        # fix suddenly renamed SO species
        if specie=='SO;v=0;#1':
            specie = specie.replace('#1','#2')

        if ('#1' in specie) or ('#2' in specie):
            spx = specie
        else:
            spx = specie+';'

        ncomp = len(data['velocity']['median'])
        molfit += f"""
{spx}   {ncomp}"""

        for c in np.arange(ncomp):
            T = data['temperature']['median'][c]
            v = data['velocity']['median'][c]
            w = data['linewidth']['median'][c]
            N = data['column density']['median'][c]

            if specie==refitspecies:
                molfit += f"""
n  0.01  10  1   n  25  1000  {'{:.1f}'.format(T)}   y  {'{:.2E}'.format(N/100)}  {'{:.2E}'.format(N*100)}  {'{:.2E}'.format(N)}   y  {'{:.1f}'.format(w*0.8)}  {'{:.1f}'.format(w*1.2)}  {'{:.1f}'.format(w)}   y  {'{:.1f}'.format(v-10.)}  {'{:.1f}'.format(v+10.)}  {'{:.1f}'.format(v)}   c"""
            else:
                molfit += f"""
n  0.01  10  1   n  50  250  130   n  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}   n  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}   n  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}   c"""

    # save to disk
    savepath = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)

    # keep savepath
    refitfiles[SSCnum]['molfit'] = savepath


for SSCnum in tqdm(refitfiles.keys()):
    create_molfit_file(SSCnum)


print("##############################################################")
print("# manually modify molfit file to get the desired components! #")
print("##############################################################")


###################################################################################################
# create XCLASS script
###################################################################################################

def create_XCLASS_script(num, SSCnum):

    alg_file   = join(basescriptdir,'NGC253','paper_20a','algorithm.xml')
    obs_file   = refitfiles[SSCnum][str(num)]['obs']
    mol_file   = refitfiles[SSCnum]['molfit']
    resultpath = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'fit','results'))
    XCLASS = f"""import os
AlgorithmXMLFile = '{alg_file}'
experimentalData = '{obs_file}'
MolfitsFileName  = '{mol_file}'

newmolfit, modeldata, JobDir = myXCLASSFit()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')
"""

    # save to disk
    savepath = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'fit','XCLASS_fit.py'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(XCLASS)

    refitfiles[SSCnum][str(num)]['fit'] = savepath


for num in tqdm(nums):
    for SSCnum in refitfiles.keys():
        create_XCLASS_script(num, SSCnum)


###################################################################################################
# save settings
###################################################################################################

fnpickle(refitfiles, escape_fname(join(refitdir, refitspecies, 'refitfiles.pickle')))


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSCnum):
    XCL_files = [refitfiles[SSCnum][str(num)]['fit'] for num in nums]
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

commands = [CASA_command(SSCnum) for SSCnum in refitfiles.keys()]
pool = Pool(len(list(refitfiles.keys())))
pool.map(run_in_casa, commands)
pool.close()
pool.join()


execfile('/home/krieger/scripts/NGC253/paper_20a/771.referee_XCLASS_final_spectra.py')


###################################################################################################
#
###################################################################################################
